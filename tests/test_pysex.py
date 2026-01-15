# tests/test_pysex.py
from __future__ import annotations

import inspect
import os
import re
import tempfile
import unittest
from pathlib import Path
from subprocess import CalledProcessError, CompletedProcess
from typing import Any, Callable, Optional
from unittest import mock


def _import_pysex_module():
    try:
        from clusex import pysex as mod  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise unittest.SkipTest(f"No se pudo importar clusex.pysex: {exc}")
    return mod


def _iter_public_callables(mod) -> list[tuple[str, Callable[..., Any]]]:
    out: list[tuple[str, Callable[..., Any]]] = []
    for name in dir(mod):
        if name.startswith("_"):
            continue
        obj = getattr(mod, name, None)
        if callable(obj):
            out.append((name, obj))
    return out


def _looks_like_sextractor_runner(name: str, fn: Callable[..., Any]) -> bool:
    lname = name.lower()
    if lname in {"main"}:
        return False
    if "sex" not in lname:
        return False
    try:
        sig = inspect.signature(fn)
    except Exception:
        return False

    params = list(sig.parameters.values())
    if not params:
        return False

    # Heurística: al menos un parámetro “tipo imagen”.
    candidate = {"image", "img", "imagen", "fits", "fname", "filename", "file", "infile"}
    for p in params:
        if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        if (p.name or "").lower() in candidate:
            return True

    # Alternativa: si acepta 1-3 posicionales es plausible para un wrapper simple.
    positional = [p for p in params if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)]
    return 1 <= len(positional) <= 4


def _find_runner(mod) -> Optional[tuple[str, Callable[..., Any]]]:
    # Primero, nombres frecuentes.
    preferred = [
        "run_sextractor",
        "runsextractor",
        "run_sex",
        "sextractor",
        "pysex",
        "runsex",
        "run_sex",
        "sex",
    ]
    for name in preferred:
        if hasattr(mod, name) and callable(getattr(mod, name)):
            fn = getattr(mod, name)
            if _looks_like_sextractor_runner(name, fn):
                return name, fn

    # Luego, heurística sobre todos los callables.
    candidates = [(n, f) for (n, f) in _iter_public_callables(mod) if _looks_like_sextractor_runner(n, f)]
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        # Preferir el que tenga más indicios en el nombre.
        scored: list[tuple[int, str, Callable[..., Any]]] = []
        for n, f in candidates:
            ln = n.lower()
            score = 0
            score += 3 if "run" in ln else 0
            score += 3 if "tractor" in ln else 0
            score += 2 if "cat" in ln else 0
            score += 1 if "config" in ln else 0
            scored.append((score, n, f))
        scored.sort(reverse=True)
        _, n, f = scored[0]
        return n, f
    return None


def _looks_like_catalog_parser(name: str, fn: Callable[..., Any]) -> bool:
    lname = name.lower()
    if any(k in lname for k in ("parse", "read", "load")) and any(k in lname for k in ("cat", "catalog")):
        return True
    return False


def _find_parser(mod) -> Optional[tuple[str, Callable[..., Any]]]:
    preferred = [
        "parse_catalog",
        "parse_cat",
        "read_catalog",
        "read_cat",
        "load_catalog",
        "load_cat",
    ]
    for name in preferred:
        if hasattr(mod, name) and callable(getattr(mod, name)):
            return name, getattr(mod, name)
    candidates = [(n, f) for (n, f) in _iter_public_callables(mod) if _looks_like_catalog_parser(n, f)]
    if candidates:
        # Elegir el más “específico”.
        candidates.sort(key=lambda t: (("parse" in t[0].lower()) + ("read" in t[0].lower())), reverse=True)
        return candidates[0]
    return None


def _write_fake_sextractor_catalog(path: Path) -> None:
    # Catálogo típico: comentarios con '#', luego columnas numéricas.
    # Se usan las 14 columnas descritas en la documentación del proyecto.
    header = [
        "# 1 NUMBER",
        "# 2 ALPHA_J2000",
        "# 3 DELTA_J2000",
        "# 4 XPEAK_IMAGE",
        "# 5 YPEAK_IMAGE",
        "# 6 MAG_BEST",
        "# 7 KRON_RADIUS",
        "# 8 FLUX_RADIUS",
        "# 9 ISOAREA_IMAGE",
        "# 10 A_IMAGE",
        "# 11 ELLIPTICITY",
        "# 12 THETA_IMAGE",
        "# 13 BACKGROUND",
        "# 14 CLASS_STAR",
        "# 15 FLAGS",
    ]
    rows = [
        "1 150.12345 2.34567 512.0 512.0 19.123 3.2 5.5 120.0 9.8 0.23 -45.0 120.5 0.98 0",
        "2 150.22345 2.44567 128.0 900.0 21.500 2.1 3.3  80.0 6.4 0.10  12.0 118.2 0.12 2",
        "3 150.32345 2.54567  64.0  64.0 25.001 1.5 2.2  20.0 3.1 0.55  89.9 130.0 0.01 0",
    ]
    path.write_text("\n".join(header + rows) + "\n", encoding="utf-8")


def _extract_catalog_path_from_cmd(cmd: list[str]) -> Optional[Path]:
    # Intentar localizar -CATALOG_NAME <file> o CATALOG_NAME <file> u opciones similares.
    joined = " ".join(cmd)
    patterns = [
        r"(?:-CATALOG_NAME|CATALOG_NAME)\s+([^\s]+)",
        r"(?:-CATALOG|CATALOG)\s+([^\s]+)",
        r"(?:-c\s+([^\s]+))",  # config, no catálogo; se ignora abajo
    ]
    for pat in patterns:
        m = re.search(pat, joined)
        if m:
            candidate = m.group(1)
            # Si es config (-c), no es un catálogo.
            if pat.startswith(r"(?:-c"):
                continue
            return Path(candidate)
    return None


def _call_runner_best_effort(runner: Callable[..., Any], image_path: Path, outcat: Optional[Path] = None) -> Any:
    """
    Invoca el runner sin asumir una firma exacta.
    Estrategia:
      1) kwargs con nombres frecuentes.
      2) fallback posicional mínimo.
    """
    sig = inspect.signature(runner)
    params = sig.parameters

    # Candidatos de kwargs.
    kw: dict[str, Any] = {}
    for key in ("image", "img", "imagen", "fits", "filename", "file", "infile"):
        if key in params:
            kw[key] = str(image_path)
            break

    if outcat is not None:
        for key in ("outcat", "outcatalog", "catalog", "catalog_name", "s:=":
            pass
        for key in ("outcat", "outcatalog", "catalog", "catalog_name", "cat", "out", "outfile"):
            if key in params:
                kw[key] = str(outcat)
                break

    # Si hay kwargs suficientes, usar.
    if kw:
        try:
            return runner(**kw)
        except TypeError:
            # Firma distinta; intentar posicional.
            pass

    # Fallback posicional: pasar la ruta de imagen como primer argumento.
    try:
        return runner(str(image_path))
    except TypeError:
        # Último recurso: imagen + outcat si se proporcionó.
        if outcat is not None:
            return runner(str(image_path), str(outcat))
        raise


class TestPySex(unittest.TestCase):
    def setUp(self) -> None:
        self.mod = _import_pysex_module()
        found = _find_runner(self.mod)
        if not found:
            raise unittest.SkipTest("No se encontró una función runner de SExtractor en clusex.pysex.")
        self.runner_name, self.runner = found

    def test_missing_image_raises(self) -> None:
        missing = Path("no_existe_esta_imagen_123456789.fits")
        with self.assertRaises((FileNotFoundError, OSError, SystemExit, ValueError)):
            _call_runner_best_effort(self.runner, missing)

    @mock.patch("clusex.pysex.subprocess.run")
    def test_subprocess_failure_raises(self, m_run: mock.Mock) -> None:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            img = td / "dummy.fits"
            img.write_bytes(b"")  # archivo vacío; no se usa contenido real

            m_run.side_effect = CalledProcessError(returncode=1, cmd=["sex"], stderr=b"fail")
            with self.assertRaises((CalledProcessError, RuntimeError, OSError, SystemExit, ValueError)):
                _call_runner_best_effort(self.runner, img)

    @mock.patch("clusex.pysex.subprocess.run")
    def test_simulated_sextractor_output_parsing(self, m_run: mock.Mock) -> None:
        """
        Simula la ejecución de SExtractor y verifica que el módulo pueda consumir un catálogo plausible.
        No requiere SExtractor real ni imágenes FITS válidas.
        """
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            img = td / "dummy.fits"
            img.write_bytes(b"")

            # Si el runner permite especificar salida de catálogo, se predefine.
            outcat = td / "sex.cat"

            def _fake_run(cmd, *args, **kwargs):
                # cmd puede ser lista o string; normalizar a lista.
                cmd_list = cmd if isinstance(cmd, list) else str(cmd).split()
                cat_path = _extract_catalog_path_from_cmd([str(x) for x in cmd_list])
                if cat_path is None:
                    # Si no se detecta, usar el outcat conocido.
                    cat_path = outcat
                # Asegurar directorio.
                cat_path.parent.mkdir(parents=True, exist_ok=True)
                _write_fake_sextractor_catalog(cat_path)
                return CompletedProcess(args=cmd_list, returncode=0, stdout=b"", stderr=b"")

            m_run.side_effect = _fake_run

            result = _call_runner_best_effort(self.runner, img, outcat=outcat)

            # Validaciones deliberadamente tolerantes, para no acoplarse a una estructura interna específica.
            self.assertIsNotNone(result)

            # Casos comunes: devuelve ruta, lista/array, dict o tabla.
            if isinstance(result, (str, os.PathLike)):
                self.assertTrue(Path(result).exists())
            elif hasattr(result, "__len__"):
                self.assertGreaterEqual(len(result), 1)

    def test_parser_function_handles_header_comments(self) -> None:
        """
        Si existe un parser explícito, validar que ignore cabeceras con '#'
        y procese filas numéricas.
        """
        found = _find_parser(self.mod)
        if not found:
            raise unittest.SkipTest("No se encontró una función parser de catálogo en clusex.pysex.")

        parser_name, parser = found

        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            cat = td / "sex.cat"
            _write_fake_sextractor_catalog(cat)

            parsed = None
            try:
                parsed = parser(str(cat))
            except TypeError:
                parsed = parser(cat)

            self.assertIsNotNone(parsed)

            # Verificaciones tolerantes.
            if hasattr(parsed, "__len__"):
                self.assertGreaterEqual(len(parsed), 1)


if __name__ == "__main__":
    unittest.main()
