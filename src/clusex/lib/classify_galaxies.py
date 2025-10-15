#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import base64
import csv
import json
import mimetypes
import os
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd
from tqdm import tqdm
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type


def ensure_api_key():
    key = os.environ.get("OPENAI_API_KEY")
    if not key:
        raise RuntimeError(
            "OPENAI_API_KEY it is not defined\n"
            "Export key before continue:\n"
            '  export OPENAI_API_KEY="sk-..."'
        )



SCHEMA = {
    "type": "json_schema",
    "json_schema": {
        "name": "galaxy_morphology",
        "strict": True,
        "schema": {
            "type": "object",
            "properties": {
                "class": {
                    "type": "string",
                    "description": "Main morphological class",
                    "enum": ["E", "S0", "Sa", "Sb", "Sc", "Irr", "Uncertain"]
                },
                "bar": {"type": "boolean", "description": "Bar present?"},
                "ring": {"type": "boolean", "description": "Ring present?"},
                "inclination_deg": {
                    "type": "number",
                    "description": "Estimated disk inclination in degrees (0=face-on, 90=edge-on)"
                },
                "confidence": {
                    "type": "number",
                    "minimum": 0.0,
                    "maximum": 1.0,
                    "description": "Global confidence (0–1)"
                },
                "notes": {
                    "type": "string",
                    "description": "Brief rationale (e.g., low S/N, edge-on dust lane, interacting)"
                }
            },
            "required": ["class", "bar", "ring", "inclination_deg", "confidence", "notes"],
            "additionalProperties": False
        }
    }
}


def load_paths(txt_file: str) -> List[str]:
    with open(txt_file, "r") as f:
        return [line.strip() for line in f if line.strip()]


def file_to_data_url(path: str) -> str:
    mime, _ = mimetypes.guess_type(path)
    if mime is None:
        ext = Path(path).suffix.lower()
        if ext in [".jpg", ".jpeg"]:
            mime = "image/jpeg"
        elif ext == ".png":
            mime = "image/png"
        elif ext == ".webp":
            mime = "image/webp"
        elif ext == ".gif":
            mime = "image/gif"
        else:
            raise ValueError(f"No se reconoce el tipo MIME de: {path}")
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("utf-8")
    return f"data:{mime};base64,{b64}"


def chunked(lst: List[Any], n: int):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]



def build_messages_for_batch(image_paths: List[str]) -> List[Dict[str, Any]]:
    content = []
    for p in image_paths:
        data_url = file_to_data_url(p)
        content.append({"type": "input_image", "image_url": data_url})

    instruction = (
        "Task: morphological classification of each galaxy image.\n"
        "Output requirements (CRITICAL):\n"
        "• Return ONLY a strict JSON list with one object per image, in the same order.\n"
        "• Each object MUST include ALL fields: class, bar, ring, inclination_deg, confidence, notes.\n"
        "• If information is insufficient or ambiguous, use class='Uncertain', set bar=false, ring=false, "
        "  set confidence ≤ 0.4, and explain briefly in notes (e.g., 'low S/N', 'very small', 'saturated').\n"
        "• Do NOT leave any field empty. Do NOT output extra text outside the JSON.\n"
        "• Classes allowed: E, S0, Sa, Sb, Sc, Irr, Uncertain.\n"
        "• inclination_deg: numeric in [0,90].\n"
        "• Keep answers concise.\n"
    )
    content.append({"type": "input_text", "text": instruction})
    return [{"role": "user", "content": content}]


def ensure_list_of_dicts(obj: Any, expected_len: int) -> List[Dict[str, Any]]:
    if not isinstance(obj, list):
        raise ValueError("La salida no es una lista JSON.")
    if len(obj) != expected_len:
        raise ValueError(f"Tamaño de salida inesperado: {len(obj)} != {expected_len}")
    for i, it in enumerate(obj):
        if not isinstance(it, dict):
            raise ValueError(f"Elemento {i} no es un objeto JSON.")
    return obj


# --- Cliente tardío para permitir tests sin importar si falta la key.
def _get_client():
    from openai import OpenAI
    return OpenAI()


@retry(
    reraise=True,
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=1, max=20),
    retry=retry_if_exception_type(Exception)
)
def call_api(client, model: str, messages: List[Dict[str, Any]], response_format: Dict[str, Any]) -> Dict[str, Any]:
    """
    Estrategia Opción A:
    1) Intenta con response_format (salida validada).
    2) Si el SDK no acepta ese argumento (TypeError), reintenta sin response_format.
       Luego extrae el primer bloque JSON del texto y lo parsea.
    """
    import re

    # Intento moderno
    try:
        resp = client.responses.create(
            model=model,
            input=messages,
            response_format=response_format,
            temperature=0.2,
        )
        return json.loads(resp.output_text)

    except TypeError as e:
        if "unexpected keyword argument 'response_format'" not in str(e):
            raise  # otro TypeError

        # Fallback sin response_format
        resp = client.responses.create(
            model=model,
            input=messages,
            temperature=0.2,
        )
        raw = resp.output_text or ""
        # Busca el primer objeto o arreglo JSON.
        m = re.search(r'(\{.*\}|\[.*\])', raw, re.S)
        if not m:
            raise RuntimeError(
                "No se encontró JSON en la salida del modelo (fallback sin response_format). "
                f"Primeros 200 chars: {raw[:200]!r}"
            )
        return json.loads(m.group(1))


def classify():
    parser = argparse.ArgumentParser(description="Clasificación morfológica de galaxias (Opción A: fallback sin response_format).")
    parser.add_argument("--paths-file", required=True, help="Archivo de texto con rutas de imágenes, una por línea.")
    parser.add_argument("--out", required=True, help="Ruta del CSV de salida.")
    parser.add_argument("--batch-size", type=int, default=8, help="Imágenes por solicitud (≤10 recomendado).")
    parser.add_argument("--model", default="gpt-4.1-mini", help="Modelo de visión, p.ej. gpt-4.1-mini.")
    args = parser.parse_args()

    ensure_api_key()
    client = _get_client()

    image_paths = load_paths(args.paths_file)
    if not image_paths:
        print("No se encontraron rutas en el archivo de paths.")
        return

    rows = []
    for batch_paths in tqdm(list(chunked(image_paths, args.batch_size)), desc="Procesando lotes"):
        # Construcción del mensaje
        try:
            messages = build_messages_for_batch(batch_paths)
        except Exception as e:
            for p in batch_paths:
                rows.append({
                    "path": p, "class": "", "bar": "", "ring": "",
                    "inclination_deg": "", "confidence": "", "notes": f"ERROR EN LECTURA/ENCODE: {repr(e)}"
                })
            continue

        # Llamada a la API (con fallback)
        try:
            result = call_api(client, args.model, messages, SCHEMA)
        except Exception as e:
            for p in batch_paths:
                rows.append({
                    "path": p, "class": "", "bar": "", "ring": "",
                    "inclination_deg": "", "confidence": "", "notes": f"ERROR API: {repr(e)}"
                })
            continue

        # Validación de formato y mapeo
        try:
            result = ensure_list_of_dicts(result, expected_len=len(batch_paths))
              
            for p, r in zip(batch_paths, result):
                rows.append({
                    "path": p,
                    "class": r.get("class", "Uncertain") or "Uncertain",
                    "bar": r.get("bar", False),
                    "ring": r.get("ring", False),
                    "inclination_deg": r.get("inclination_deg", 0),
                    "confidence": r.get("confidence", 0.2),
                    "notes": r.get("notes", "Model returned incomplete fields")
                })

        except Exception as e:
            for p in batch_paths:
                rows.append({
                    "path": p, "class": "", "bar": "", "ring": "",
                    "inclination_deg": "", "confidence": "", "notes": f"FORMATO INVALIDO: {repr(e)}"
                })

    # Escritura CSV
    out_cols = ["path", "class", "bar", "ring", "inclination_deg", "confidence", "notes"]
    Path(os.path.dirname(args.out) or ".").mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    # Parquet opcional
    try:
        df = pd.DataFrame(rows, columns=out_cols)
        df.to_parquet(Path(args.out).with_suffix(".parquet"))
    except Exception:
        pass

    print(f"Done. Results en: {args.out}")


if __name__ == "__main__":
    classify()
