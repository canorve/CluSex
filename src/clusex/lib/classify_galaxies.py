#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import base64
import csv
import json
import mimetypes
import os
import re
import shutil
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd
from tqdm import tqdm
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type


def ensure_api_key():
    key = os.environ.get("OPENAI_API_KEY")
    if not key:
        raise RuntimeError(
            "OPENAI_API_KEY is not defined.\n"
            "Export the key before run:\n"
            '  export OPENAI_API_KEY="sk-..."'
        )


def load_textfile(path: str) -> str:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"prompt file was not found: {path}")
    txt = p.read_text(encoding="utf-8").strip()
    if not txt:
        raise ValueError("prompt file is empty")
    return txt


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
                    "description": "Morphologic class",
                    "enum": ["E","S0","S0-","S0/a","Sa","Sab","Sb","Sbc",
                             "Sc","Scd","Sd","Sdm","Sm","Im","Unknown",
                             "S","Irr","Uncertain"]
                },
                "bar": {"type": "boolean"},
                "ring": {"type": "boolean"},
                "inclination_deg": {"type": "number"},
                "confidence": {"type": "number", "minimum": 0.0, "maximum": 1.0},
                "notes": {"type": "string"}
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
            raise ValueError(f"MIME type is not recognized: {path}")
    b = Path(path).read_bytes()
    import base64 as _b64
    return f"data:{mime};base64,{_b64.b64encode(b).decode('utf-8')}"


def chunked(lst: List[Any], n: int):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def build_messages_for_batch(image_paths: List[str], instruction_text: str) -> List[Dict[str, Any]]:
    content = []
    for p in image_paths:
        data_url = file_to_data_url(p)
        content.append({"type": "input_image", "image_url": data_url})
    content.append({"type": "input_text", "text": instruction_text})
    return [{"role": "user", "content": content}]


def ensure_list_of_dicts(obj: Any, expected_len: int) -> List[Dict[str, Any]]:
    if not isinstance(obj, list):
        raise ValueError("output is not a JSON list.")
    if len(obj) != expected_len:
        raise ValueError(f"output size unexpected: {len(obj)} != {expected_len}")
    for i, it in enumerate(obj):
        if not isinstance(it, dict):
            raise ValueError(f"Element {i} is not an JSON object.")
    return obj


def _get_client():
    from openai import OpenAI
    return OpenAI()


@retry(
    reraise=True,
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=1, max=20),
    retry=retry_if_exception_type(Exception)
)
def call_api(client, model: str, messages: List[Dict[str, Any]], response_format: Dict[str, Any], temp: float) -> Dict[str, Any]:
    """
    Option:
      1) Try with response_format (esqueme validation).
      2) If SDK does not admit these argument, fallback without response_format + parsing JSON.
    """
    try:
        resp = client.responses.create(
            model=model,
            input=messages,
            response_format=response_format,
            temperature=temp,
        )
        return json.loads(resp.output_text)

    except TypeError as e:
        if "unexpected keyword argument 'response_format'" not in str(e):
            raise
        # Fallback sin response_format
        resp = client.responses.create(
            model=model,
            input=messages,
            temperature=temp,
        )
        raw = resp.output_text or ""
        m = re.search(r'(\{.*\}|\[.*\])', raw, re.S)
        if not m:
            raise RuntimeError(
                "JSON was not found in output (fallback without response_format). "
                f"First 200 chars: {raw[:200]!r}"
            )
        return json.loads(m.group(1))


# -------- Copia/renombrado --------

def sanitize_label(label: str) -> str:
    """Convert the class in a fragment for file names"""
    lab = re.sub(r'[^A-Za-z0-9]+', '-', str(label).strip())
    return lab.strip('-') or "Unknown"


def build_name_suffix(cls: str, bar: bool, include_bar: bool) -> str:
    """
    Returns sufix for the name: _<CLASS> o _<CLASS>_<BAR|NOBAR>.
    """
    base = sanitize_label(cls)
    if include_bar:
        return f"{base}_{'BAR' if bar else 'NOBAR'}"
    return base


def copy_with_label(src_path: str, out_dir: Path, cls: str, bar: bool, include_bar: bool) -> Path:
    """
    Copy the image to out_dir with suffix _<CLASS>[_BAR|_NOBAR] before the file extension.
    Avoid overwriting.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    p = Path(src_path)
    stem = p.stem
    ext = p.suffix
    suffix = build_name_suffix(cls, bar, include_bar)
    target = out_dir / f"{stem}_{suffix}{ext}"
    idx = 1
    while target.exists():
        target = out_dir / f"{stem}_{suffix}_{idx}{ext}"
        idx += 1
    shutil.copy2(p, target)
    return target


# -------- main function --------

def classify():
    parser = argparse.ArgumentParser(
        description="Morphologic classification (extern prompt, fallback without response_format) + copy/rename optional."
    )
    parser.add_argument("--paths-file", required=True, help="text file with images paths, one per line")
    parser.add_argument("--prompt-file", required=True, help="File .txt with the prompt instruction.")
    parser.add_argument("--out", required=True, help="CSV file output.")
    parser.add_argument("--batch-size", type=int, default=8, help="Images per batch (≤10 recomended).")
    parser.add_argument("--temperature", type=float, default=0.2, help="temperature. Less temperature less variation (default=0.2).")
    parser.add_argument("--model", default="gpt-4.1-mini", help="Vision model, example: gpt-4.1-mini.")
    parser.add_argument("--copy-dir", default=None, help="If activated, it will copy each image renamed to this directory.")
    parser.add_argument("--bar-in-name", action="store_true",
                        help="If it is used with --copy-dir, adds _BAR or _NOBAR to the name of the copy.")
    args = parser.parse_args()

    ensure_api_key()
    instruction_text = load_textfile(args.prompt_file)
    client = _get_client()

    image_paths = load_paths(args.paths_file)
    if not image_paths:
        print("Paths were not found in the file.")
        return

    
    copy_dir = Path(args.copy_dir) if args.copy_dir else None   # <-- CORRECTA
    # Corrección del guion en acceso:
    copy_dir = Path(args.copy_dir) if args.copy_dir else None

    rows = []
    for batch_paths in tqdm(list(chunked(image_paths, args.batch_size)), desc="Procesing batch"):
        # Construcción del mensaje
        try:
            messages = build_messages_for_batch(batch_paths, instruction_text)
        except Exception as e:
            for p in batch_paths:
                row = {
                    "path": p, "class": "Unknown", "bar": False, "ring": False,
                    "inclination_deg": 0, "confidence": 0.2, "notes": f"ERROR IN LECTURE/ENCODE: {repr(e)}"
                }
                rows.append(row)
                if copy_dir:
                    try:
                        copy_with_label(p, copy_dir, row["class"], row["bar"], args.bar_in_name)
                    except Exception:
                        pass
            continue

        # Llamada a la API (con fallback)
        try:
            result = call_api(client, args.model, messages, SCHEMA, args.temperature)
        except Exception as e:
            for p in batch_paths:
                row = {
                    "path": p, "class": "Unknown", "bar": False, "ring": False,
                    "inclination_deg": 0, "confidence": 0.2, "notes": f"ERROR API: {repr(e)}"
                }
                rows.append(row)
                if copy_dir:
                    try:
                        copy_with_label(p, copy_dir, row["class"], row["bar"], args.bar_in_name)
                    except Exception:
                        pass
            continue

        # Validación de formato, mapeo y copia
        try:
            result = ensure_list_of_dicts(result, expected_len=len(batch_paths))
            for p, r in zip(batch_paths, result):
                row = {
                    "path": p,
                    "class": (r.get("class") or "Unknown"),
                    "bar": bool(r.get("bar", False)),
                    "ring": bool(r.get("ring", False)),
                    "inclination_deg": float(r.get("inclination_deg", 0) or 0),
                    "confidence": float(r.get("confidence", 0.2) or 0.2),
                    "notes": r.get("notes", "empty fields filled")
                }
                rows.append(row)

                if copy_dir:
                    try:
                        copy_with_label(p, copy_dir, row["class"], row["bar"], args.bar_in_name)
                    except Exception as ce:
                        row["notes"] = (row.get("notes","") + f" | COPY WARN: {repr(ce)}").strip()

        except Exception as e:
            for p in batch_paths:
                row = {
                    "path": p, "class": "Unknown", "bar": False, "ring": False,
                    "inclination_deg": 0, "confidence": 0.2, "notes": f"INVALID FORMAT: {repr(e)}"
                }
                rows.append(row)
                if copy_dir:
                    try:
                        copy_with_label(p, copy_dir, row["class"], row["bar"], args.bar_in_name)
                    except Exception:
                        pass

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

    print(f"Done. Results in: {args.out}")
    if copy_dir:
        print(f"Copied/renamed images in: {copy_dir.resolve()}")


if __name__ == "__main__":
    classify()
