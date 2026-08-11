"""
NIST-format ``.MSP`` spectral-library parser (text records: ``Name:`` …
``Num Peaks:`` followed by mass/intensity pairs separated by ';' or
whitespace). Covers the NIST chunk exports, MassBank/MoNA exports, and
in-house libraries following the same convention.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import re
from typing import Iterator

import numpy as np

__all__ = ["read_msp", "iter_msp"]

_PAIR = re.compile(r"([\d.eE+-]+)[\s,]+([\d.eE+-]+)")


def iter_msp(path: str, encoding: str = "latin1") -> Iterator[dict]:
    """Yield one dict per library record: name, fields, mz, intensity."""
    record: dict = {}
    peaks: list = []
    expected = None

    def flush():
        nonlocal record, peaks, expected
        if record.get("name") or peaks:
            arr = np.asarray(peaks, dtype=np.float64).reshape(-1, 2)
            if expected is not None and arr.shape[0] != expected:
                raise ValueError(
                    f"{path}: record {record.get('name')!r} declares "
                    f"{expected} peaks but lists {arr.shape[0]}"
                )
            yield_rec = {
                "name": record.pop("name", ""),
                "fields": record,
                "mz": arr[:, 0],
                "intensity": arr[:, 1],
            }
            record, peaks, expected = {}, [], None
            return yield_rec
        record, peaks, expected = {}, [], None
        return None

    with open(path, "r", encoding=encoding, errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                rec = flush()
                if rec:
                    yield rec
                continue
            m = re.match(r"^([A-Za-z][A-Za-z0-9 _#/.\-]*?):\s*(.*)$", line)
            if m and expected is None:
                key, val = m.group(1).strip().lower(), m.group(2).strip()
                if key == "name":
                    rec = flush()
                    if rec:
                        yield rec
                    record["name"] = val
                elif key in ("num peaks", "numpeaks"):
                    expected = int(val)
                else:
                    record[key] = val
            else:
                for pm in _PAIR.finditer(line):
                    peaks.append((float(pm.group(1)), float(pm.group(2))))
    rec = flush()
    if rec:
        yield rec


def read_msp(path: str, encoding: str = "latin1") -> list:
    """Read a whole ``.MSP`` library into a list of records."""
    return list(iter_msp(path, encoding=encoding))
