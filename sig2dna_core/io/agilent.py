"""
Agilent/HP ChemStation ``data.ms`` reader (GC single-quadrupole, file
version "2" — magic ``\\x012`` + "GC / MS Data File").

Binary layout (big-endian throughout; format facts cross-checked against the
MIT-licensed Chemplexity Chromatography Toolbox and verified on LNE files):

Header
    0x000  pascal string  file version ("2")
    0x004  pascal string  description ("GC / MS Data File")
    0x018  pascal string  sample name
    0x0B2  pascal string  acquisition date/time (e.g. "17 Jan 22  05:23 pm")
    0x0D0  pascal string  instrument ("GC MSD")
    0x0E4  pascal string  method name
    0x10A  uint16         data start = value * 2 - 2 (byte offset)
    0x118  uint16         number of scans

Scan record (at ``pos``)
    +0x00  uint16  record length in 16-bit words → next record at pos + 2*len
    +0x02  uint32  retention time in milliseconds
    +0x06  6 bytes skipped
    +0x0C  uint16  number of peaks n_peaks (payload = 2 × n_peaks uint16)
    +0x0E  4 bytes skipped, then payload:
           n_peaks × (uint16 mass_x20, uint16 packed_intensity)
    next_pos-4  uint32  total ion current of the scan

Decoding
    mass      = mass_x20 / 20
    intensity = (packed & 0x3FFF) * 8 ** (packed >> 14)

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import os
import struct

import numpy as np

from .model import GCMSRun

__all__ = ["read_agilent_ms", "read_agilent_d"]

_HDR_STRINGS = {
    0x018: "sample_name",
    0x0B2: "datetime",
    0x0D0: "instrument",
    0x0E4: "method",
}


def _pascal(buf: bytes, offset: int) -> str:
    n = buf[offset]
    return buf[offset + 1 : offset + 1 + n].decode("latin1", "replace").strip()


def read_agilent_ms(path: str) -> GCMSRun:
    """Read a ChemStation ``data.ms`` (version 2) into a :class:`GCMSRun`."""
    with open(path, "rb") as f:
        buf = f.read()

    version = _pascal(buf, 0x0)
    description = _pascal(buf, 0x4)
    if version != "2" or "MS Data File" not in description:
        raise ValueError(
            f"{path}: unsupported Agilent file (version={version!r}, "
            f"description={description!r}); only GC-MS version 2 is implemented"
        )

    metadata = {k: _pascal(buf, off) for off, k in _HDR_STRINGS.items()}
    metadata["file_version"] = version

    (start_words,) = struct.unpack_from(">H", buf, 0x10A)
    (n_scans,) = struct.unpack_from(">H", buf, 0x118)
    pos = start_words * 2 - 2

    scan_times = np.empty(n_scans, dtype=np.float64)
    tic = np.empty(n_scans, dtype=np.float64)
    point_count = np.empty(n_scans, dtype=np.int64)
    payloads = []

    for i in range(n_scans):
        (rec_words,) = struct.unpack_from(">H", buf, pos)
        next_pos = pos + rec_words * 2
        (time_ms,) = struct.unpack_from(">I", buf, pos + 2)
        scan_times[i] = time_ms / 1000.0  # seconds
        (n_peaks,) = struct.unpack_from(">H", buf, pos + 12)
        data_start = pos + 12 + 2 + 4
        payload = np.frombuffer(
            buf, dtype=">u2", count=2 * n_peaks, offset=data_start
        ).astype(np.uint16)
        payloads.append(payload.reshape(-1, 2))
        point_count[i] = n_peaks
        (tic_i,) = struct.unpack_from(">I", buf, next_pos - 4)
        tic[i] = float(tic_i)
        pos = next_pos

    if payloads:
        stacked = np.concatenate(payloads, axis=0)
        mass_values = stacked[:, 0].astype(np.float64) / 20.0
        packed = stacked[:, 1].astype(np.uint32)
        intensity_values = (packed & 0x3FFF).astype(np.float64) * np.power(
            8.0, (packed >> 14).astype(np.float64)
        )
    else:
        mass_values = np.empty(0)
        intensity_values = np.empty(0)

    metadata["n_scans"] = int(n_scans)
    return GCMSRun(
        source=str(path),
        format="agilent-ms",
        scan_times=scan_times,
        mass_values=mass_values,
        intensity_values=intensity_values,
        scan_index=np.concatenate(([0], np.cumsum(point_count[:-1]))),
        point_count=point_count,
        tic=tic,
        metadata=metadata,
    )


def read_agilent_d(path: str) -> GCMSRun:
    """Read an Agilent ``.D`` acquisition folder (locates ``data.ms``)."""
    for candidate in ("data.ms", "DATA.MS", "Data.ms"):
        ms_path = os.path.join(path, candidate)
        if os.path.isfile(ms_path):
            run = read_agilent_ms(ms_path)
            run.metadata["acquisition_folder"] = str(path)
            return run
    raise FileNotFoundError(f"{path}: no data.ms found in .D folder")
