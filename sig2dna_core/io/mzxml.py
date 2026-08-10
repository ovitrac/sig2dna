"""
mzXML reader (versions 2.x/3.x), streaming and lossless.

Important: on the LNE/SCION acquisition chain the mzXML files carry a plain
``.xml`` extension. The dispatcher in ``sig2dna_core.io`` sniffs content, so
both ``.mzXML`` and ``.xml`` reach this reader.

Peaks are base64-encoded (m/z, intensity) pairs, big-endian ("network") float32
or float64, optionally zlib-compressed.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import base64
import re
import zlib
import xml.etree.ElementTree as ET

import numpy as np

from .model import GCMSRun

__all__ = ["read_mzxml"]

_ISO_DURATION = re.compile(r"^PT(?:(\d+(?:\.\d+)?)M)?(?:(\d+(?:\.\d+)?)S)?$")


def _parse_time_s(value: str) -> float:
    """Parse an xsd:duration retention time ('PT150.124S', 'PT2M30S') to seconds."""
    m = _ISO_DURATION.match(value.strip())
    if not m:
        raise ValueError(f"unsupported retentionTime format: {value!r}")
    minutes = float(m.group(1) or 0.0)
    seconds = float(m.group(2) or 0.0)
    return 60.0 * minutes + seconds


def _localname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _decode_peaks(text: str, precision: int, byte_order: str, compression: str):
    raw = base64.b64decode(text)
    if compression and compression.lower() == "zlib":
        raw = zlib.decompress(raw)
    dtype = {32: ">f4", 64: ">f8"}[int(precision)]
    if byte_order not in ("network", "big", None, ""):
        dtype = dtype.replace(">", "<")
    pairs = np.frombuffer(raw, dtype=dtype).astype(np.float64)
    if pairs.size % 2:
        raise ValueError("mzXML peaks payload has an odd number of values")
    pairs = pairs.reshape(-1, 2)
    return pairs[:, 0], pairs[:, 1]


def read_mzxml(path: str) -> GCMSRun:
    """Stream-parse an mzXML file losslessly into a :class:`GCMSRun`."""
    scan_times, masses, intensities, point_count, tic = [], [], [], [], []
    metadata: dict = {}
    n_declared = None

    context = ET.iterparse(str(path), events=("start", "end"))
    _, root = next(context)  # grab root to free memory as we go
    scan_attrs = None
    peak_attrs = None
    for event, elem in context:
        name = _localname(elem.tag)
        if event == "start":
            if name == "msRun":
                n_declared = elem.get("scanCount")
            elif name == "scan":
                scan_attrs = dict(elem.attrib)
            continue
        # end events
        if name == "parentFile":
            metadata.setdefault("parent_files", []).append(
                elem.get("fileName") or elem.get("fileURI") or ""
            )
        elif name in ("msInstrument", "msManufacturer", "msModel", "software"):
            for k, v in elem.attrib.items():
                metadata.setdefault(name, {})[k] = v
        elif name == "peaks":
            peak_attrs = dict(elem.attrib)
            peak_text = elem.text or ""
            if scan_attrs is None:
                raise ValueError(f"{path}: <peaks> outside of <scan>")
            if int(scan_attrs.get("msLevel", "1")) != 1:
                peak_attrs = None  # MSn not expected on this chain; skip
                continue
            mz, inten = _decode_peaks(
                peak_text,
                precision=int(peak_attrs.get("precision", "32")),
                byte_order=peak_attrs.get("byteOrder", "network"),
                compression=peak_attrs.get("compressionType", ""),
            )
            declared = scan_attrs.get("peaksCount")
            if declared is not None and int(declared) != len(mz):
                raise ValueError(
                    f"{path}: scan {scan_attrs.get('num')}: peaksCount="
                    f"{declared} but decoded {len(mz)} peaks"
                )
            scan_times.append(_parse_time_s(scan_attrs["retentionTime"]))
            masses.append(mz)
            intensities.append(inten)
            point_count.append(len(mz))
            tic.append(float(scan_attrs.get("totIonCurrent", "nan")))
        elif name == "scan":
            scan_attrs = None
            elem.clear()
            root.clear()

    if not scan_times:
        raise ValueError(f"{path}: no MS1 scans found (is this an mzXML file?)")

    order = np.argsort(scan_times, kind="stable")
    scan_times_a = np.asarray(scan_times, dtype=np.float64)[order]
    masses = [masses[i] for i in order]
    intensities = [intensities[i] for i in order]
    point_count_a = np.asarray(point_count, dtype=np.int64)[order]
    tic_a = np.asarray(tic, dtype=np.float64)[order]

    metadata.update({"n_scans": int(len(scan_times_a)), "declared_scan_count": n_declared})
    return GCMSRun(
        source=str(path),
        format="mzxml",
        scan_times=scan_times_a,
        mass_values=np.concatenate(masses) if masses else np.empty(0),
        intensity_values=np.concatenate(intensities) if intensities else np.empty(0),
        scan_index=np.concatenate(([0], np.cumsum(point_count_a[:-1]))),
        point_count=point_count_a,
        tic=None if np.all(np.isnan(tic_a)) else tic_a,
        metadata=metadata,
    )
