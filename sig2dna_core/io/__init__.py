"""
sig2dna_core.io — lossless ingestion of GC-MS acquisitions.

Formats covered (the complete set used in the PhD of J. Kermorvant,
https://theses.hal.science/tel-04194172):

- ANDI/AIA netCDF ``.cdf``            → :func:`read_cdf`
- mzXML (``.mzXML`` **or plain ``.xml``** — the LNE/SCION chain exports
  mzXML with an ``.xml`` extension)   → :func:`read_mzxml`
- Matlab ``chrom`` archives ``.mat``  → :func:`read_mat`
- Agilent ChemStation ``data.ms`` / ``.D`` folders → :func:`read_agilent_ms`,
  :func:`read_agilent_d`
- NIST-format ``.MSP`` spectral libraries → :func:`read_msp`

``read(path)`` dispatches on extension and content sniffing and returns a
:class:`GCMSRun` (or a library list for ``.msp``). Regularization to unit-mass
channels on a regular RT grid is provided by :func:`bin_uniform`.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import os

from .model import GCMSRun
from .cdf import read_cdf
from .mzxml import read_mzxml
from .matlab import read_mat
from .agilent import read_agilent_ms, read_agilent_d
from .binning import bin_uniform
from .msp import read_msp, iter_msp

__all__ = [
    "GCMSRun",
    "read",
    "read_cdf",
    "read_mzxml",
    "read_mat",
    "read_agilent_ms",
    "read_agilent_d",
    "read_msp",
    "iter_msp",
    "bin_uniform",
]

__version__ = "0.1.0"


def _sniff_xml_is_mzxml(path: str) -> bool:
    with open(path, "rb") as f:
        head = f.read(2048)
    return b"mzXML" in head


def read(path: str, **kwargs) -> GCMSRun:
    """Read any supported chromatogram file/folder into a :class:`GCMSRun`.

    Dispatch rules:

    - directory ending in ``.D``  → Agilent acquisition folder
    - ``.cdf``                    → ANDI netCDF
    - ``.mzxml``                  → mzXML
    - ``.xml``                    → mzXML if the header mentions it, else error
    - ``.mat``                    → Matlab ``chrom`` archive
    - ``.ms``                     → Agilent ChemStation data.ms
    """
    p = str(path)
    if os.path.isdir(p):
        if p.rstrip("/\\").lower().endswith(".d"):
            return read_agilent_d(p, **kwargs)
        raise ValueError(f"{p}: directory is not an Agilent .D folder")
    ext = os.path.splitext(p)[1].lower()
    if ext == ".cdf":
        return read_cdf(p, **kwargs)
    if ext == ".mzxml":
        return read_mzxml(p, **kwargs)
    if ext == ".xml":
        if _sniff_xml_is_mzxml(p):
            return read_mzxml(p, **kwargs)
        raise ValueError(
            f"{p}: .xml file is not mzXML (on the LNE chain raw mzXML files "
            "carry the .xml extension; this file has another schema)"
        )
    if ext == ".mat":
        return read_mat(p, **kwargs)
    if ext == ".ms":
        return read_agilent_ms(p, **kwargs)
    raise ValueError(f"{p}: unsupported extension {ext!r}")
