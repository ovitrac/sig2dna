"""
Reader for the Matlab ``chrom`` archives produced during the PhD of
J. Kermorvant (MAT-file v5, one ``chrom`` struct per acquisition).

The struct holds the *regularized* chromatogram — the post-ingestion canonical
object of the thesis (§4.2.3): unit-mass channels, a regular retention-time
grid, the Kovats retention-index axis, and optionally extracted spectra:

- ``mz``      (n_mz,)  unit masses (uint16)
- ``rt``      (n_rt,)  seconds, regular grid
- ``ri``      (n_rt,)  Kovats retention indices (from the C7–C40 ladder)
- ``y``       (n_mz, n_rt) ion-chromatogram matrix
- ``spectra`` optional array of structs (extracted peak spectra)
- ``rifile``  provenance of the RI calibration

These files serve both as a fast ingestion path and as the *oracle* against
which the raw readers (cdf, mzXML) are validated.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import mat_struct

from .model import GCMSRun

__all__ = ["read_mat"]


def _struct_to_dict(obj):
    """Recursively convert scipy mat_struct objects to plain dicts/arrays."""
    if isinstance(obj, mat_struct):
        return {f: _struct_to_dict(getattr(obj, f)) for f in obj._fieldnames}
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [_struct_to_dict(v) for v in obj.ravel()]
    return obj


def read_mat(path: str, varname: str = "chrom") -> GCMSRun:
    """Read a Kermorvant ``chrom`` MAT-file into a :class:`GCMSRun`."""
    d = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    if varname not in d:
        candidates = [k for k in d if not k.startswith("__")]
        raise ValueError(
            f"{path}: no '{varname}' struct (variables: {candidates}); "
            "this reader targets the thesis chromatogram archives"
        )
    c = d[varname]
    fields = set(getattr(c, "_fieldnames", []))
    needed = {"mz", "rt", "y"}
    if not needed.issubset(fields):
        raise ValueError(f"{path}: '{varname}' lacks {needed - fields}")

    y = np.asarray(c.y, dtype=np.float64)
    mz = np.asarray(c.mz, dtype=np.float64).ravel()
    rt = np.asarray(c.rt, dtype=np.float64).ravel()
    if y.shape != (mz.size, rt.size):
        if y.shape == (rt.size, mz.size):
            y = y.T
        else:
            raise ValueError(
                f"{path}: y shape {y.shape} inconsistent with mz ({mz.size}) "
                f"and rt ({rt.size})"
            )

    ri = np.asarray(c.ri, dtype=np.float64).ravel() if "ri" in fields else None
    spectra = _struct_to_dict(c.spectra) if "spectra" in fields else None
    if spectra is not None and not isinstance(spectra, list):
        spectra = [spectra]

    metadata = {
        "n_scans": int(rt.size),
        "n_mz": int(mz.size),
        "varname": varname,
    }
    if "rifile" in fields:
        metadata["rifile"] = str(c.rifile)

    return GCMSRun(
        source=str(path),
        format="mat",
        mz=mz,
        rt=rt,
        y=y,
        ri=ri,
        spectra=spectra,
        metadata=metadata,
    )
