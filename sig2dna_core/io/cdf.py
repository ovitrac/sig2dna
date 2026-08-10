"""
ANDI/AIA netCDF (``.cdf``) reader — the standard exchange format for
low-resolution GC-MS chromatograms (ASTM E1947).

Files are NetCDF-3 classic (magic ``CDF\\x01``), read with
``scipy.io.netcdf_file`` — no external netCDF library required.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import numpy as np
from scipy.io import netcdf_file

from .model import GCMSRun

__all__ = ["read_cdf"]

_REQUIRED = ("scan_acquisition_time", "mass_values", "intensity_values")


def _decode(value):
    if isinstance(value, bytes):
        return value.decode("latin1", "replace").strip("\x00").strip()
    return value


def read_cdf(path: str) -> GCMSRun:
    """Read an ANDI-MS ``.cdf`` file losslessly into a :class:`GCMSRun`."""
    nc = netcdf_file(str(path), "r", mmap=False)
    try:
        var = nc.variables
        missing = [k for k in _REQUIRED if k not in var]
        if missing:
            raise ValueError(f"{path}: not an ANDI-MS file (missing {missing})")

        def arr(name, dtype=np.float64):
            return np.asarray(var[name][:], dtype=dtype) if name in var else None

        scan_times = arr("scan_acquisition_time")
        mass_values = arr("mass_values")
        intensity_values = arr("intensity_values")

        # ANDI allows a scale factor attribute on intensities
        scale = getattr(var["intensity_values"], "scale_factor", None)
        if scale is not None and float(scale) not in (0.0, 1.0):
            intensity_values = intensity_values * float(scale)

        n_scans = len(scan_times)
        point_count = arr("point_count", np.int64)
        scan_index = arr("scan_index", np.int64)
        if point_count is None:
            if scan_index is None:
                raise ValueError(f"{path}: neither point_count nor scan_index present")
            point_count = np.diff(np.append(scan_index, len(mass_values)))
        if scan_index is None:
            scan_index = np.concatenate(([0], np.cumsum(point_count[:-1])))

        metadata = {k: _decode(v) for k, v in nc._attributes.items()}
        metadata.update(
            {
                "n_scans": int(n_scans),
                "mass_range": (
                    float(mass_values.min()) if mass_values.size else None,
                    float(mass_values.max()) if mass_values.size else None,
                ),
            }
        )

        return GCMSRun(
            source=str(path),
            format="cdf",
            scan_times=scan_times,
            mass_values=mass_values,
            intensity_values=intensity_values,
            scan_index=scan_index,
            point_count=point_count,
            tic=arr("total_intensity"),
            metadata=metadata,
        )
    finally:
        nc.close()
