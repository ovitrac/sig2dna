"""
Regularization of sparse GC-MS acquisitions: unit-mass binning + interpolation
on a regular retention-time grid.

Implements the reduction of the thesis (§4.2.3): raw acquisitions are
irregular in both m/z and time; low-resolution single-quadrupole data are
binned to unit masses over the acquisition range and linearly interpolated on
a regular RT grid whose step is a low percentile of the raw time steps
(1st percentile by default), so that no scan is undersampled.

Unit-mass binning is intensity-conserving (sum of the binned matrix equals the
sum of raw intensities); the time interpolation is the only lossy step, and
its grid is chosen to make that loss negligible.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .model import GCMSRun

__all__ = ["bin_uniform"]


def bin_uniform(
    run: GCMSRun,
    mz_min: Optional[int] = None,
    mz_max: Optional[int] = None,
    dt: Optional[float] = None,
    dt_percentile: float = 1.0,
    inplace: bool = True,
) -> GCMSRun:
    """Fill ``run.mz / run.rt / run.y`` from the sparse representation.

    Parameters
    ----------
    run : GCMSRun with sparse data (cdf / mzxml / agilent readers).
    mz_min, mz_max : acquisition mass range; defaults to the rounded data range.
    dt : regular grid step in seconds; default = ``dt_percentile`` percentile
        of the raw scan-time steps (thesis §4.2.3 uses the 1st percentile).
    inplace : mutate ``run`` (default) or return a regularized copy.
    """
    if not run.has_sparse:
        raise ValueError(f"{run.source}: no sparse data to regularize")

    mass_bins = np.rint(run.mass_values).astype(np.int64)
    lo = int(mz_min) if mz_min is not None else int(mass_bins.min())
    hi = int(mz_max) if mz_max is not None else int(mass_bins.max())
    keep = (mass_bins >= lo) & (mass_bins <= hi)

    n_scans = run.n_scans
    scan_of_point = np.repeat(np.arange(n_scans), run.point_count)

    n_mz = hi - lo + 1
    per_scan = np.zeros((n_mz, n_scans), dtype=np.float64)
    np.add.at(
        per_scan,
        (mass_bins[keep] - lo, scan_of_point[keep]),
        run.intensity_values[keep],
    )

    t = np.asarray(run.scan_times, dtype=np.float64)
    steps = np.diff(t)
    steps = steps[steps > 0]
    if steps.size == 0:
        raise ValueError(f"{run.source}: degenerate time axis")
    step = float(dt) if dt is not None else float(np.percentile(steps, dt_percentile))
    rt = np.arange(t[0], t[-1] + 0.5 * step, step)

    y = np.empty((n_mz, rt.size), dtype=np.float64)
    for i in range(n_mz):
        y[i, :] = np.interp(rt, t, per_scan[i, :])

    target = run if inplace else GCMSRun(
        source=run.source,
        format=run.format,
        scan_times=run.scan_times,
        mass_values=run.mass_values,
        intensity_values=run.intensity_values,
        scan_index=run.scan_index,
        point_count=run.point_count,
        tic=run.tic,
        metadata=dict(run.metadata),
    )
    target.mz = np.arange(lo, hi + 1, dtype=np.float64)
    target.rt = rt
    target.y = y
    target.metadata.update(
        {
            "binning": {
                "mz_min": lo,
                "mz_max": hi,
                "dt": step,
                "dt_percentile": None if dt is not None else dt_percentile,
                "n_rt": int(rt.size),
                "intensity_conservation": float(
                    per_scan.sum() / run.intensity_values[keep].sum()
                )
                if run.intensity_values[keep].sum() > 0
                else 1.0,
            }
        }
    )
    return target
