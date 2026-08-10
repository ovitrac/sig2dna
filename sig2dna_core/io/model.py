"""
Common data model for GC-MS ingestion (sig2dna_core.io).

`GCMSRun` is the single canonical container produced by every reader in this
subpackage. It preserves the acquisition **losslessly** in sparse/centroid form
(`scan_times`, `mass_values`, `intensity_values`, `scan_index`, `point_count`
— the ANDI/AIA layout, adopted here as the common denominator across formats)
and optionally carries the **regularized** form (`mz`, `rt`, `y`, `ri`) as
produced by unit-mass binning + regular-grid interpolation (see
`sig2dna_core.io.binning`) or as read back from a Matlab `chrom` archive.

Bridges `to_signal()` / `to_signal_collection()` hand results to
`sig2dna_core.signomics` without modifying it.

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np

__all__ = ["GCMSRun"]


@dataclass
class GCMSRun:
    """One GC-MS acquisition, format-agnostic.

    Sparse (lossless, native) representation — present for cdf / mzXML /
    Agilent sources, ``None`` for Matlab archives (which store only the
    regularized form):

    - ``scan_times``        (n_scans,) float64, seconds
    - ``mass_values``       (n_points,) float64, flat concatenation of scans
    - ``intensity_values``  (n_points,) float64
    - ``scan_index``        (n_scans,) int64, start offset of each scan
    - ``point_count``       (n_scans,) int64, points per scan
    - ``tic``               (n_scans,) float64, instrument TIC when recorded

    Regularized representation (filled by ``binning.bin_uniform`` or by the
    Matlab reader):

    - ``mz``  (n_mz,) unit masses
    - ``rt``  (n_rt,) seconds, regular grid
    - ``y``   (n_mz, n_rt) ion-chromatogram matrix
    - ``ri``  (n_rt,) Kovats retention-index axis, when available
    """

    source: str
    format: str
    scan_times: Optional[np.ndarray] = None
    mass_values: Optional[np.ndarray] = None
    intensity_values: Optional[np.ndarray] = None
    scan_index: Optional[np.ndarray] = None
    point_count: Optional[np.ndarray] = None
    tic: Optional[np.ndarray] = None
    mz: Optional[np.ndarray] = None
    rt: Optional[np.ndarray] = None
    y: Optional[np.ndarray] = None
    ri: Optional[np.ndarray] = None
    spectra: Optional[list] = None
    metadata: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- checks
    @property
    def has_sparse(self) -> bool:
        return self.mass_values is not None and self.scan_index is not None

    @property
    def has_matrix(self) -> bool:
        return self.y is not None and self.rt is not None and self.mz is not None

    @property
    def n_scans(self) -> int:
        if self.scan_times is not None:
            return int(len(self.scan_times))
        if self.rt is not None:
            return int(len(self.rt))
        return 0

    # ------------------------------------------------------------------ TIC
    def tic_from_scans(self) -> np.ndarray:
        """Total ion current recomputed from the sparse scans (lossless)."""
        if not self.has_sparse:
            if self.y is not None:
                return np.asarray(self.y).sum(axis=0)
            raise ValueError(f"{self.source}: no data to compute a TIC from")
        out = np.zeros(self.n_scans, dtype=np.float64)
        np.add.at(
            out,
            np.repeat(np.arange(self.n_scans), self.point_count),
            self.intensity_values,
        )
        return out

    def get_tic(self, prefer_instrument: bool = True) -> np.ndarray:
        if prefer_instrument and self.tic is not None:
            return np.asarray(self.tic, dtype=np.float64)
        return self.tic_from_scans()

    def eic(self, mass: float, tol: float = 0.5) -> np.ndarray:
        """Extracted-ion chromatogram at ``mass`` ± ``tol`` from sparse data."""
        if self.has_matrix and float(mass) in set(np.asarray(self.mz).tolist()):
            row = int(np.argmin(np.abs(np.asarray(self.mz) - mass)))
            return np.asarray(self.y)[row, :]
        if not self.has_sparse:
            raise ValueError(f"{self.source}: no sparse data for EIC extraction")
        sel = np.abs(self.mass_values - mass) <= tol
        out = np.zeros(self.n_scans, dtype=np.float64)
        scan_of_point = np.repeat(np.arange(self.n_scans), self.point_count)
        np.add.at(out, scan_of_point[sel], self.intensity_values[sel])
        return out

    # ----------------------------------------------------------- provenance
    def sha256(self) -> str:
        """SHA-256 of the source file (lazy; cached in metadata)."""
        cached = self.metadata.get("sha256")
        if cached:
            return cached
        h = hashlib.sha256()
        with open(self.source, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        self.metadata["sha256"] = h.hexdigest()
        return self.metadata["sha256"]

    # ------------------------------------------------------------ signomics
    def to_signal(self, which: str = "tic", name: Optional[str] = None) -> Any:
        """Bridge to ``sig2dna_core.signomics.signal``.

        which: 'tic' (default) or a numeric m/z for a single-ion signal.
        """
        from sig2dna_core.signomics import signal  # deferred, heavy import

        if which == "tic":
            if self.has_matrix:
                x, y = np.asarray(self.rt), np.asarray(self.y).sum(axis=0)
            else:
                x, y = np.asarray(self.scan_times), self.get_tic()
            label = "TIC"
        else:
            mass = float(which)
            if self.has_matrix:
                row = int(np.argmin(np.abs(np.asarray(self.mz) - mass)))
                x, y = np.asarray(self.rt), np.asarray(self.y)[row, :]
            else:
                x, y = np.asarray(self.scan_times), self.eic(mass)
            label = f"m/z {mass:g}"
        return signal(
            x=x,
            y=y,
            name=name or f"{os.path.basename(self.source)} [{label}]",
            type="GC-MS",
            x_label="retention time",
            x_unit="s",
            y_label="intensity",
            y_unit="counts",
            metadata={"source": self.source, "format": self.format, "channel": label},
            source=self.format,
        )

    def to_signal_collection(self, masses: Optional[list] = None, n: int = 1024) -> Any:
        """Per-ion ``signal_collection`` (defaults to every unit mass in ``mz``)."""
        from sig2dna_core.signomics import signal_collection  # deferred

        if not self.has_matrix:
            raise ValueError(
                f"{self.source}: regularize first (io.binning.bin_uniform) "
                "before building a per-ion collection"
            )
        mzv = np.asarray(self.mz)
        rows = (
            range(len(mzv))
            if masses is None
            else [int(np.argmin(np.abs(mzv - m))) for m in masses]
        )
        signals = [self.to_signal(which=float(mzv[i])) for i in rows]
        return signal_collection(*signals, n=n)
