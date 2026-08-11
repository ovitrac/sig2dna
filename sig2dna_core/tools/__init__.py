"""
sig2dna_core.tools — signal-analysis toolchain (peaks, baseline, deconvolution).

Python port of the reference Matlab chain of the thesis of J. Kermorvant
(https://theses.hal.science/tel-04194172):

- ``monotone`` / ``filtzero``     → :mod:`.segments` (segment analysis,
  zero-phase smoothing, letter coding, control plots)
- ``monotonepeak``                → :mod:`.peaks` (peak detection)
- ``removepeaks`` / ``ndf``       → :mod:`.baseline` (baseline rebuild)
- ``monotonepeakfit``             → :mod:`.peakfit` (Gaussian/Lorentzian
  deconvolution)

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

from __future__ import annotations

from .baseline import ndf, removepeaks
from .peakfit import FitResult, monotonepeakfit, plot_fit
from .peaks import PeakTable, monotonepeak, plot_peaks
from .segments import (DEFAULT_LETTER_RULES, MonotoneFull, MonotoneSegments,
                       filtzero, monotone, monotone_full, plot_classes,
                       plot_segments, to_letters)

__all__ = [
    "DEFAULT_LETTER_RULES",
    "FitResult",
    "MonotoneFull",
    "MonotoneSegments",
    "PeakTable",
    "filtzero",
    "monotone",
    "monotone_full",
    "monotonepeak",
    "monotonepeakfit",
    "ndf",
    "plot_classes",
    "plot_fit",
    "plot_peaks",
    "plot_segments",
    "removepeaks",
    "to_letters",
]
