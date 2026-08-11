"""
Matlab-oracle parity tests for sig2dna_core.tools.

Fixtures were produced by the reference Matlab chain (monotone.m,
filtzero.m, monotonepeak.m, removepeaks.m, monotonepeakfit.m — "MS"
toolbox, O. Vitrac) on the synthetic signals of
``tools_oracle_input.json``; indices in ``tools_oracle.json`` are 1-based
(Matlab) and converted here.

Parity levels:
- segments (monotone): exact (integer indices, float heights to 1e-12)
- filtzero: numerical (both sides implement the same filtfilt algorithm)
- peaks / baseline / deconvolution: added with their modules

Author: Olivier Vitrac, PhD, HDR — Adservio Innovation Lab — Adservio Group —
olivier.vitrac@gmail.com
"""

import json
import os

import numpy as np
import pytest

from sig2dna_core.tools.segments import filtzero, monotone, monotone_full

HERE = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def oracle():
    with open(os.path.join(HERE, "fixtures", "tools_oracle.json")) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def inputs():
    with open(os.path.join(HERE, "fixtures", "tools_oracle_input.json")) as f:
        return json.load(f)


def _arr(v):
    return np.atleast_1d(np.asarray(v, dtype=float))


def assert_segments_match(seg, ref):
    """Python MonotoneSegments vs Matlab {pos,width,height} (1-based)."""
    np.testing.assert_array_equal(seg.start, _arr(ref["pos"]).astype(int) - 1)
    np.testing.assert_array_equal(seg.width, _arr(ref["width"]).astype(int))
    np.testing.assert_allclose(seg.height, _arr(ref["height"]), rtol=0, atol=1e-12)


def assert_full_match(full, ref):
    """Python MonotoneFull vs Matlab full2struct output (1-based)."""
    assert full.nclass == int(ref["nclass"])
    np.testing.assert_array_equal(full.start, _arr(ref["start"]).astype(int) - 1)
    np.testing.assert_array_equal(full.stop, _arr(ref["stop"]).astype(int) - 1)
    np.testing.assert_array_equal(full.width, _arr(ref["width"]).astype(int))
    np.testing.assert_allclose(full.height, _arr(ref["height"]), rtol=0, atol=1e-9)
    np.testing.assert_allclose(full.area, _arr(ref["area"]), rtol=0, atol=1e-9)
    np.testing.assert_array_equal(full.seg_class, _arr(ref["class"]).astype(int) - 1)
    np.testing.assert_array_equal(
        full.fullclass, _arr(ref["fullclass"]).astype(int) - 1
    )
    np.testing.assert_array_equal(
        full.attributes.to_numpy(),
        np.atleast_2d(np.asarray(ref["attributes"], dtype=bool)),
    )


# ---------------------------------------------------------------------------
# monotone parity — doc example
# ---------------------------------------------------------------------------
class TestMonotoneParity:
    def test_plus(self, oracle, inputs):
        y = _arr(inputs["docexample"]["y"])
        assert_segments_match(monotone(y, "+"), oracle["docexample"]["plus"])

    def test_minus(self, oracle, inputs):
        y = _arr(inputs["docexample"]["y"])
        assert_segments_match(monotone(y, "-"), oracle["docexample"]["minus"])

    def test_pm(self, oracle, inputs):
        y = _arr(inputs["docexample"]["y"])
        assert_segments_match(monotone(y, "+-"), oracle["docexample"]["pm"])

    def test_full(self, oracle, inputs):
        y = _arr(inputs["docexample"]["y"])
        assert_full_match(monotone_full(y), oracle["docexample"]["full"])


# ---------------------------------------------------------------------------
# monotone parity — issue-history edge cases
# ---------------------------------------------------------------------------
class TestEdgeParity:
    @pytest.mark.parametrize("i", range(1, 7))
    def test_edge(self, oracle, inputs, i):
        x = _arr(inputs["edges"][i - 1])
        assert_full_match(monotone_full(x), oracle["edges"][f"edge_{i}"])


# ---------------------------------------------------------------------------
# filtzero parity
# ---------------------------------------------------------------------------
class TestFiltzeroParity:
    @pytest.mark.parametrize("m,key", [(5, "filtzero_m5"), (21, "filtzero_m21")])
    def test_docexample(self, oracle, inputs, m, key):
        y = _arr(inputs["docexample"]["y"])
        np.testing.assert_allclose(
            filtzero(y, m), _arr(oracle["docexample"][key]), rtol=1e-9, atol=1e-12
        )

    def test_chromatogram(self, oracle, inputs):
        y = _arr(inputs["chromatogram"]["y"])
        np.testing.assert_allclose(
            filtzero(y, 9),
            _arr(oracle["chromatogram"]["filtzero_m9"]),
            rtol=1e-9,
            atol=1e-12,
        )
