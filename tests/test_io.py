"""
Unit tests for sig2dna_core.io — every fixture is synthetic and generated
in-test, so the public repository carries no real data.

Run: pytest -q tests/test_io.py
"""

from __future__ import annotations

import base64
import struct

import numpy as np
import pytest
from scipy.io import netcdf_file, savemat

from sig2dna_core.io import (
    bin_uniform,
    read,
    read_agilent_ms,
    read_cdf,
    read_mat,
    read_msp,
    read_mzxml,
)

# ----------------------------------------------------------------------------
# Reference synthetic acquisition: 3 scans, unit masses, known intensities
# ----------------------------------------------------------------------------
SCAN_TIMES = np.array([10.0, 10.5, 11.0])           # seconds
SCANS = [                                            # (mass, intensity) pairs
    [(50.0, 100.0), (51.0, 40.0)],
    [(50.0, 120.0), (52.0, 8.0)],
    [(51.0, 60.0)],
]
POINT_COUNT = np.array([2, 2, 1])
TIC = np.array([140.0, 128.0, 60.0])


def _flat():
    mass = np.array([m for scan in SCANS for m, _ in scan])
    inten = np.array([i for scan in SCANS for _, i in scan])
    return mass, inten


def _check_run(run):
    np.testing.assert_allclose(run.scan_times, SCAN_TIMES)
    np.testing.assert_array_equal(run.point_count, POINT_COUNT)
    mass, inten = _flat()
    np.testing.assert_allclose(run.mass_values, mass)
    np.testing.assert_allclose(run.intensity_values, inten)
    np.testing.assert_allclose(run.tic_from_scans(), TIC)


# ------------------------------------------------------------------- fixtures
def make_cdf(path):
    nc = netcdf_file(str(path), "w")
    mass, inten = _flat()
    nc.createDimension("scan_number", len(SCAN_TIMES))
    nc.createDimension("point_number", len(mass))
    for name, dim, data, tc in [
        ("scan_acquisition_time", "scan_number", SCAN_TIMES, "d"),
        ("total_intensity", "scan_number", TIC, "d"),
        ("point_count", "scan_number", POINT_COUNT, "i"),
        ("scan_index", "scan_number", np.array([0, 2, 4]), "i"),
        ("mass_values", "point_number", mass, "d"),
        ("intensity_values", "point_number", inten, "d"),
    ]:
        v = nc.createVariable(name, tc, (dim,))
        v[:] = data
    nc.close()


def make_mzxml(path, precision=32):
    scans = []
    for t, scan, tic in zip(SCAN_TIMES, SCANS, TIC):
        pairs = np.array(scan, dtype=np.float64).ravel()
        payload = base64.b64encode(
            pairs.astype(">f4" if precision == 32 else ">f8").tobytes()
        ).decode()
        scans.append(
            f'<scan num="{len(scans)+1}" msLevel="1" peaksCount="{len(scan)}" '
            f'retentionTime="PT{t}S" totIonCurrent="{tic}">'
            f'<peaks precision="{precision}" byteOrder="network" '
            f'contentType="m/z-int">{payload}</peaks></scan>'
        )
    xml = (
        '<?xml version="1.0" encoding="ISO-8859-1"?>'
        '<mzXML xmlns="http://sashimi.sourceforge.net/schema_revision/mzXML_2.0">'
        f'<msRun scanCount="{len(SCAN_TIMES)}">'
        '<parentFile fileName="file://TEST-PC/synthetic.xml" fileType="RAWData"/>'
        + "".join(scans)
        + "</msRun></mzXML>"
    )
    with open(path, "w") as f:
        f.write(xml)


def make_agilent_ms(path):
    """Craft a minimal ChemStation v2 file per the documented layout."""
    header = bytearray(1024)

    def pascal(off, text):
        b = text.encode("latin1")
        header[off] = len(b)
        header[off + 1 : off + 1 + len(b)] = b

    pascal(0x000, "2")
    pascal(0x004, "GC / MS Data File")
    pascal(0x018, "SYNTH")
    pascal(0x0B2, "01 Jan 26  12:00 pm")
    pascal(0x0D0, "GC MSD")
    pascal(0x0E4, "TEST.M")
    start = len(header)
    struct.pack_into(">H", header, 0x10A, (start + 2) // 2)  # start = v*2-2
    struct.pack_into(">H", header, 0x118, len(SCANS))

    records = bytearray()
    for t, scan, tic in zip(SCAN_TIMES, SCANS, TIC):
        payload = bytearray()
        for m, i in scan:
            packed = int(i) & 0x3FFF  # exponent 0 (intensities < 16384)
            payload += struct.pack(">HH", int(round(m * 20)), packed)
        body_len = 2 + 4 + 6 + 2 + 4 + len(payload) + 4
        rec = bytearray()
        rec += struct.pack(">H", body_len // 2)          # record length (words)
        rec += struct.pack(">I", int(t * 1000))           # time in ms
        rec += b"\x00" * 6
        rec += struct.pack(">H", len(payload) // 4)       # number of peaks
        rec += b"\x00" * 4
        rec += payload
        rec += struct.pack(">I", int(tic))                # TIC at record end
        assert len(rec) == body_len
        records += rec
    with open(path, "wb") as f:
        f.write(bytes(header) + bytes(records))


def make_mat(path):
    mz = np.array([50.0, 51.0, 52.0])
    rt = SCAN_TIMES.copy()
    y = np.array(
        [[100.0, 120.0, 0.0], [40.0, 0.0, 60.0], [0.0, 8.0, 0.0]]
    )
    savemat(
        str(path),
        {"chrom": {"mz": mz, "rt": rt, "ri": rt * 10, "y": y, "rifile": "synthetic"}},
    )


MSP_TEXT = """Name: Synthetic A
CAS#: 000-00-0
MW: 100
Num Peaks: 3
50 999; 51 400; 52 80;

Name: Synthetic B
Num Peaks: 2
60 100
61 50
"""


# ----------------------------------------------------------------------- tests
def test_cdf_roundtrip(tmp_path):
    p = tmp_path / "run.cdf"
    make_cdf(p)
    run = read_cdf(str(p))
    _check_run(run)
    np.testing.assert_allclose(run.tic, TIC)
    assert run.format == "cdf"


def test_mzxml_roundtrip(tmp_path):
    for precision in (32, 64):
        p = tmp_path / f"run{precision}.xml"
        make_mzxml(p, precision=precision)
        run = read_mzxml(str(p))
        _check_run(run)
        np.testing.assert_allclose(run.tic, TIC)
        assert run.metadata["parent_files"] == ["file://TEST-PC/synthetic.xml"]


def test_agilent_roundtrip(tmp_path):
    p = tmp_path / "data.ms"
    make_agilent_ms(p)
    run = read_agilent_ms(str(p))
    _check_run(run)
    np.testing.assert_allclose(run.tic, TIC)
    assert run.metadata["sample_name"] == "SYNTH"
    assert run.metadata["method"] == "TEST.M"


def test_agilent_packed_intensity_exponent(tmp_path):
    """Intensities >= 16384 use the 2-bit exponent (value * 8**exp)."""
    p = tmp_path / "data.ms"
    make_agilent_ms(p)
    raw = bytearray(p.read_bytes())
    # patch first payload intensity to packed value: exp=2, mantissa=1000
    packed = (2 << 14) | 1000
    struct.pack_into(">H", raw, 1024 + 2 + 4 + 6 + 2 + 4 + 2, packed)
    p.write_bytes(bytes(raw))
    run = read_agilent_ms(str(p))
    assert run.intensity_values[0] == 1000 * 64


def test_mat_roundtrip(tmp_path):
    p = tmp_path / "run.mat"
    make_mat(p)
    run = read_mat(str(p))
    assert run.has_matrix and not run.has_sparse
    np.testing.assert_allclose(run.rt, SCAN_TIMES)
    np.testing.assert_allclose(run.y.sum(axis=0), TIC)
    np.testing.assert_allclose(run.ri, SCAN_TIMES * 10)


def test_dispatcher(tmp_path):
    make_cdf(tmp_path / "a.cdf")
    make_mzxml(tmp_path / "b.xml")
    make_mat(tmp_path / "c.mat")
    for fname, fmt in [("a.cdf", "cdf"), ("b.xml", "mzxml"), ("c.mat", "mat")]:
        assert read(str(tmp_path / fname)).format == fmt
    (tmp_path / "notmz.xml").write_text("<html><body>nope</body></html>")
    with pytest.raises(ValueError):
        read(str(tmp_path / "notmz.xml"))
    d = tmp_path / "run.D"
    d.mkdir()
    make_agilent_ms(d / "data.ms")
    run = read(str(d))
    assert run.format == "agilent-ms"
    assert run.metadata["acquisition_folder"].endswith("run.D")


def test_binning_conserves_intensity(tmp_path):
    make_cdf(tmp_path / "a.cdf")
    run = bin_uniform(read_cdf(str(tmp_path / "a.cdf")))
    assert run.has_matrix
    np.testing.assert_array_equal(run.mz, [50.0, 51.0, 52.0])
    # binning itself conserves total intensity (before interpolation)
    assert abs(run.metadata["binning"]["intensity_conservation"] - 1.0) < 1e-12
    # binned matrix at original scan times equals the mat oracle
    make_mat(tmp_path / "c.mat")
    oracle = read_mat(str(tmp_path / "c.mat"))
    for i, t in enumerate(SCAN_TIMES):
        j = int(np.argmin(np.abs(run.rt - t)))
        np.testing.assert_allclose(run.y[:, j], oracle.y[:, i], atol=1e-9)


def test_cross_format_equivalence(tmp_path):
    """The three raw encodings of the same acquisition must agree exactly."""
    make_cdf(tmp_path / "a.cdf")
    make_mzxml(tmp_path / "b.xml")
    d = tmp_path / "run.D"
    d.mkdir()
    make_agilent_ms(d / "data.ms")
    runs = [read(str(tmp_path / "a.cdf")), read(str(tmp_path / "b.xml")), read(str(d))]
    for other in runs[1:]:
        np.testing.assert_allclose(other.scan_times, runs[0].scan_times)
        np.testing.assert_allclose(other.mass_values, runs[0].mass_values)
        np.testing.assert_allclose(other.intensity_values, runs[0].intensity_values)


def test_msp(tmp_path):
    p = tmp_path / "lib.msp"
    p.write_text(MSP_TEXT)
    lib = read_msp(str(p))
    assert [r["name"] for r in lib] == ["Synthetic A", "Synthetic B"]
    np.testing.assert_allclose(lib[0]["mz"], [50, 51, 52])
    np.testing.assert_allclose(lib[0]["intensity"], [999, 400, 80])
    assert lib[0]["fields"]["cas#"] == "000-00-0"
    np.testing.assert_allclose(lib[1]["intensity"], [100, 50])


def test_to_signal_bridge(tmp_path):
    make_cdf(tmp_path / "a.cdf")
    run = read_cdf(str(tmp_path / "a.cdf"))
    s = run.to_signal()
    assert s.type == "GC-MS"
    np.testing.assert_allclose(np.asarray(s.y), TIC)
