"""Tests for :mod:`selene_base.data.pds4_table`."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from selene_base.data.pds4_table import (
    FieldSpec,
    TableSpec,
    parse_pds4_label,
    read_pds4_table,
)

REAL_XML = Path("data/raw/diviner/dlre_prp_south.xml")


@pytest.mark.skipif(not REAL_XML.exists(), reason="PRP XML not downloaded")
class TestParsePDS4Label:
    def test_top_level_metadata(self) -> None:
        spec = parse_pds4_label(REAL_XML)
        assert spec["header_offset"] == 210
        assert spec["record_count"] == 2_880_000
        assert spec["record_length"] == 210
        assert len(spec["fields"]) == 15

    def test_known_field_offsets_match_spec(self) -> None:
        spec = parse_pds4_label(REAL_XML)
        by_name = {f["name"]: f for f in spec["fields"]}
        # Spec values from the prompt; offsets are converted to 0-indexed.
        assert by_name["tri_clon"]["offset_0"] == 126
        assert by_name["tri_clon"]["length"] == 15
        assert by_name["tri_clat"]["offset_0"] == 143
        assert by_name["tri_clat"]["length"] == 15
        assert by_name["temp_avg"]["offset_0"] == 177
        assert by_name["temp_avg"]["length"] == 8
        assert by_name["temp_max"]["offset_0"] == 187
        assert by_name["temp_max"]["length"] == 8
        assert by_name["ice_depth"]["offset_0"] == 197
        assert by_name["ice_depth"]["length"] == 11
        assert -999.0 in by_name["ice_depth"]["invalid_constants"]


def _synthetic_spec() -> TableSpec:
    return TableSpec(
        header_offset=10,
        record_count=3,
        record_length=22,
        fields=[
            FieldSpec(
                name="x",
                offset_0=0,
                length=8,
                data_type="ASCII_Real",
                unit=None,
                invalid_constants=[],
            ),
            FieldSpec(
                name="y",
                offset_0=8,
                length=8,
                data_type="ASCII_Real",
                unit=None,
                invalid_constants=[],
            ),
            FieldSpec(
                name="z",
                offset_0=16,
                length=4,
                data_type="ASCII_Real",
                unit=None,
                invalid_constants=[-999.0],
            ),
        ],
    )


def _write_synthetic_tab(tmp_path: Path) -> Path:
    spec = _synthetic_spec()
    # x: 8 chars, y: 8 chars, z: 4 chars, then CRLF -> 22 bytes total.
    rows = [
        b"   1.000   2.0003.00\r\n",
        b"   4.500  -1.500-999\r\n",
        b"   0.000   0.0000.50\r\n",
    ]
    for row in rows:
        assert len(row) == spec["record_length"]
    payload = b"HEADER____" + b"".join(rows)  # 10-byte header + 3 records
    tab = tmp_path / "synthetic.tab"
    tab.write_bytes(payload)
    return tab


class TestReadPDS4Table:
    def test_round_trip_synthetic(self, tmp_path: Path) -> None:
        spec = _synthetic_spec()
        path = _write_synthetic_tab(tmp_path)
        df = read_pds4_table(path, spec)
        assert list(df.columns) == ["x", "y", "z"]
        assert len(df) == 3
        np.testing.assert_allclose(df["x"], [1.0, 4.5, 0.0])
        np.testing.assert_allclose(df["y"], [2.0, -1.5, 0.0])
        # row 2: z=-999 sentinel -> NaN
        assert np.isnan(df["z"].iloc[1])
        np.testing.assert_allclose(df["z"].iloc[[0, 2]], [3.0, 0.5])

    def test_subset_fields(self, tmp_path: Path) -> None:
        spec = _synthetic_spec()
        path = _write_synthetic_tab(tmp_path)
        df = read_pds4_table(path, spec, fields=["y"])
        assert list(df.columns) == ["y"]
        assert len(df) == 3

    def test_unknown_field_rejected(self, tmp_path: Path) -> None:
        spec = _synthetic_spec()
        path = _write_synthetic_tab(tmp_path)
        with pytest.raises(ValueError, match="unknown field"):
            read_pds4_table(path, spec, fields=["bogus"])

    def test_chunked_read_matches_unchunked(self, tmp_path: Path) -> None:
        spec = _synthetic_spec()
        path = _write_synthetic_tab(tmp_path)
        df_one = read_pds4_table(path, spec, chunk_rows=10)
        df_two = read_pds4_table(path, spec, chunk_rows=1)
        # Both should produce identical output regardless of chunk size.
        for col in df_one.columns:
            np.testing.assert_array_equal(
                df_one[col].fillna(-1).to_numpy(),
                df_two[col].fillna(-1).to_numpy(),
            )
