"""Streaming reader for PDS4 ``Table_Character`` fixed-width ASCII tables.

Used by the Diviner Polar Resource Product (``dlre_prp_south.tab``) —
2.88 M records of 210 bytes each, ~605 MB on disk. Tested against that
specific schema; should generalise to any PDS4 character table whose
label declares one ``Record_Character`` with N ``Field_Character``
children.

Two functions:

- :func:`parse_pds4_label` reads the XML label and returns a
  parsing-ready spec.
- :func:`read_pds4_table` streams the ``.tab`` file and returns a
  pandas ``DataFrame`` of the requested columns. Invalid sentinels
  declared in ``Special_Constants/invalid_constant`` are converted to
  ``NaN``.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import TypedDict

import numpy as np
import pandas as pd

PDS4_NS = "http://pds.nasa.gov/pds4/pds/v1"
NS = {"pds": PDS4_NS}


class FieldSpec(TypedDict):
    name: str
    offset_0: int  # zero-indexed byte offset within a record
    length: int  # bytes
    data_type: str
    unit: str | None
    invalid_constants: list[float]


class TableSpec(TypedDict):
    header_offset: int  # bytes to skip before the first record
    record_count: int
    record_length: int  # includes the trailing record delimiter (e.g. \r\n)
    fields: list[FieldSpec]


def parse_pds4_label(xml_path: Path) -> TableSpec:
    """Parse a PDS4 XML label into a :class:`TableSpec`.

    Args:
        xml_path: Path to the ``.xml`` label.

    Returns:
        :class:`TableSpec` with one ``FieldSpec`` per declared
        ``Field_Character``. Byte offsets are converted from PDS4's
        1-indexed convention to 0-indexed Python slices.

    Raises:
        ValueError: If no ``Table_Character`` is found in the label.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    table = root.find(".//pds:Table_Character", NS)
    if table is None:
        raise ValueError(f"no Table_Character in {xml_path}")

    header_offset = int(_text(table, "offset"))
    record_count = int(_text(table, "records"))
    rec = table.find("pds:Record_Character", NS)
    if rec is None:
        raise ValueError(f"no Record_Character in {xml_path}")
    record_length = int(_text(rec, "record_length"))

    fields: list[FieldSpec] = []
    for field in rec.findall("pds:Field_Character", NS):
        name = _text(field, "name")
        offset_1 = int(_text(field, "field_location"))
        length = int(_text(field, "field_length"))
        data_type = _text(field, "data_type")
        unit_el = field.find("pds:unit", NS)
        unit = unit_el.text if unit_el is not None and unit_el.text else None
        invalid_constants = [
            float(s.text)
            for s in field.findall(".//pds:Special_Constants/pds:invalid_constant", NS)
            if s.text is not None
        ]
        fields.append(
            FieldSpec(
                name=name,
                offset_0=offset_1 - 1,
                length=length,
                data_type=data_type,
                unit=unit,
                invalid_constants=invalid_constants,
            )
        )

    return TableSpec(
        header_offset=header_offset,
        record_count=record_count,
        record_length=record_length,
        fields=fields,
    )


def _text(element: ET.Element, tag: str) -> str:
    el = element.find(f"pds:{tag}", NS)
    if el is None or el.text is None:
        raise ValueError(f"missing element <{tag}> in {element.tag}")
    return el.text


def _parse_field_block(
    block: bytes, field: FieldSpec, record_length: int, n_records: int
) -> np.ndarray:
    """Slice one field's column out of a contiguous binary record block."""
    arr = np.frombuffer(block, dtype=np.uint8).reshape(n_records, record_length)
    end = field["offset_0"] + field["length"]
    column_bytes = arr[:, field["offset_0"] : end]
    text_view = column_bytes.tobytes().decode("ascii")
    width = field["length"]
    pieces = [text_view[i * width : (i + 1) * width] for i in range(n_records)]
    out = np.fromiter(
        (float(p) if p.strip() else np.nan for p in pieces),
        dtype=np.float64,
        count=n_records,
    )
    if field["invalid_constants"]:
        for sentinel in field["invalid_constants"]:
            out[out == sentinel] = np.nan
    return out


def read_pds4_table(
    tab_path: Path,
    spec: TableSpec,
    *,
    fields: list[str] | None = None,
    chunk_rows: int = 100_000,
) -> pd.DataFrame:
    """Stream-parse a PDS4 fixed-width ASCII table.

    The reader walks the file in fixed-byte chunks of ``chunk_rows``
    records and slices each requested field as a 2-D byte view; this
    keeps peak memory bounded even for large tables (the PRP south is
    ~605 MB raw).

    Args:
        tab_path: Path to the ``.tab`` data file.
        spec: Parsing spec from :func:`parse_pds4_label`.
        fields: Subset of column names to extract. ``None`` reads
            every declared column.
        chunk_rows: Records read per pass.

    Returns:
        :class:`pandas.DataFrame` with one column per requested field.
        Cells matching any of the field's declared invalid constants
        are :py:data:`numpy.nan`.

    Raises:
        ValueError: If a requested field name is not declared in
            ``spec``.
    """
    requested = [f for f in spec["fields"] if fields is None or f["name"] in fields]
    if fields is not None:
        missing = set(fields) - {f["name"] for f in spec["fields"]}
        if missing:
            raise ValueError(f"unknown field name(s): {sorted(missing)}")

    record_length = spec["record_length"]
    n_total = spec["record_count"]

    columns: dict[str, list[np.ndarray]] = {f["name"]: [] for f in requested}
    with tab_path.open("rb") as fh:
        fh.seek(spec["header_offset"])
        rows_left = n_total
        while rows_left > 0:
            n = min(chunk_rows, rows_left)
            block = fh.read(n * record_length)
            actual = len(block) // record_length
            if actual == 0:
                break
            for field in requested:
                columns[field["name"]].append(
                    _parse_field_block(block, field, record_length, actual)
                )
            rows_left -= actual

    return pd.DataFrame({name: np.concatenate(parts) for name, parts in columns.items()})
