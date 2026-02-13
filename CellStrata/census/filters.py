from __future__ import annotations
from .config import ObsFilters


def _q(s: str) -> str:
    return "'" + str(s).replace("'", "\\'") + "'"


def _in(field: str, values: list[str]) -> str:
    values = [v for v in (values or []) if v is not None and str(v) != ""]
    if not values:
        return ""
    if len(values) == 1:
        return f"{field} == {_q(values[0])}"
    return f"{field} in [{', '.join(_q(v) for v in values)}]"


def build_value_filter(f: ObsFilters) -> str:
    parts: list[str] = []

    parts.append(f"is_primary_data == {str(bool(f.is_primary_data))}")

    if f.dataset_id:
        parts.append(f"dataset_id == {_q(f.dataset_id)}")

    for field, values in [
        ("sex", f.sex),
        ("disease", f.disease),
        ("tissue_general", f.tissue_general),
        ("cell_type", f.cell_type),
        ("assay", f.assay),
    ]:
        expr = _in(field, values)
        if expr:
            parts.append(expr)

    if f.extra_value_filter:
        parts.append(f"({f.extra_value_filter})")

    return " and ".join(parts) if parts else "True"
