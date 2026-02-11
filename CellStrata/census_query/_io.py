"""
_io.py - I/O utilities, streaming, and TileDB context management.

This module has no internal dependencies. It deals only with Arrow tables,
file paths, and optional TileDB configuration.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    Optional,
    Sequence,
    Tuple,
)

import pyarrow as pa
import pyarrow.parquet as pq

# Optional tiledbsoma for context tuning
try:
    import tiledbsoma as soma
except ImportError:
    soma = None

logger = logging.getLogger(__name__)


def _ensure_writable_parent(path: Path) -> None:
    """Ensure the parent directory exists and is writable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if not os.access(path.parent, os.W_OK):
        raise PermissionError(f"Output directory is not writable: {path.parent}")


def _resolve_outpath(outpath: str, overwrite: bool) -> Path:
    """
    Resolve and validate the output path.

    Args:
        outpath: Target output path.
        overwrite: Whether to overwrite existing files.

    Returns:
        Resolved Path object.

    Raises:
        FileExistsError: If file exists and overwrite is False.
        PermissionError: If file cannot be written/overwritten.
    """
    p = Path(outpath)
    _ensure_writable_parent(p)
    if p.exists():
        if overwrite:
            if not os.access(p, os.W_OK):
                raise PermissionError(f"Cannot overwrite existing file: {p}")
            p.unlink()
        else:
            raise FileExistsError(f"Output exists and overwrite=False: {p}")
    return p


_DEFAULT_S3_CONFIG: Dict[str, Any] = {
    "vfs.s3.connect_timeout_ms": 60_000,       # 60 s connect timeout
    "vfs.s3.request_timeout_ms": 600_000,       # 10 min request timeout
    "vfs.s3.connect_max_tries": 5,              # retry up to 5 times
    "vfs.s3.connect_scale_factor": 25,          # exponential backoff factor
}


def _make_soma_context(tiledb_config: Dict[str, Any]):
    """
    Create a SOMATileDBContext for tuning S3 timeouts/retries.

    Applies sensible default S3 timeout settings even when no user config
    is provided, to prevent stalls on slow or congested S3 connections.

    Args:
        tiledb_config: TileDB configuration dictionary (user overrides).

    Returns:
        SOMATileDBContext, or None if tiledbsoma is not available.
    """
    if soma is None:
        if tiledb_config:
            raise RuntimeError(
                "tiledbsoma is not importable, but tiledb_config was provided."
            )
        return None

    merged = dict(_DEFAULT_S3_CONFIG)
    if tiledb_config:
        merged.update(tiledb_config)

    logger.info("TileDB context: %s", {k: v for k, v in merged.items()})
    return soma.SOMATileDBContext(tiledb_config=merged)


def stream_obs_tables(
    exp,
    value_filter: str,
    column_names: Sequence[str],
) -> Iterator[pa.Table]:
    """
    Stream obs results as Arrow Tables.

    Args:
        exp: Census experiment object.
        value_filter: SOMA value_filter string.
        column_names: Columns to include.

    Yields:
        Arrow Table for each batch of results.
    """
    it = exp.obs.read(value_filter=value_filter, column_names=list(column_names))
    for tbl in it:
        if isinstance(tbl, pa.RecordBatch):
            yield pa.Table.from_batches([tbl])
        elif isinstance(tbl, pa.Table):
            yield tbl
        else:
            yield pa.Table.from_batches(tbl.to_batches())


def write_parquet_stream(
    tables: Iterable[pa.Table],
    outpath: Path,
    compression: str = "zstd",
) -> Tuple[int, int]:
    """
    Write an iterable of Arrow Tables to a single Parquet file.

    Args:
        tables: Iterable of Arrow Tables.
        outpath: Output file path.
        compression: Parquet compression codec.

    Returns:
        Tuple of (rows_written, batches_written).
    """
    writer: Optional[pq.ParquetWriter] = None
    rows = 0
    batches = 0

    try:
        for tbl in tables:
            if tbl.num_rows == 0:
                continue
            if writer is None:
                writer = pq.ParquetWriter(str(outpath), tbl.schema, compression=compression)
            writer.write_table(tbl)
            rows += tbl.num_rows
            batches += 1
            if batches % 10 == 0:
                logger.info(f"Progress: {batches} batches, {rows:,} rows written")
    finally:
        if writer is not None:
            writer.close()

    return rows, batches


def write_parquet_parts(
    tables: Iterable[pa.Table],
    outdir: Path,
    compression: str = "zstd",
    log_every_n: int = 10,
) -> Tuple[int, int]:
    """
    Write an iterable of Arrow Tables to a directory of Parquet part files.

    Each non-empty table is written as a separate file named
    ``part-00000.parquet``, ``part-00001.parquet``, etc.

    Args:
        tables: Iterable of Arrow Tables (e.g. from stream_obs_tables).
        outdir: Target directory for part files. Created if it does not exist.
        compression: Parquet compression codec.
        log_every_n: Log progress every N parts written.

    Returns:
        Tuple of (rows_written, parts_written).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    rows = 0
    parts = 0

    for tbl in tables:
        if tbl.num_rows == 0:
            continue
        part_path = outdir / f"part-{parts:05d}.parquet"
        pq.write_table(tbl, str(part_path), compression=compression)
        rows += tbl.num_rows
        parts += 1
        if log_every_n and parts % log_every_n == 0:
            logger.info(f"Progress: {parts} parts, {rows:,} rows written")

    return rows, parts
