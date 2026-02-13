from __future__ import annotations
import pandas as pd
import cellxgene_census as cxc

from .config import CensusQueryConfig
from .filters import build_value_filter


def validate_obs_cols(census, organism: str, obs_cols: list[str]) -> None:
    human = census["census_data"][organism]
    schema_names = set(human.obs.keys())
    missing = [c for c in obs_cols if c not in schema_names]
    if missing:
        raise ValueError(f"Missing obs columns in schema: {missing}")


def fetch_obs_df(cfg: CensusQueryConfig) -> tuple[pd.DataFrame, str]:
    value_filter = build_value_filter(cfg.obs_filters)

    with cxc.open_soma(census_version=cfg.target.census_version) as census:
        validate_obs_cols(census, cfg.target.organism, cfg.obs_cols)

        human = census["census_data"][cfg.target.organism]
        df = (
            human.obs.read(value_filter=value_filter, column_names=cfg.obs_cols)
            .concat()
            .to_pandas()
        )

    return df, value_filter


def summarize_by_dataset(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("dataset_id", sort=True)
        .size()
        .reset_index(name="cell_count")
    )
