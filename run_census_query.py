"""
Standalone script to run a census_query.
Run: python run_census_query.py

#Disclaimer AI was used to help write this script. Please review carefully.
"""

from CellStrata.census_query import (
    CensusTarget,
    ObsFilters,
    OutputSpec,
    QuerySpec,
    run_query,
)

# --- Configure your query here ---
spec = QuerySpec(
    target=CensusTarget(
        census_version="stable",
        organism="homo_sapiens",
    ),
    obs_filters=ObsFilters(
        is_primary_data=True,
        suspension_type="cell",
        disease_ontology_term_ids=["PATO:0000461"],  # healthy/normal
        tissue_general_labels=["lung"],
        cell_type_labels=["mast cell"],
    ),
    output=OutputSpec(mode="pandas"),
)

# --- Run ---
print("Running query...")
result = run_query(spec)

print(f"\nReturned {len(result):,} rows")
print(f"Columns: {list(result.columns)}")
print(f"\nFirst 5 rows:\n{result.head()}")
print(f"\nCell type breakdown:\n{result['cell_type'].value_counts()}")
