#!/usr/bin/env python3
"""
Function to load a state's presidential results and standardize to dem/rep/oth columns.
"""

import os
import gzip
import zipfile
import tempfile
import geopandas as gpd
import pandas as pd
from pathlib import Path

from columns import STATE_COLUMNS


def find_shapefile_in_zip(zip_path, temp_dir):
    """Extract and find the shapefile in a zip archive."""
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(temp_dir)

    # Find .shp file in the extracted directory (search recursively)
    shp_files = list(Path(temp_dir).rglob("*.shp"))
    if not shp_files:
        return None
    return str(shp_files[0])


def _validate_result(result, state_name=None):
    """
    Validate the result GeoDataFrame for NaN and negative values.

    Parameters:
    -----------
    result : gpd.GeoDataFrame
        The result GeoDataFrame with dem, rep, oth columns
    state_name : str, optional
        State name for error messages

    Raises:
    --------
    ValueError: If any column contains NaN or negative values
    """
    state_prefix = f"{state_name}: " if state_name else ""


    columns = ["dem", "rep", "oth"]

    total = result[columns].sum().sum()
    num_rows = result.shape[0]

    # Check for NaN values
    for col in columns:
        if col in result.columns and result[col].isna().any():
            missing_count = result[col].isna().sum()
            if missing_count > 0.01 * num_rows:
                raise ValueError(
                    f"{state_prefix}Found {missing_count} NaN values in {col} column"
                )
            # small enough to ignore
            result[col][result[col].isna()] = 0

    # Check for negative values
    for col in columns:
        if col in result.columns and (result[col] < 0).any():
            count_negative = (result[col] < 0).sum()
            total_negative = result[col][result[col] < 0].sum()
            if total_negative > 0.0001 * total:
                raise ValueError(
                    f"{state_prefix}Found {count_negative} negative values in {col} column "
                    f"(total: {total_negative:.2f})"
                )
            # small enough to ignore
            result[col][result[col] < 0] = 0


def _load_from_columns(
    gdf, dem_col, rep_col, oth_col=None, total_col=None, *, state_name
):
    """
    Common case: Extract presidential votes from explicitly named columns.

    Parameters:
    -----------
    gdf : gpd.GeoDataFrame
        The GeoDataFrame with election data
    dem_col : str
        Column name for Democratic votes
    rep_col : str
        Column name for Republican votes
    oth_col : str, list, or None
        Column name(s) for Other votes. If a list, will sum all columns.
        If None, total_col must be provided to calculate oth = total - dem - rep
    total_col : str, optional
        Column name for Total votes (required if oth_col is None)
    state_name : str
        State name for error messages

    Returns:
    --------
    gpd.GeoDataFrame with columns: STATE, STATE_ABBR, dem, rep, oth, geometry

    Raises:
    --------
    ValueError: If columns are missing, contain NaN, or contain invalid data
    """
    assert state_name is not None
    # Validate that we have either oth_col or total_col
    if oth_col is None and total_col is None:
        raise ValueError(
            f"Must provide either 'oth_col' or 'total_col'. "
            f"Available columns: {list(gdf.columns)}"
        )

    # Validate columns exist
    missing_cols = []
    if dem_col not in gdf.columns:
        missing_cols.append(dem_col)
    if rep_col not in gdf.columns:
        missing_cols.append(rep_col)

    # Handle oth_col - can be string, list, or None
    if oth_col is not None:
        if isinstance(oth_col, str):
            oth_cols_to_check = [oth_col]
        else:  # list
            oth_cols_to_check = oth_col
        for col in oth_cols_to_check:
            if col not in gdf.columns:
                missing_cols.append(col)

    if total_col and total_col not in gdf.columns:
        missing_cols.append(total_col)

    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}. "
            f"Available columns: {list(gdf.columns)}"
        )

    # Create result GeoDataFrame
    result = gpd.GeoDataFrame()
    result["geometry"] = gdf.geometry
    result.crs = gdf.crs

    # Extract Democratic votes
    result["dem"] = pd.to_numeric(gdf[dem_col], errors="raise")

    # Extract Republican votes
    result["rep"] = pd.to_numeric(gdf[rep_col], errors="raise")

    # Extract or calculate Other votes
    if oth_col is not None:
        if isinstance(oth_col, str):
            # Single column
            result["oth"] = pd.to_numeric(gdf[oth_col], errors="raise")
        else:
            # List of columns - sum them
            result["oth"] = pd.Series(0, index=gdf.index, dtype="float64")
            for col in oth_col:
                col_data = pd.to_numeric(gdf[col], errors="raise")
                result["oth"] += col_data
    else:
        # Calculate oth from total
        total = pd.to_numeric(gdf[total_col], errors="raise")
        result["oth"] = total - result["dem"] - result["rep"]

    # Validate at the end
    _validate_result(result, state_name)

    return result


def merge_new_hampshire(gdf, df):
    """
    Merge New Hampshire shapefile with CSV data.
    """
    gdf["ident"] = gdf["COUNTY"] + "::" + gdf["WARD"]
    size = len(gdf)
    gdf = gdf.dissolve(by="ident")
    # After dissolve, 'ident' becomes the index, so reset it to a column
    gdf = gdf.reset_index()
    assert size - 1 == len(gdf), "Dissolve should have reduced the number of rows by 1"
    df["ident"] = df["County"] + "::" + df["Precinct"]
    return gdf.merge(df, on="ident", how="left")


def load_gdf(state_name, state_path):
    # Special handling for Pennsylvania: load GeoJSON directly (it has all the data)
    if state_name == "Pennsylvania" or state_name == "Michigan":
        [geojson_file] = list(state_path.glob("*.geojson.gz"))
        with gzip.open(geojson_file, "rb") as f:
            return gpd.read_file(f)

    # Find zip files in the state directory
    [zip_file] = list(state_path.glob("*.zip"))
    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        # Extract and find shapefile
        shp_path = find_shapefile_in_zip(zip_file, temp_dir)
        if not shp_path:
            return None

        # Read the shapefile
        return gpd.read_file(shp_path)

def load_table_with_possible_csv(state_name, state_path):
    gdf = load_gdf(state_name, state_path)
    if gdf is None:
        return None

    # Check for CSV file to merge (e.g., Arizona has AZ24.csv)
    csv_files = list(state_path.glob("*.csv"))
    if not csv_files:
        return gdf
    # Use the first CSV file found
    csv_path = csv_files[0]
    df = pd.read_csv(csv_path)
    # Try to merge on common key columns
    # Try common join keys in order of preference

    if state_name == "New Hampshire":
        return merge_new_hampshire(gdf, df)

    if state_name == "Maine":
        df["JoinField"] = df["COUNTY"] + "::" + df["TOWN"]
        gdf["JoinField"] = gdf["COUNTY"] + "::" + gdf["TOWN"]
    if state_name == "Rhode Island":
        gdf["JoinField"] = gdf["DISTRICTN"]
        df["JoinField"] = df["GEOID"]
    join_keys = ["JoinField", "PCTNUM", "Precinct", "PRECINCT"]
    for key in join_keys:
        if key in gdf.columns and key in df.columns:
            in_gdf = set(gdf[key])
            in_df = set(df[key])
            if state_name == "Utah":
                # these have 0 votes anyway
                in_df -= {'CACHE: LOG24:U', 'CACHE: SMI01:U2'}
            if in_gdf == in_df:
                return gdf.merge(df, on=key, how="left")
            else:
                extra_in_gdf = in_gdf - in_df
                extra_in_df = in_df - in_gdf
                raise ValueError(f"Common key {key} has different values in gdf and df: {extra_in_gdf} in gdf but not in df, {extra_in_df} in df but not in gdf")


    raise ValueError(f"No common key found for {state_name}")


def load_state_presidential(state_dir, state_name):
    """
    Load a state's presidential results and standardize to dem/rep/oth columns.

    Parameters:
    -----------
    state_dir : str or Path
        Path to the state directory (will check additional_data/ first if available)
    state_name : str
        Name of the state

    Returns:
    --------
    gpd.GeoDataFrame with columns: STATE, STATE_ABBR, dem, rep, oth, geometry
    Returns None if loading fails.
    """
    # Check additional_data/ first, then fall back to provided state_dir
    project_root = (
        Path(state_dir).parent.parent
        if Path(state_dir).parent.name == "states"
        else Path(state_dir).parent
    )
    additional_data_path = project_root / "additional_data" / state_name

    if additional_data_path.exists() and additional_data_path.is_dir():
        state_path = additional_data_path
    else:
        state_path = Path(state_dir)

    gdf = load_table_with_possible_csv(state_name, state_path)

    # Route to state-specific handler or use common case
    assert state_name in STATE_COLUMNS, f"State {state_name} not found in STATE_COLUMNS"

    cols = STATE_COLUMNS[state_name]
    # Validate that if oth is None, total must be provided
    if cols["oth"] is None and cols["total"] is None:
        raise ValueError(
            f"State '{state_name}' configuration error: "
            f"must provide either 'oth' or 'total' column. "
            f"Available columns: {list(gdf.columns)}"
        )
    result = _load_from_columns(
        gdf,
        dem_col=cols["dem"],
        rep_col=cols["rep"],
        oth_col=cols["oth"],
        total_col=cols["total"],
        state_name=state_name,
    )

    return result


def load_all_states(project_root):
    """Test loading all states and report which ones work."""
    states_dir = Path(project_root) / "states"

    if not states_dir.exists():
        print(f"Error: States directory not found at {states_dir}")
        return

    state_dirs = sorted([d for d in states_dir.iterdir() if d.is_dir()])

    print(f"Testing {len(state_dirs)} states...\n")

    successful = {}
    failed = {}

    for state_dir in state_dirs:
        state_name = state_dir.name
        print(f"Testing {state_name}...", end=" ")

        try:
            gdf = load_state_presidential(state_dir, state_name)
        except Exception as e:
            print(f"✗ Failed to load: {e}")
            failed[state_name] = e
            continue

        if gdf is not None:
            print(
                f"✓ Success ({len(gdf)} features, dem={gdf['dem'].sum():.0f}, rep={gdf['rep'].sum():.0f})"
            )
            successful[state_name] = gdf
        else:
            print(f"✗ Failed to load")
            failed.append(state_name)

    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Successful: {len(successful)}/{len(state_dirs)}")
    print(f"  Failed: {len(failed)}/{len(state_dirs)}")

    if successful:
        print(f"\nSuccessful states:")
        for state in successful:
            print(f"  - {state}")

    if failed:
        print(f"\nFailed states:")
        for state in failed:
            print(f"  - {state}")
            print(f"    - {failed[state]}")

    return successful, failed

def merge_shapefiles(project_root):
    successful, failed = load_all_states(project_root)
    assert not failed, "Failed to load some states"
    frames = []
    for state_name, gdf in successful.items():
        gdf["state_name"] = state_name
        gdf = gdf.to_crs(4326)
        gdf = gdf.reset_index(drop=True)
        gdf["precinct_no"] = gdf.index
        frames.append(gdf)
    gdf = pd.concat(frames)
    return gdf

if __name__ == "__main__":
    project_root = Path(__file__).parent
    gdf = merge_shapefiles(project_root)
    os.makedirs("output", exist_ok=True)
    gdf.to_file("output/all.shp", driver="ESRI Shapefile")
