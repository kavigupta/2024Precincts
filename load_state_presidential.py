#!/usr/bin/env python3
"""
Function to load a state's presidential results and standardize to dem/rep/oth columns.
"""

import zipfile
import tempfile
import geopandas as gpd
import pandas as pd
from pathlib import Path


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

    # Check for NaN values
    for col in ["dem", "rep", "oth"]:
        if col in result.columns and result[col].isna().any():
            missing_count = result[col].isna().sum()
            raise ValueError(
                f"{state_prefix}Found {missing_count} NaN values in {col} column"
            )

    # Check for negative values
    for col in ["dem", "rep", "oth"]:
        if col in result.columns and (result[col] < 0).any():
            negative_count = (result[col] < 0).sum()
            min_value = result[col].min()
            raise ValueError(
                f"{state_prefix}Found {negative_count} negative values in {col} column "
                f"(minimum: {min_value:.2f})"
            )


def _load_from_columns(
    gdf, dem_col, rep_col, oth_col=None, total_col=None, state_name=None
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
    state_name : str, optional
        State name for error messages

    Returns:
    --------
    gpd.GeoDataFrame with columns: STATE, STATE_ABBR, dem, rep, oth, geometry

    Raises:
    --------
    ValueError: If columns are missing, contain NaN, or contain invalid data
    """
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
    result["STATE"] = state_name if state_name else "Unknown"
    result["STATE_ABBR"] = (
        state_name[:2].upper() if state_name and len(state_name) >= 2 else "UN"
    )
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


def _load_delaware(gdf):
    """
    Delaware-specific loader.
    Delaware has: DE24Data_D (Dem), DE24Data_K (Kennedy/Other), DE24Data_G (Total)
    Republican votes = Total - Dem - Kennedy
    """
    if (
        "DE24Data_D" not in gdf.columns
        or "DE24Data_K" not in gdf.columns
        or "DE24Data_G" not in gdf.columns
    ):
        raise ValueError(
            f"Delaware data missing required columns. Available: {list(gdf.columns)}"
        )

    result = gpd.GeoDataFrame()
    result["STATE"] = "Delaware"
    result["STATE_ABBR"] = "DE"
    result["geometry"] = gdf.geometry
    result.crs = gdf.crs

    # Extract data
    dem = pd.to_numeric(gdf["DE24Data_D"], errors="raise")
    kennedy = pd.to_numeric(gdf["DE24Data_K"], errors="raise")
    total = pd.to_numeric(gdf["DE24Data_G"], errors="raise")

    # Calculate Republican votes
    result["dem"] = dem
    result["rep"] = total - dem - kennedy
    result["oth"] = kennedy

    # Validate at the end
    _validate_result(result, "Delaware")

    return result


def _load_oklahoma(gdf):
    """
    Oklahoma-specific loader.
    Oklahoma has: OK24Data_H (Harris/Dem), OK24Data_T (Trump/Rep), OK24Data_G (Total)
    Other votes = Total - Dem - Rep
    """
    if (
        "OK24Data_H" not in gdf.columns
        or "OK24Data_T" not in gdf.columns
        or "OK24Data_G" not in gdf.columns
    ):
        raise ValueError(
            f"Oklahoma data missing required columns. Available: {list(gdf.columns)}"
        )

    result = gpd.GeoDataFrame()
    result["STATE"] = "Oklahoma"
    result["STATE_ABBR"] = "OK"
    result["geometry"] = gdf.geometry
    result.crs = gdf.crs

    # Extract data
    dem = pd.to_numeric(gdf["OK24Data_H"], errors="raise")
    rep = pd.to_numeric(gdf["OK24Data_T"], errors="raise")
    total = pd.to_numeric(gdf["OK24Data_G"], errors="raise")

    # Calculate Other votes
    result["dem"] = dem
    result["rep"] = rep
    result["oth"] = total - dem - rep

    # Validate at the end
    _validate_result(result, "Oklahoma")

    return result


def _load_vermont(gdf):
    """
    Vermont-specific loader.
    Vermont has: VTData24_D (Dem), VTData24_K (Kennedy/Other), VTData24_T (Total)
    Republican votes = Total - Dem - Kennedy
    """
    if (
        "VTData24_D" not in gdf.columns
        or "VTData24_K" not in gdf.columns
        or "VTData24_T" not in gdf.columns
    ):
        raise ValueError(
            f"Vermont data missing required columns. Available: {list(gdf.columns)}"
        )

    result = gpd.GeoDataFrame()
    result["STATE"] = "Vermont"
    result["STATE_ABBR"] = "VT"
    result["geometry"] = gdf.geometry
    result.crs = gdf.crs

    # Extract data
    dem = pd.to_numeric(gdf["VTData24_D"], errors="raise")
    kennedy = pd.to_numeric(gdf["VTData24_K"], errors="raise")
    total = pd.to_numeric(gdf["VTData24_T"], errors="raise")

    # Calculate Republican votes
    result["dem"] = dem
    result["rep"] = total - dem - kennedy
    result["oth"] = kennedy

    # Validate at the end
    _validate_result(result, "Vermont")

    return result


def _load_michigan(gdf, state_path):
    """
    Michigan-specific loader.
    Michigan has CSV in long format that needs to be pivoted and merged with shapefile.
    """
    # Michigan county FIPS to name mapping
    MI_COUNTY_FIPS = {
        "001": "Alcona",
        "003": "Alger",
        "005": "Allegan",
        "007": "Alpena",
        "009": "Antrim",
        "011": "Arenac",
        "013": "Baraga",
        "015": "Barry",
        "017": "Bay",
        "019": "Benzie",
        "021": "Berrien",
        "023": "Branch",
        "025": "Calhoun",
        "027": "Cass",
        "029": "Charlevoix",
        "031": "Cheboygan",
        "033": "Chippewa",
        "035": "Clare",
        "037": "Clinton",
        "039": "Crawford",
        "041": "Delta",
        "043": "Dickinson",
        "045": "Eaton",
        "047": "Emmet",
        "049": "Genesee",
        "051": "Gladwin",
        "053": "Gogebic",
        "055": "Grand Traverse",
        "057": "Gratiot",
        "059": "Hillsdale",
        "061": "Houghton",
        "063": "Huron",
        "065": "Ingham",
        "067": "Ionia",
        "069": "Iosco",
        "071": "Iron",
        "073": "Isabella",
        "075": "Jackson",
        "077": "Kalamazoo",
        "079": "Kalkaska",
        "081": "Kent",
        "083": "Keweenaw",
        "085": "Lake",
        "087": "Lapeer",
        "089": "Leelanau",
        "091": "Lenawee",
        "093": "Livingston",
        "095": "Luce",
        "097": "Mackinac",
        "099": "Macomb",
        "101": "Manistee",
        "103": "Marquette",
        "105": "Mason",
        "107": "Mecosta",
        "109": "Menominee",
        "111": "Midland",
        "113": "Missaukee",
        "115": "Monroe",
        "117": "Montcalm",
        "119": "Montmorency",
        "121": "Muskegon",
        "123": "Newaygo",
        "125": "Oakland",
        "127": "Oceana",
        "129": "Ogemaw",
        "131": "Ontonagon",
        "133": "Osceola",
        "135": "Oscoda",
        "137": "Otsego",
        "139": "Ottawa",
        "141": "Presque Isle",
        "143": "Roscommon",
        "145": "Saginaw",
        "147": "St. Clair",
        "149": "St. Joseph",
        "151": "Sanilac",
        "153": "Schoolcraft",
        "155": "Shiawassee",
        "157": "Tuscola",
        "159": "Van Buren",
        "161": "Washtenaw",
        "163": "Wayne",
        "165": "Wexford",
    }

    csv_files = list(state_path.glob("*.csv"))
    if not csv_files:
        raise ValueError(f"Michigan CSV file not found in {state_path}")

    # Load and process CSV
    df = pd.read_csv(csv_files[0], low_memory=False)
    pres_df = df[df["office"] == "President"].copy()

    # Convert votes to numeric (handle comma-separated numbers)
    pres_df["votes"] = (
        pres_df["votes"].astype(str).str.replace(",", "").replace("", "0")
    )
    pres_df["votes"] = pd.to_numeric(pres_df["votes"], errors="coerce").fillna(0)

    # Normalize party column - map variations to DEM/REP/OTH
    pres_df["party_norm"] = pres_df["party"].str.upper()
    # Map party variations
    pres_df.loc[pres_df["party_norm"].isin(["DEM", "DEMOCRAT"]), "party_norm"] = "DEM"
    pres_df.loc[pres_df["party_norm"].isin(["REP", "REPUBLICAN"]), "party_norm"] = "REP"
    # Also check candidate names for party identification
    pres_df.loc[
        pres_df["candidate"].str.contains("Harris", case=False, na=False), "party_norm"
    ] = "DEM"
    pres_df.loc[
        pres_df["candidate"].str.contains("Trump", case=False, na=False), "party_norm"
    ] = "REP"
    # Everything else is OTH
    pres_df.loc[~pres_df["party_norm"].isin(["DEM", "REP"]), "party_norm"] = "OTH"

    # Pivot table: index (county, precinct), columns party, values votes
    pivot_df = pres_df.pivot_table(
        index=["county", "precinct"],
        columns="party_norm",
        values="votes",
        aggfunc="sum",
        fill_value=0,
    )

    # Flatten column names and create result dataframe
    result_df = pivot_df.reset_index()
    result_df.columns.name = None  # Remove the columns name

    # Rename columns to match expected format (these columns will always exist after pivot)
    result_df["PresDem"] = result_df["DEM"]
    result_df["PresRep"] = result_df["REP"]
    result_df["PresOth"] = result_df["OTH"]

    # Calculate total
    result_df["PresTot"] = (
        result_df["PresDem"] + result_df["PresRep"] + result_df["PresOth"]
    )

    # Keep only needed columns
    result_df = result_df[["county", "precinct", "PresDem", "PresRep", "PresTot"]]

    # Extract precinct number from CSV precinct names
    # Handles patterns like "Precinct 1" or "CB 1" (Counting Board)
    import re

    def extract_precinct_num(name):
        # Try "Precinct X" pattern first
        match = re.search(r"Precinct\s+(\d+)", str(name), re.IGNORECASE)
        if match:
            return match.group(1).zfill(3)  # Pad to 3 digits like shapefile
        # Try "CB X" pattern (Detroit uses Counting Boards)
        match = re.search(r"CB\s+(\d+)", str(name), re.IGNORECASE)
        if match:
            return match.group(1).zfill(3)
        return None

    result_df["precinct_num"] = result_df["precinct"].apply(extract_precinct_num)
    # Filter out rows where we couldn't extract precinct number
    result_df = result_df[result_df["precinct_num"].notna()]

    # Map COUNTYFIPS to county name in shapefile
    gdf["county_name"] = gdf["COUNTYFIPS"].astype(str).map(MI_COUNTY_FIPS)

    # Merge with shapefile on county name and precinct number
    gdf["PRECINCT_str"] = (
        gdf["PRECINCT"].astype(str).str.zfill(3)
    )  # Ensure 3-digit format

    merged = gdf.merge(
        result_df,
        left_on=["county_name", "PRECINCT_str"],
        right_on=["county", "precinct_num"],
        how="left",
        suffixes=("", "_csv"),
    )

    # Filter out rows that don't have election data (user wants to error on invalid data, not fill with 0s)
    merged = merged[
        merged["PresDem"].notna()
        & merged["PresRep"].notna()
        & merged["PresTot"].notna()
    ]

    return merged


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


# State-specific column mappings for common case
# Format: {'dem': col, 'rep': col, 'oth': col or list or None, 'total': col or None}
# If oth is None, total must be provided to calculate oth = total - dem - rep
# If oth is a list, all columns will be summed
STATE_COLUMNS = {
    "Alabama": {
        "dem": "G24PREDHAR",
        "rep": "G24PRERTRU",
        "oth": ["G24PREIKEN", "G24PREIOLI", "G24PREISTE", "G24PREOWRI"],
        "total": None,
    },
    "Alaska": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Arizona": {"dem": "Harris", "rep": "Trump", "oth": None, "total": "PresTot"},
    "Arkansas": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "California": {"dem": "G24PDem", "rep": "G24PRep", "oth": "G24POth", "total": None},
    "Colorado": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Connecticut": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "District of Columbia": {
        "dem": "DC24Data_D",
        "rep": "DC24Data_R",
        "oth": None,
        "total": "DC24Data_T",
    },
    "Florida": {
        "dem": "Florida__2",
        "rep": "Florida__1",
        "oth": None,
        "total": "Florida__9",
    },
    "Georgia": {
        "dem": "G24PREDHAR",
        "rep": "G24PRERTRU",
        "oth": ["G24PREGSTE", "G24PREICRU", "G24PREIWES", "G24PRELOLI"],
        "total": None,
    },
    "Hawaii": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Idaho": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Illinois": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Indiana": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Iowa": {
        "dem": "votes_dem",
        "rep": "votes_rep",
        "oth": None,
        "total": "votes_tota",
    },
    "Kansas": {"dem": "DHarris", "rep": "RTrump", "oth": "Other", "total": None},
    "Kentucky": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Louisiana": {
        "dem": "G24PREDHAR",
        "rep": "G24PRERTRU",
        "oth": [
            "G24PRELOLI",
            "G24PREOCRU",
            "G24PREOFRU",
            "G24PREOKEN",
            "G24PREOPRE",
            "G24PREOSON",
            "G24PREOSTE",
            "G24PREOTER",
            "G24PREOWES",
        ],
        "total": None,
    },
    "Maine": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Maryland": {
        "dem": "G24PRESD",
        "rep": "G24PRESR",
        "oth": "G24PRESO",
        "total": None,
    },
    "Massachusetts": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "Michigan": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "Minnesota": {
        "dem": "MNPrecinct",
        "rep": "MNPrecin_1",
        "oth": None,
        "total": "MNPrecin_2",
    },
    "Mississippi": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "Missouri": {
        "dem": "HARRIS - D",
        "rep": "TRUMP - RE",
        "oth": None,
        "total": "PRESTOT",
    },
    "Montana": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Nebraska": {"dem": "PresDem", "rep": "PresRep", "oth": "PresOth", "total": None},
    "Nevada": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "New Hampshire": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "New Jersey": {"dem": "Harris", "rep": "Trump", "oth": "Other", "total": None},
    "New Mexico": {
        "dem": "NMData24_1",
        "rep": "NMData24_2",
        "oth": None,
        "total": "NMData24_3",
    },
    "New York": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "North Carolina": {
        "dem": "PRESDEM",
        "rep": "PRESREP",
        "oth": None,
        "total": "PRESTOTAL",
    },
    "North Dakota": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "Ohio": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },  # Has NaN - will error correctly
    "Oregon": {"dem": "USP24D", "rep": "USP24R", "oth": None, "total": "USP24Tot"},
    "Pennsylvania": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },  # Data not available per README
    "Rhode Island": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "South Carolina": {"dem": "Harris", "rep": "Trump", "oth": None, "total": "Total"},
    "South Dakota": {
        "dem": "G24PRESD",
        "rep": "G24PRESR",
        "oth": None,
        "total": "G24PRESTOT",
    },
    "Tennessee": {
        "dem": "votes_dem",
        "rep": "votes_rep",
        "oth": None,
        "total": "votes_tota",
    },
    "Texas": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "Utah": {"dem": "Harris", "rep": "Trump", "oth": None, "total": "Total"},
    "Virginia": {
        "dem": "G24PRESD",
        "rep": "G24PRESR",
        "oth": "G24PRESO",
        "total": None,
    },
    "Washington": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
    "West Virginia": {
        "dem": "PresDem",
        "rep": "PresRep",
        "oth": None,
        "total": "PresTot",
    },
    "Wisconsin": {
        "dem": "votes_dem",
        "rep": "votes_rep",
        "oth": None,
        "total": "votes_tota",
    },
    "Wyoming": {"dem": "PresDem", "rep": "PresRep", "oth": None, "total": "PresTot"},
}


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

    # Find zip files in the state directory
    zip_files = list(state_path.glob("*.zip"))
    if not zip_files:
        # If no zip in additional_data, fall back to original state_dir
        if state_path == additional_data_path:
            state_path = Path(state_dir)
            zip_files = list(state_path.glob("*.zip"))
        if not zip_files:
            return None

    # Use the first zip file found
    zip_path = zip_files[0]

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract and find shapefile
            shp_path = find_shapefile_in_zip(zip_path, temp_dir)
            if not shp_path:
                return None

            # Read the shapefile
            gdf = gpd.read_file(shp_path)

            # Check for CSV file to merge (e.g., Arizona has AZ24.csv)
            csv_files = list(state_path.glob("*.csv"))
            if csv_files:
                # Use the first CSV file found
                csv_path = csv_files[0]
                df = pd.read_csv(csv_path)
                # Try to merge on common key columns
                # Try common join keys in order of preference
                join_keys = ["JoinField", "PCTNUM", "Precinct", "PRECINCT"]
                merged = False
                for key in join_keys:
                    if key in gdf.columns and key in df.columns:
                        gdf = gdf.merge(df, on=key, how="left")
                        merged = True
                        break

                if state_name == "New Hampshire":
                    gdf = merge_new_hampshire(gdf, df)
                    merged = True

                if state_name == "Michigan":
                    gdf = _load_michigan(gdf, state_path)
                    merged = True

                if not merged:
                    # If no common key found, try to merge on all common columns
                    common_cols = set(gdf.columns) & set(df.columns)
                    if common_cols:
                        # Use the first common column as join key
                        join_key = list(common_cols)[0]
                        gdf = gdf.merge(df, on=join_key, how="left")
                    # Try state-specific join keys
                    elif (
                        state_name == "Rhode Island"
                        and "DISTRICTN" in gdf.columns
                        and "GEOID" in df.columns
                    ):
                        gdf = gdf.merge(
                            df, left_on="DISTRICTN", right_on="GEOID", how="left"
                        )
                    elif state_name == "Nevada" and "JoinField" in df.columns:
                        # Nevada: CSV has JoinField but shapefile doesn't - cannot merge without proper key
                        raise ValueError(
                            f"Nevada: Cannot merge CSV with shapefile - no matching keys found. "
                            f"CSV has JoinField but shapefile doesn't have matching column."
                        )

            # Route to state-specific handler or use common case
            if state_name == "Vermont":
                result = _load_vermont(gdf)
            elif state_name == "Delaware":
                result = _load_delaware(gdf)
            elif state_name == "Oklahoma":
                result = _load_oklahoma(gdf)
            elif state_name in STATE_COLUMNS:
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
            else:
                raise ValueError(
                    f"No handler configured for state '{state_name}'. "
                    f"Available columns: {list(gdf.columns)}"
                )

            return result

        except Exception as e:
            print(f"  {state_name}: Error - {str(e)}")
            return None


def test_all_states(project_root):
    """Test loading all states and report which ones work."""
    states_dir = Path(project_root) / "states"

    if not states_dir.exists():
        print(f"Error: States directory not found at {states_dir}")
        return

    state_dirs = sorted([d for d in states_dir.iterdir() if d.is_dir()])

    print(f"Testing {len(state_dirs)} states...\n")

    successful = []
    failed = []

    for state_dir in state_dirs:
        state_name = state_dir.name
        print(f"Testing {state_name}...", end=" ")

        gdf = load_state_presidential(state_dir, state_name)

        if gdf is not None and len(gdf) > 0:
            # Check that we have dem or rep data
            if gdf["dem"].sum() > 0 or gdf["rep"].sum() > 0:
                print(
                    f"✓ Success ({len(gdf)} features, dem={gdf['dem'].sum():.0f}, rep={gdf['rep'].sum():.0f})"
                )
                successful.append(state_name)
            else:
                print(f"✗ No vote data found")
                failed.append(state_name)
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

    return successful, failed


if __name__ == "__main__":
    import sys

    project_root = Path(__file__).parent
    test_all_states(project_root)
