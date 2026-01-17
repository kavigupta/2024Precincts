# Additional Data Sources

This directory contains additional data sources that complement the main state-level data in the `states/` directory.

## Michigan

**CSV Source:** [OpenElections Michigan](https://github.com/openelections/openelections-data-mi/blob/master/2024/20241105__mi__general__precinct.csv)

**Shapefile Source:** [Michigan Open Data Portal](https://gis-michigan.opendata.arcgis.com/maps/Michigan::2024-voting-precincts/about)

**Description:**
- **CSV**: Precinct-level election results from OpenElections project, containing results for U.S. President and other races. Data is in long format with columns: `county`, `precinct`, `office`, `candidate`, `party`, `votes`.
- **Shapefile**: 2024 Voting Precinct boundaries from Michigan's Open Data Portal. Contains 4,340 precincts with columns including `COUNTYFIPS`, `PRECINCT`, `NAME`, `WARD`, etc.

**Files:**
- `20241105__mi__general__precinct.csv` (15 MB) - Election results in long format
- `2024_voting_precincts.zip` (12 MB) - Shapefile containing precinct boundaries

**Data Structure:**
- CSV has ~4,298 unique precincts across 83 counties
- Presidential office is labeled as "President"
- Main candidates: Kamala D. Harris (DEM), Donald J. Trump (REP), plus other candidates
- Shapefile has 4,340 precincts with geometry and metadata

**Note:** The CSV is in long format and needs to be pivoted/aggregated to get presidential results by precinct. The shapefile can be joined with the CSV using `county`/`COUNTYFIPS` and `precinct`/`PRECINCT` identifiers.
