import geopandas as gpd

margins_from_county_level = {'California': 0.20138532703724044,
 'Texas': -0.1369360827514333,
 'Florida': -0.13100296900541208,
 'New York': 0.12689486117138934,
 'Pennsylvania': -0.0170973099167127,
 'Illinois': 0.1097538645525473,
 'Ohio': -0.11212750977143691,
 'Georgia': -0.021923613255271805,
 'North Carolina': -0.03222876351294367,
 'Michigan': -0.014146215172651532,
 'New Jersey': 0.059095308029419154,
 'Virginia': 0.058071519590521155,
 'Washington': 0.18336913462611262,
 'Arizona': -0.055286032385856845,
 'Massachusetts': 0.2534393996473657,
 'Tennessee': -0.29719883731480545,
 'Indiana': -0.18977536898796735,
 'Maryland': 0.2875091605458192,
 'Missouri': -0.18418884846184982,
 'Wisconsin': -0.00860766224537093,
 'Colorado': 0.10991286710564786,
 'Minnesota': 0.042564200985153566,
 'South Carolina': -0.17867738821257853,
 'Alabama': -0.3058937612571087,
 'Louisiana': -0.2200500753621743,
 'Kentucky': -0.30552657611576456,
 'Oregon': 0.143916085943417,
 'Oklahoma': -0.3426275385924797,
 'Connecticut': 0.14509257979708015,
 'Utah': -0.21590227097532097,
 'Iowa': -0.1326258458073125,
 'Nevada': -0.030985156649874713,
 'Arkansas': -0.30636962278764407,
 'Mississippi': -0.2288877596888618,
 'Kansas': -0.16115580777513544,
 'New Mexico': 0.060007385724326144,
 'Nebraska': -0.2056898577746712,
 'Idaho': -0.365019473658616,
 'West Virginia': -0.4187449992785839,
 'Hawaii': 0.2310636922582676,
 'New Hampshire': 0.027875590066569082,
 'Maine': 0.06937302661253938,
 'Rhode Island': 0.13854646642867327,
 'Montana': -0.19930410323684836,
 'Delaware': 0.1473665079138628,
 'South Dakota': -0.2919458549573109,
 'North Dakota': -0.36755154646235166,
 'Alaska': -0.13797412503046885,
 'District of Columbia': 0.8380944489963733,
 'Vermont': 0.31768420995171714,
 'Wyoming': -0.46219115234294167}

margins_to_check_against = margins_from_county_level.copy()
# manually from wiki
margins_to_check_against["Alabama"] = -30.47e-2
margins_to_check_against["Alaska"] = -13.13e-2
margins_to_check_against["Delaware"] = 14.70e-2
margins_to_check_against["Illinois"] = 10.90e-2
margins_to_check_against["Maryland"] = 28.54e-2
margins_to_check_against["Massachusetts"] = 25.20e-2
margins_to_check_against["New York"] = 12.60e-2
margins_to_check_against["North Dakota"] = -36.45e-2
margins_to_check_against["Rhode Island"] = 13.78e-2
margins_to_check_against["Tennessee"] = -29.72e-2
margins_to_check_against["Vermont"] = 31.51e-2
margins_to_check_against["Wyoming"] = -45.76e-2

def load_precincts():
    gdf = gpd.read_file("output/all.shp")
    return gdf

def check_margins():
    gdf = load_precincts()
    gdf = gdf[["state_name", "dem", "rep", "oth"]].groupby("state_name").sum()
    gdf["total"] = gdf.sum(axis=1)
    gdf = (gdf.T / gdf["total"]).T
    margin = gdf["dem"] - gdf["rep"]
    for state_name in gdf.index:
        diff = margin[state_name] - margins_to_check_against[state_name]
        if abs(diff) > 0.001:
            print(f"{state_name:20s}: {margin[state_name]:10.2%} {margins_to_check_against[state_name]:10.2%} {diff:10.2%}")

if __name__ == "__main__":
    check_margins()