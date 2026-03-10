import sys
sys.path.append(r'C:\Users\cmg0530\Code Library\lodes-dataset-construction\main')
from analysis import *
sys.path.append(r'C:\Users\cmg0530\Code Library\gis_suite')
from census_pulls import *
import pandas as pd
import geopandas as gpd
import os


spath = r"C:\Users\cmg0530\Projects\lodes_package\lodes_tx_slim.db"

con,cur = connect_to_od(spath=spath)

#get time of day for dfw
total_employment = census_pull(tables=["C24020"],geom="tract",year="2023")

#dfw counties
dfw_counties = ['085','113','121','139','143',
                '221','231','251','257','349',
                '363','367','397','425','439','497']

dfw_tdata = total_employment.query("county in @dfw_counties")

#get crosswalked blocks
cw = pd.read_sql("select tabblk2020, trct,trctname from tx_xwalk",con=con)
dfw_cw_blks = cw.merge(dfw_tdata,left_on='trct',right_on=['GEOID'])[['tabblk2020','GEOID']]
cw_mapdict = dict(zip(dfw_cw_blks['tabblk2020'],dfw_cw_blks['GEOID']))

#get healthcare concentrations

q = generate_query(data_type='rac',
                       job_type='primary',
                       year='2023',
                       geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
rac_df = pull_data(query=q,crsr=cur,rename=True)

q = generate_query(data_type='wac',
                       job_type='primary',
                       year='2023',
                       geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
wac_df = pull_data(query=q,crsr=cur,rename=True)


b20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geocode from blocks_2020_geom",con=con)
b20 = gpd.GeoDataFrame(b20,geometry=gpd.GeoSeries.from_wkt(b20['geom_wkt']),crs='epsg:4326')
b20['county'] = b20['geocode'].str[2:5]
dfw20 = b20.query("county in @dfw_counties")

hc_gdf = dfw20.merge(hc_df,left_on='geocode',right_on='w_geocode',how='left').fillna(0)
hc_gdf.drop(columns='geom_wkt',inplace=True)

#define majority healthcare tracts plus adjacent tracts with over 100+ employees
#spatial autocorrelation
def major_industry(gdf=hc_gdf,industry='Health_62'):
    import libpysal as lp
    import esda
    import matplotlib.pyplot  as plt
    wq =  lp.weights.Queen.from_dataframe(gdf)
    wq.transform = 'r'
    y = gdf[industry]
    ylag = lp.weights.lag_spatial(wq, y)
    
    li = esda.moran.Moran_Local(y, wq)
    (li.p_sim < 0.05).sum()

    sig = li.p_sim < 0.05
    hotspot = sig * li.q == 1
    coldspot = sig * li.q == 3
    doughnut = sig * li.q == 2
    diamond = sig * li.q == 4

    spots = ["n.sig.", "hot spot"]
    labels = [spots[i] for i in hotspot * 1]

    df = gdf
    from matplotlib import colors

    hmap = colors.ListedColormap(["red", "lightgrey"])
    f, ax = plt.subplots(1, figsize=(9, 9))
    df = df.to_crs("EPSG:3857")
    df = df.assign(cl=labels)
    df.plot(
        column="cl",
        categorical=True,
        k=2,
        cmap=hmap,
        linewidth=0.1,
        ax=ax,
        edgecolor="white",
        legend=True,
    )
    xmin, ymin, xmax, ymax = [-10847756.368917,3834478.418943,-10737203.292594,3918201.362709]
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")
    plt.show()

    df.to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
               layer='HealthCare_Zones_v1')


    gdf['Health_62'].describe()


