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

q = generate_query(data_type='wac',
                       job_type='primary',
                       year='2023',
                       geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
hc_df = pull_data(query=q,crsr=cur,rename=True)


b20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geocode from blocks_2020_geom",con=con)
b20 = gpd.GeoDataFrame(b20,geometry=gpd.GeoSeries.from_wkt(b20['geom_wkt']),crs='epsg:4326')
b20['county'] = b20['geocode'].str[2:5]
dfw20 = b20.query("county in @dfw_counties")
hc_gdf = dfw20.merge(hc_df,left_on='geocode',right_on='w_geocode',how='left').fillna(0)
hc_gdf.drop(columns='geom_wkt',inplace=True)

#define majority healthcare tracts plus adjacent tracts with over 100+ employees
#spatial autocorrelation
def major_industry(gdf=hc_gdf,industry='tot'):
    import libpysal as lp
    import esda
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
    df.assign(cl=labels).plot(
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



    gdf['Health_62'].describe()

#adjusted destiantion for acs
adj_dict = dict(zip(dfw_tdata['GEOID'],dfw_tdata['C24020_001E']))
od_pivot['origin_total_acs2024'] = od_pivot['orig_tract'].map(adj_dict)

#example of utsw analysis
utsw = ['48113010001','48113000409','48113000401']

med_dist_commute = od_pivot.query("dest_tract in @utsw")
med_dist_commute['expected_from_2012_proportions_val2024'] = med_dist_commute['proportion_of_origin_to_dest_2012'].fillna(0) * med_dist_commute['origin_total_acs2024'] 
med_dist_commute['expected_from_2022_proportions_val2024'] = med_dist_commute['proportion_of_origin_to_dest_2022'].fillna(0) * med_dist_commute['origin_total_acs2024'] 
med_dist_commute['difference_2012_2022_val2024'] = med_dist_commute['expected_from_2022_proportions_val2024'] - med_dist_commute['expected_from_2012_proportions_val2024']
med_dist_gb = med_dist_commute.groupby('orig_tract').agg({'total_2012':'sum',
                                            'total_2022':'sum',
                                            'expected_from_2012_proportions_val2024':'sum',
                                            'expected_from_2022_proportions_val2024':'sum',
                                            'difference_2012_2022_val2024':'sum'
                                            })

t20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geoid from tracts_2020_geom",con=con)
t20 = gpd.GeoDataFrame(t20,geometry=gpd.GeoSeries.from_wkt(t20['geom_wkt']),crs='epsg:4326')
t20.drop(columns='geom_wkt',inplace=True)

sdata = t20.merge(med_dist_gb,left_on='GEOID',right_on='orig_tract')
sdata = sdata.to_crs("EPSG:3857")


import mapclassify
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import contextily as ctx

def quickplot():
    fig, ax = plt.subplots(1, 1, figsize=(19, 19))
    print("Plotting...")
    sdata.query("difference_2012_2022_val2024 > -25").plot(
        ax=ax,
        column="difference_2012_2022_val2024",
        cmap="magma",
        alpha=0.8,
        linewidth=2,
        zorder=4,
        legend=True,
        legend_kwds={"shrink": 0.3} ,)
    
    sdata.query('GEOID in @utsw').dissolve().plot(
        ax=ax,
        cmap="viridis",
        alpha=0.8,
        linewidth=1,
        facecolor=None,
        edgecolor='#FF0000',
        zorder=7,
        legend=True,
        legend_kwds={"shrink": 0.3} ,)


    
    #xmin, ymin, xmax, ymax = sdata.geometry.total_bounds
    #farther
    #xmin, ymin, xmax, ymax = [-10869027.941352,3821862.042989,-10711378.159285,3932572.978240]
    #close
    xmin, ymin, xmax, ymax = [-10847756.368917,3834478.418943,-10737203.292594,3918201.362709]
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis("off")
    fig.tight_layout()
    ctx.add_basemap(ax=ax, source=ctx.providers.CartoDB.PositronOnlyLabels, zorder=8)
    ctx.add_basemap(ax=ax, source=ctx.providers.CartoDB.PositronNoLabels, zorder=3)
    plt.show()


quickplot(orig_geoid=f"48113013500")


