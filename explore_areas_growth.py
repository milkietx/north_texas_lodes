import sys
sys.path.append(r'C:\Users\cmg0530\Code Library\lodes-dataset-construction\main')
from analysis import *
sys.path.append(r'C:\Users\cmg0530\Code Library\gis_suite')
from census_pulls import *
import pandas as pd
import os


spath = r"C:\Users\cmg0530\Projects\lodes_package\lodes_tx_slim.db"

con,cur = connect_to_od(spath=spath)

#get time of day for dfw
total_employment = census_pull(tables=["C24020"],geom="tract",year="2022")

#dfw counties
dfw_counties = ['085','113','121','139','143',
                '221','231','251','257','349',
                '363','367','397','425','439','497']

dfw_tdata = total_employment.query("county in @dfw_counties")

#get crosswalked blocks
cw = pd.read_sql("select tabblk2020, trct,trctname from tx_xwalk",con=con)
dfw_cw_blks = cw.merge(dfw_tdata,left_on='trct',right_on=['GEOID'])[['tabblk2020','GEOID']]
cw_mapdict = dict(zip(dfw_cw_blks['tabblk2020'],dfw_cw_blks['GEOID']))
#get od data for the dfw_cw_blks in 2022 and 2012

dfs = []
for y in ['2012','2017','2022']:
    q = generate_query(data_type='od',
                       perspective='home',
                       job_type='primary',
                       year=y,
                       geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
    od_df = pull_data(query=q,crsr=cur,rename=True)
    retype(od_df)
    od_df = od_df.fillna(0)
    dfs.append(od_df)

od_primary = pd.concat(dfs)
od_primary['dest_tract'] = od_primary['w_geocode'].map(cw_mapdict)
od_primary['orig_tract'] = od_primary['h_geocode'].map(cw_mapdict)


od_set = od_primary.groupby(['year',
                             'orig_tract',
                             'dest_tract']).agg({'total':'sum',
                                                 'goods':'sum',
                                                 'transp':'sum',
                                                 'other':'sum',}).reset_index()

#get the percentage of total commuters for each OD pair year
ods_set_list = []
for each_year in od_set['year'].unique():
    ods = od_set[od_set['year'] == each_year]
    tm = ods.groupby('orig_tract').agg({'total':'sum'}).reset_index()
    tm_dict = dict(zip(tm['orig_tract'],tm['total']))
    ods['total_origin_lodes'] = ods['orig_tract'].map(tm_dict)
    ods['proportion_of_origin_to_dest'] = ods.apply(lambda x: x['total']/tm_dict[x['orig_tract']],axis=1).fillna(0)
    ods_set_list.append(ods)

#pair with each for value
od_set = pd.concat(ods_set_list)

#pivot by year
od_pivot = od_set.pivot(index=['orig_tract','dest_tract'],
             columns=['year'],
             values=['total','goods','transp','other',
                     'total_origin_lodes','proportion_of_origin_to_dest'])

rev_columns = [f"{x[0]}_{x[1]}" for x in od_pivot.columns]
od_pivot.columns = rev_columns
od_pivot = od_pivot.reset_index()

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


