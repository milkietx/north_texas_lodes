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
time_of_day = census_pull(tables=["B08302"],geom="tract",year="2022")

#dfw counties
dfw_counties = ['085','113','121','139','143',
                '221','231','251','257','349',
                '363','367','397','425','439','497']

dfw_tdata = time_of_day.query("county in @dfw_counties")

#get crosswalked blocks
cw = pd.read_sql("select tabblk2020, trct,trctname from tx_xwalk",con=con)
dfw_cw_blks = cw.merge(dfw_tdata,left_on='trct',right_on=['GEOID'])[['tabblk2020','GEOID']]
cw_mapdict = dict(zip(dfw_cw_blks['tabblk2020'],dfw_cw_blks['GEOID']))
#get od data for the dfw_cw_blks in 2022 and 2012

dfs = []
for y in ['2012','2017','2022']:
    q = generate_query(data_type='od',perspective='home',
                       job_type='primary',year=y,geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
    od_df = pull_data(query=q,crsr=cur,rename=True)
    retype(od_df)
    dfs.append(od_df)

od_primary = pd.concat(dfs)
od_primary['dest_tract'] = od_primary['w_geocode'].map(cw_mapdict)
od_primary['orig_tract'] = od_primary['h_geocode'].map(cw_mapdict)

od_set = od_primary.groupby(['year','orig_tract','dest_tract']).agg({'total':'sum',
                                                           'goods':'sum',
                                                           'transp':'sum',
                                                           'other':'sum',}).reset_index()

#get the proportion of each pair to the total for the tract
ods_set_list = []
for each_year in od_set['year'].unique():
    ods = od_set[od_set['year'] == each_year]
    tm = ods.groupby('orig_tract').agg({'total':'sum'}).reset_index()
    tm_dict = dict(zip(tm['orig_tract'],tm['total']))
    ods['total_orig'] = ods['orig_tract'].map(tm_dict)
    ods['proportion_of_orig_to_dest'] = ods.apply(lambda x: x['total']/tm_dict[x['orig_tract']],axis=1)
    ods_set_list.append(ods)

od_set = pd.concat(ods_set_list)

#adjusted destiantion for acs
adj_dict = dict(zip(time_of_day['GEOID'],time_of_day['B08302_001E']))

#time_blocks 
time_blocks = {'t12am_t5am':['B08302_002E'],
               't5am_t8am':['B08302_003E','B08302_004E','B08302_005E',
                            'B08302_006E','B08302_007E','B08302_008E'],
               't9am_t12pm':['B08302_009E','B08302_010E','B08302_011E',
                            'B08302_012E','B08302_013E'],
               't12pm_t4pm':['B08302_014E'],
               't4pm_t12am':['B08302_015E']}



merge_ts = []
for k in time_blocks.keys():
    time_of_day[k] = time_of_day[time_blocks[k]].sum(axis=1)
    time_of_day[f"proportion_{k}"] = time_of_day[k]/time_of_day['B08302_001E']
    merge_ts.append(f"proportion_{k}")

od_set['b08302_tot'] = od_set['orig_tract'].map(adj_dict)
od_set['adjusted_dest'] = od_set['b08302_tot'] * od_set['proportion_of_orig_to_dest']

for k in merge_ts:
    adj_dict = dict(zip(time_of_day['GEOID'],time_of_day[k]))
    od_set[f'{k}'] = ods['orig_tract'].map(adj_dict)
    od_set[f"{k.replace("proportion_","total_")}"] = od_set[k] * od_set['adjusted_dest']

ods_matrix = od_set[['orig_tract','dest_tract','year','total','adjusted_dest',
     'total_t12am_t5am','total_t5am_t8am','total_t9am_t12pm','total_t12pm_t4pm','total_t4pm_t12am']]

ods_matrix.query("year == '2022'")



#if employment growth followed the same spatial patterns of 2012 then the 2022 adjusted_dest should equal the 2012 adjusted_dest
ods_matrix['kv_tup'] = list(zip(ods_matrix['orig_tract'],ods_matrix['dest_tract']))
y_2012 = ods_matrix.query("year == 2012")
y_2022 = ods_matrix.query("year == 2022")

exp_2012 = dict(zip(y_2012['kv_tup'],y_2012['adjusted_dest']))
y_2022 = y_2022.fillna(0)
y_2022['exp_2012'] = y_2022['kv_tup'].map(exp_2012)
y_2022['actual_2022'] = y_2022['adjusted_dest'] 
y_2022['diff_exp_act'] = y_2022['actual_2022'] - y_2022['exp_2012']

y_2022.query("orig_tract == '48113020500'").sort_values(by='diff_exp_act',ascending=True)

#get tract geometry
t20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geoid from tracts_2020_geom",con=con)
t20 = gpd.GeoDataFrame(t20,geometry=gpd.GeoSeries.from_wkt(t20['geom_wkt']),crs='epsg:4326')
t20.drop(columns='geom_wkt',inplace=True)

sdata = t20.merge(y_2022,left_on='GEOID',right_on='dest_tract')
sdata = sdata.to_crs("EPSG:3857")

import mapclassify
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import contextily as ctx

def quickplot(orig_geoid='48113013202'):
    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    print("Plotting...")
    sdata.query('orig_tract == @orig_geoid').plot(
        ax=ax,
        column="diff_exp_act",
        cmap="magma_r",
        alpha=0.8,
        linewidth=2,
        zorder=4,
        legend=True,
        legend_kwds={"shrink": 0.3} ,)
    
    sdata.query('dest_tract == @orig_geoid').iloc[0:1].plot(
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

#get uni pairs
a = ods_matrix['orig_tract'].tolist()
b = ods_matrix['dest_tract'].tolist()
c = pd.concat([pd.DataFrame([a,b]).T,pd.DataFrame([b,a]).T])
uni_pairs = c.drop_duplicates()


## 3 get a buffer around a city and then query it
#get a random example from houston
import geopandas as gpd
example = gpd.read_file(gpd.datasets.get_path('naturalearth_cities')).query("name.str.contains('Houston')")
example = example.to_crs("EPSG:2278")
example.geometry = example.buffer(5 * 5280)

#get the wkt
wkts = transform_to_wkt(gdf=example)

#get the blocks as a list that intersect it 
isect = id_intersections(wkt=wkts.iloc[0], spath=spath)['geocode'].tolist()

#run the query
qb = generate_query(data_type='wac',job_type='all',year='2020',geocodes = isect)
dfb = pull_data(qb,spath=spath,rename=True)



## 4 quickly generate a plot of the destination for workers from within a 2 mile buffer of houston
#get a random example from houston
import geopandas as gpd
example = gpd.read_file(gpd.datasets.get_path('naturalearth_cities'))
htown = example.query("name.str.contains('Houston')")
htown = htown.to_crs("EPSG:2278")
#five miles from downtown houston
htown.geometry = htown.buffer(5 * 5280)

#get the wkt
wkts = transform_to_wkt(gdf=htown)

#get the blocks as a list that intersect it 
isect_gdf = id_intersections(wkt=wkts.iloc[0], spath=spath,return_geom=True)
isect = isect_gdf['geocode'].tolist()

#run the query
qb = generate_query(data_type='od',perspective='work',job_type='primary',year='2021',geocodes = isect)
dfb = pull_data(qb,spath=spath,rename=True)
totals = dfb.groupby('h_geocode').agg({'total':'sum'}).reset_index()


#pull geometries within 50 miles 
htown50 = example.query("name.str.contains('Houston')")
htown50 = htown50.to_crs("EPSG:2278")
#five miles from downtown houston
htown50.geometry = htown50.buffer(5 * 5280)
wkts50 = transform_to_wkt(gdf=htown50)
b50mi = id_intersections(wkt=wkts50.iloc[0], spath=spath,return_geom=True)

#add in the data
m = b50mi.merge(totals,left_on='geocode',right_on='h_geocode',how='left').fillna(0)

m.plot(column='total',legend=True,figsize=(10,10))

#could easily dissolve into tracts
m['tract_id'] = m['geocode'].str[:11]
m2 = m[['tract_id','geometry','total']].dissolve(by='tract_id',aggfunc='sum')
m2.plot(column='total',legend=True,figsize=(10,10))