import sys
sys.path.append(r'C:\Users\cmg0530\Code Library\lodes-dataset-construction\main')
from analysis import *
sys.path.append(r'C:\Users\cmg0530\Code Library\gis_suite')
from census_pulls import *
import pandas as pd
import os
import geopandas as gpd

#connect to od
spath = r"C:\Users\cmg0530\Projects\lodes_package\lodes_tx_slim.db"
con,cur = connect_to_od(spath=spath)

def prep_employment_data():
    #get total residential employment for 2023
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

    #get city crosswalk
    cw = pd.read_sql("select tabblk2020,trct,stplcname,cty from tx_xwalk",con=con)
    cw['county_id'] = cw['cty'].str[2:5]
    dfw_cw_blks = cw[cw['county_id'].isin(dfw_counties)]
    city_mapdict = dict(zip(dfw_cw_blks['tabblk2020'],dfw_cw_blks['stplcname']))


    #get all the od data for 2013, 2018, and 2023
    dfs = []
    for y in ['2013','2018','2023']:
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
    return od_primary,cw_mapdict,city_mapdict,dfw_cw_blks,dfw_tdata

def assign_tracts(od_primary,cw_mapdict,city_mapdict):
    #assign tracts and origin city
    od_primary['dest_tract'] = od_primary['w_geocode'].map(cw_mapdict)
    od_primary['orig_tract'] = od_primary['h_geocode'].map(cw_mapdict)
    od_primary['orig_city'] = od_primary['h_geocode'].map(city_mapdict)
    return od_primary


#read in a study area and assign to od_primary
def assign_study_area(od_primary,layer_name="UTSW_v1"):
    study_area = gpd.read_file(r"C:\Users\cmg0530\Data Storage\GIS Data\LODES_StudyAreas.gdb", 
                                layer=layer_name)
    study_area_dict = dict(zip(study_area['geocode'],[layer_name] * len(study_area['geocode'])))
    od_primary['study_area'] = od_primary['w_geocode'].map(study_area_dict)
    od_primary['study_area'] = od_primary['study_area'].fillna('not study area')
    return od_primary


#aggregate to tract level
def aggregate_to_origin(od_primary,geography='tract'):
    od_set = od_primary.groupby(['year',
                                f'orig_{geography}',
                                'study_area']).agg({'total':'sum',
                                                    'goods':'sum',
                                                    'transp':'sum',
                                                    'other':'sum',}).reset_index()
    return od_set


#get the percentage of total commuters for each OD pair year
def adjust_od_set_by_year(od_set,geography='tract'):
    ods_set_list = []
    for each_year in od_set['year'].unique():
        ods = od_set[od_set['year'] == each_year]
        tm = ods.groupby(f'orig_{geography}').agg({'total':'sum'}).reset_index()
        tm_dict = dict(zip(tm[f'orig_{geography}'],tm['total']))
        ods['total_origin_lodes'] = ods[f'orig_{geography}'].map(tm_dict)
        ods['proportion_of_origin_to_dest'] = ods.apply(lambda x: x['total']/tm_dict[x[f'orig_{geography}']],axis=1).fillna(0)
        ods_set_list.append(ods)

    #pair with each for value
    od_set = pd.concat(ods_set_list)

    #pivot by year
    od_pivot = od_set.pivot(index=[f'orig_{geography}','study_area'],
                columns=['year'],
                values=['total','total_origin_lodes','proportion_of_origin_to_dest'])

    rev_columns = [f"{x[0]}_{x[1]}" for x in od_pivot.columns]
    od_pivot.columns = rev_columns
    od_pivot = od_pivot.reset_index().fillna(0)

    return od_pivot


def adjust_to_acs(od_pivot,dfw_tdata):
    #adjusted destiantion for acs
    adj_dict = dict(zip(dfw_tdata['GEOID'],dfw_tdata['C24020_001E']))
    od_pivot['origin_total_acs_2023'] = od_pivot['orig_tract'].map(adj_dict)

    #expected by different year 
    od_pivot['expected_from_2013_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2013'].fillna(0) * od_pivot['origin_total_acs_2023'] 
    od_pivot['expected_from_2018_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2018'].fillna(0) * od_pivot['origin_total_acs_2023'] 
    od_pivot['expected_from_2023_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2023'].fillna(0) * od_pivot['origin_total_acs_2023'] 

    return od_pivot

def merge_in_geometry(od_pivot):
    #bring in geometry for origin
    t20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geoid from tracts_2020_geom",con=con)
    t20 = gpd.GeoDataFrame(t20,geometry=gpd.GeoSeries.from_wkt(t20['geom_wkt']),crs='epsg:4326')
    t20.drop(columns='geom_wkt',inplace=True)
    sdata = t20.merge(od_pivot,left_on='GEOID',right_on='orig_tract',how='left').fillna(0)
    sdata = sdata.to_crs("EPSG:3857")
    return sdata

#read in employment data and prepare dictioanries and OD data
od_primary,cw_mapdict,city_mapdict,dfw_cw_blks,dfw_tdata = prep_employment_data()

#assign tracts and cities
od_primary = assign_tracts(od_primary,cw_mapdict,city_mapdict)

#pull in and assign study areas
od_primary_utsw = assign_study_area(od_primary,layer_name="UTSW_v1")
od_primary_medcity = assign_study_area(od_primary,layer_name="MedicalCity_v1")
od_primary_methodistftw= assign_study_area(od_primary,layer_name="MethodistFTW_v1")
od_primary_presbydenton = assign_study_area(od_primary,layer_name="PresbyterianDenton_v1")

#process data into dataset
od_frames = {}
for stuare in [od_primary_utsw,
               od_primary_medcity,
               od_primary_methodistftw,
               od_primary_presbydenton]:
    od_set = aggregate_to_origin(stuare, geography='tract')
    od_pivot = adjust_od_set_by_year(od_set, geography='tract')
    od_pivot = adjust_to_acs(od_pivot,dfw_tdata)
    sdata = merge_in_geometry(od_pivot)
    x = [y for y in list(od_pivot['study_area'].unique()) if y != 'not study area'][0]
    od_frames[f"{x}_tract"] = sdata
    
    od_setc = aggregate_to_origin(stuare, geography='city')
    od_pivotc = adjust_od_set_by_year(od_setc, geography='city')
    od_frames[f"{x}_city"] = od_pivotc

#save out
for x in od_frames.keys():
    if 'tract' in x:
        od_frames[x].to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
                    layer=x)
    if 'city' in x:
        od_frames[x].to_excel(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\CityLevelData.xlsx",
                    sheet_name=x)



#%%random extraneous mapping

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


