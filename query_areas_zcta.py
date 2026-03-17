import sys
sys.path.append(r'C:\Users\cmg0530\Code Library\lodes-dataset-construction\main')
from analysis import *
sys.path.append(r'C:\Users\cmg0530\Code Library\gis_suite')
from census_pulls import *
import pandas as pd
import os
import geopandas as gpd
import numpy as np

#connect to od
spath = r"C:\Users\cmg0530\Projects\lodes_package\lodes_tx_slim.db"
con,cur = connect_to_od(spath=spath)

def prep_employment_data():
    #get total residential employment for 2023
    total_employment = census_pull(tables=["C24020"],geom="zcta",year="2023")

    #dfw counties
    dfw_counties = ['085','113','121','139','143',
                    '221','231','251','257','349',
                    '363','367','397','425','439','497']

    #get crosswalked blocks
    cw = pd.read_sql("select tabblk2020, zcta,zctaname,cty from tx_xwalk",con=con)
    dfw_zips = cw[cw["cty"].str[-3:].isin(dfw_counties)]
    
    relev_zips = dfw_zips.merge(total_employment,left_on='zcta',right_on=['GEOID'])['zcta']
    dfw_zdata = total_employment.query("GEOID in @relev_zips")

    cw_mapdict = dict(zip(dfw_zips['tabblk2020'],dfw_zips['zcta']))

    #get city crosswalk
    cw = pd.read_sql("select tabblk2020,stplcname,cty from tx_xwalk",con=con)
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
    return od_primary,cw_mapdict,city_mapdict,dfw_cw_blks,dfw_zdata

def assign_zcta(od_primary,cw_mapdict,city_mapdict):
    #assign tracts and origin city
    od_primary['dest_zcta'] = od_primary['w_geocode'].map(cw_mapdict)
    od_primary['orig_zcta'] = od_primary['h_geocode'].map(cw_mapdict)
    od_primary['orig_city'] = od_primary['h_geocode'].map(city_mapdict)
    return od_primary


#read in a study area and assign to od_primary
def assign_study_area(od_primary,layer_name="UTSW_v1"):
    study_area = gpd.read_file(r"C:\Users\cmg0530\Data Storage\GIS Data\LODES_StudyAreas.gdb", 
                                layer=layer_name)
    study_area_dict = dict(zip(study_area['geocode'],[layer_name] * len(study_area['geocode'])))
    od_primary_2 = od_primary.copy()
    od_primary_2['study_area'] = od_primary_2['w_geocode'].map(study_area_dict)
    od_primary_2['study_area'] = od_primary_2['study_area'].fillna('not study area')
    return od_primary_2


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


def adjust_to_acs(od_pivot,dfw_tdata,geography='zcta'):
    #adjusted destiantion for acs
    adj_dict = dict(zip(dfw_tdata['GEOID'],dfw_tdata['C24020_001E']))
    od_pivot['origin_total_acs_2023'] = od_pivot[f'orig_{geography}'].map(adj_dict)

    #expected by different year 
    od_pivot['expected_from_2013_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2013'].fillna(0) * od_pivot['origin_total_acs_2023'] 
    od_pivot['expected_from_2018_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2018'].fillna(0) * od_pivot['origin_total_acs_2023'] 
    od_pivot['expected_from_2023_proportions_val_acs_2023'] = od_pivot['proportion_of_origin_to_dest_2023'].fillna(0) * od_pivot['origin_total_acs_2023'] 

    return od_pivot

def merge_in_geometry(od_pivot,geography='zcta'):
    #bring in geometry for origin
    t20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, GEOID20 as GEOID from zcta_2020_geom",con=con)
    t20 = gpd.GeoDataFrame(t20,geometry=gpd.GeoSeries.from_wkt(t20['geom_wkt']),crs='epsg:4326')
    t20.drop(columns='geom_wkt',inplace=True)
    sdata = t20.merge(od_pivot,left_on='GEOID',right_on=f'orig_{geography}',how='left').fillna(0)
    sdata = sdata.to_crs("EPSG:3857")
    return sdata

#read in employment data and prepare dictioanries and OD data
od_primary,cw_mapdict,city_mapdict,dfw_cw_blks,dfw_zdata = prep_employment_data()

#assign tracts and cities
od_primary = assign_zcta(od_primary,cw_mapdict,city_mapdict)

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
    od_set = aggregate_to_origin(od_primary=stuare, geography='zcta')
    od_pivot = adjust_od_set_by_year(od_set, geography='zcta')
    od_pivot = adjust_to_acs(od_pivot,dfw_zdata)
    sdata = merge_in_geometry(od_pivot)
    x = [y for y in list(sdata['study_area'].unique()) if (y != 'not study area') and len(str(y)) > 2][0]
    od_frames[f"{x}_zcta"] = sdata
    
    od_setc = aggregate_to_origin(stuare, geography='city')
    od_pivotc = adjust_od_set_by_year(od_setc, geography='city')
    od_frames[f"{x}_city"] = od_pivotc

#save out
for x in od_frames.keys():
    if 'zcta' in x:
        od_frames[x].to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
                    layer=f"{x}")
    if 'city' in x:
        od_frames[x].to_excel(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\CityLevelData.xlsx",
                    sheet_name=x)

#read and process lc data
lc = pd.read_csv(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\lightcast_zcta_health.csv")
lc['ZCTA'] = lc['Area'].astype(str)

def format_for_ratios(lc):
    pivotye = lc.pivot(index=['ZCTA'],columns=['Occupation',"Year"],values='Resident Workers')
    pivotye.columns = [f"{x[0]}_{x[1]}" for x in pivotye.columns]
    pivotye = pivotye.reset_index()
    return pivotye

lc2 = format_for_ratios(lc)

#read geoms
z20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, GEOID20 as GEOID from zcta_2020_geom",con=con)
z20 = gpd.GeoDataFrame(z20,geometry=gpd.GeoSeries.from_wkt(z20['geom_wkt']),crs='epsg:4326')
z20.drop(columns='geom_wkt',inplace=True)
z20dfw = z20[z20['GEOID'].isin(dfw_zdata['GEOID'])].copy()

lcsdata = z20dfw.merge(lc2,left_on='GEOID',right_on=f'ZCTA',how='inner').fillna(0)
lcsdata = lcsdata.to_crs("EPSG:3857")
for x in lcsdata.columns:
    if x not in ['GEOID','geometry']:
        lcsdata[x]= lcsdata[x].astype(float,errors='ignore').fillna(0)
lcsdata.to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
                    layer=f"zcta_lightcast_data_columnar")



def calculate_tiers(tier_1_occs,tier_2_occs,tier_3_occs,tier_name,year,lcsdata,tier_4_occs=None):
    lcsdata[f'{tier_name}_tier1_2025'] = lcsdata[[f"{y}_{year}" for y in tier_1_occs]].sum(axis=1)
    lcsdata[f'{tier_name}_tier2_2025'] = lcsdata[[f"{y}_{year}" for y in tier_2_occs]].sum(axis=1)
    lcsdata[f'{tier_name}_tier3_2025'] = lcsdata[[f"{y}_{year}" for y in tier_3_occs]].sum(axis=1)
    if tier_4_occs != None:
        lcsdata[f'{tier_name}_tier4_2025'] = lcsdata[[f"{y}_{year}" for y in tier_4_occs]].sum(axis=1)
    return lcsdata


#nursing occupation tiers
#tier 1: CNA  31-1131
#tier 2: LVN 29-2061
#tier 3: RN 29-1141
#tier 4: Nurse Practicioner 29-1171
tier_1_occs = ['31-1131']
tier_2_occs = ['29-2061']
tier_3_occs = ['29-1141']
tier_4_occs = ['29-1171']
tier_name = 'nursing'
year = '2025'
lcsdata_nursing = calculate_tiers(tier_1_occs=tier_1_occs,
                                  tier_2_occs=tier_2_occs,
                                  tier_3_occs=tier_3_occs,
                                  tier_4_occs=tier_4_occs,
                                  tier_name=tier_name,
                                  year=year,
                                  lcsdata=lcsdata)



#tier 1: medical assistant, orderlies 31-9092, 31-1132
#tier 2: cardiovascular tech + surgical technologist 29-2031, 29-2055
#tier 3: phyisician assistant 29-1071
tier_1_occs = ['31-9092','31-1132']
tier_2_occs = ['29-2031','29-2055']
tier_3_occs = ['29-1071']
tier_4_occs = None
tier_name = 'surgery'
year = '2025'
lcsdata_surgery = calculate_tiers(tier_1_occs=tier_1_occs,
                                  tier_2_occs=tier_2_occs,
                                  tier_3_occs=tier_3_occs,
                                  tier_4_occs=tier_4_occs,
                                  tier_name=tier_name,
                                  year=year,
                                  lcsdata=lcsdata)


def tier_ratio(lcsdata, tier_name,year):

    cols = [y for y in lcsdata.columns if tier_name in y]
    tiers = []
    for q in cols:
        tiers.append(q.split("_")[-2])
    #what are the ratios
    
    #ratio 1 
    lcsdata[f"{tier_name}_ratio_t2_t1_{year}"] = (lcsdata[f"{tier_name}_tier1_{year}"]/lcsdata[f"{tier_name}_tier2_{year}"]).fillna(0).replace(np.inf,0)
    #tier 2 to tier 1

    #ratio 2
    lcsdata[f"{tier_name}_ratio_t3_t1t2_{year}"] = ((lcsdata[f"{tier_name}_tier1_{year}"]+lcsdata[f"{tier_name}_tier2_{year}"])/lcsdata[f"{tier_name}_tier3_{year}"]).fillna(0).replace(np.inf,0)
    #tier 3 to tier 1 + 2
    
    #optional
    if "tier4" in tiers:
        #tier 4 to tier 1 + 2 + 3
        lcsdata[f"{tier_name}_ratio_t4_t1t2t3_{year}"] = ((lcsdata[f"{tier_name}_tier1_{year}"]+lcsdata[f"{tier_name}_tier2_{year}"]+lcsdata[f"{tier_name}_tier3_{year}"])/lcsdata[f"{tier_name}_tier4_{year}"]).fillna(0).replace(np.inf,0)
    
    return lcsdata

lcs_nursing = tier_ratio(lcsdata=lcsdata_nursing,tier_name='nursing',year='2025')
lcs_surgery = tier_ratio(lcsdata=lcsdata_surgery,tier_name='surgery',year='2025')

lcs_surgery.to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
                    layer=f"Employment_Ratios_Surgery_v1")
lcs_nursing.to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
                    layer=f"Employment_Ratios_Nursing_v1")

#%%random extraneous mapping

import mapclassify
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import contextily as ctx

def quickplot():
    fig, ax = plt.subplots(1, 1, figsize=(19, 19))
    print("Plotting...")
    lcsdata.query("difference_2012_2022_val2024 > -25").plot(
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


