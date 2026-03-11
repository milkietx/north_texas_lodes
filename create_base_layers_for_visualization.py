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
all_y = []
for y in ['2008','2013','2018','2023']:
    q = generate_query(data_type='rac',
                        job_type='primary',
                        year=y,
                        geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
    rac_df = pull_data(query=q,crsr=cur,rename=True)

    q = generate_query(data_type='wac',
                        job_type='primary',
                        year=y,
                        geocodes=dfw_cw_blks['tabblk2020'].unique().tolist())
    wac_df = pull_data(query=q,crsr=cur,rename=True)

    lodes_df = pd.merge(wac_df,
                        rac_df,
                        left_on='w_geocode',
                        right_on='h_geocode',
                        suffixes=('_wac','_rac'),
                        how='outer').fillna(0)
    lodes_df['geocode'] = lodes_df.apply(lambda x: x['w_geocode'] if x['h_geocode'] == 0 else x['h_geocode'],axis=1)

    lodes_df['year'] = y
    all_y.append(lodes_df)

all_y_single = pd.concat(all_y)
all_y_piv = all_y_single.pivot(index=['geocode'],
             columns=['year'])
            
rev_columns = [f"{x[0]}_{x[1]}" for x in all_y_piv.columns]
all_y_piv.columns = rev_columns

all_y_piv = all_y_piv.drop(columns=['w_geocode_2008',
 'w_geocode_2013',
 'w_geocode_2018',
 'w_geocode_2023',
 'h_geocode_2008',
 'h_geocode_2013',
 'h_geocode_2018',
 'h_geocode_2023',
 'year_wac_2008',
 'year_wac_2013',
 'year_wac_2018',
 'year_wac_2023',
  'year_rac_2008',
 'year_rac_2013',
 'year_rac_2018',
 'year_rac_2023'])

b20 = pd.read_sql(r"select ST_AsText(geom) as geom_wkt, geocode from blocks_2020_geom",con=con)
b20 = gpd.GeoDataFrame(b20,geometry=gpd.GeoSeries.from_wkt(b20['geom_wkt']),crs='epsg:4326')
b20['county'] = b20['geocode'].str[2:5]
dfw20 = b20.query("county in @dfw_counties")

lodes_dfw = dfw20.merge(all_y_piv,left_on='geocode',right_on='geocode',how='left').fillna(0)

#get density of jobs 
#get density of employment
lodes_dfw = lodes_dfw.to_crs("EPSG:2276")
lodes_dfw['area_sqmi'] = lodes_dfw.area/(5280*5280)
lodes_dfw['wf_per_sqmi_2023'] = lodes_dfw['tot_rac_2023']/lodes_dfw['area_sqmi']
lodes_dfw['jobs_per_sqmi_2023'] = lodes_dfw['tot_wac_2023']/lodes_dfw['area_sqmi']
lodes_dfw['wf_per_sqmi_2018'] = lodes_dfw['tot_rac_2018']/lodes_dfw['area_sqmi']
lodes_dfw['jobs_per_sqmi_2018'] = lodes_dfw['tot_wac_2018']/lodes_dfw['area_sqmi']
lodes_dfw['wf_per_sqmi_2013'] = lodes_dfw['tot_rac_2013']/lodes_dfw['area_sqmi']
lodes_dfw['jobs_per_sqmi_2013'] = lodes_dfw['tot_wac_2013']/lodes_dfw['area_sqmi']
lodes_dfw['wf_per_sqmi_2008'] = lodes_dfw['tot_rac_2008']/lodes_dfw['area_sqmi']
lodes_dfw['jobs_per_sqmi_2008'] = lodes_dfw['tot_wac_2008']/lodes_dfw['area_sqmi']

lodes_dfw = lodes_dfw.to_crs("EPSG:3857")

lodes_dfw.to_file(r"C:\Users\cmg0530\Projects\lodes_north_texas\north_texas_lodes\Data Storage\Zones.gpkg",
               layer='Background_LODES_v1')


lodes_dfw['wf_per_sqmi_2023'].quantile(.4)