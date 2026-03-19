import pandas as pd
import json

cull = pd.read_excel('data/260201_COTS-Cull-Data-Ewels.xlsx')
edna = pd.read_excel('data/eDNA data_ALL_20260225.xlsx', sheet_name='eDNA_data_ALL')

def find_spatial(cols):
    return [c for c in cols if 'lat' in str(c).lower() or 'lon' in str(c).lower() or 'site' in str(c).lower() or 'reef' in str(c).lower()]

info = {
    'cull_spatial_cols': find_spatial(cull.columns),
    'edna_spatial_cols': find_spatial(edna.columns),
    'cull_all_cols': cull.columns.tolist()[:30],
    'edna_all_cols': edna.columns.tolist()[:30]
}

with open('cols.txt', 'w') as f:
    json.dump(info, f, indent=2)
