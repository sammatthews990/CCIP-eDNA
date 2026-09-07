import pandas as pd

cull_file = "data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
edna_file = "data/eDNA data_ALL_20260528.xlsx"

df_cull = pd.read_excel(cull_file, sheet_name="Cull")
df_manta = pd.read_excel(cull_file, sheet_name="Manta")
df_rhis = pd.read_excel(cull_file, sheet_name="RHIS")

df_edna = pd.read_excel(edna_file, sheet_name="eDNA_data_ALL")

print(f"Cull shape: {df_cull.shape}")
print("Cull date col:", "SurveyDate" in df_cull.columns)
print("Cull min date:", df_cull['SurveyDate'].min(), "max date:", df_cull['SurveyDate'].max())

print(f"\nManta shape: {df_manta.shape}")
print("Manta date col:", "SurveyTime" in df_manta.columns)
print("Manta min date:", df_manta['SurveyTime'].min(), "max date:", df_manta['SurveyTime'].max())

print(f"\nRHIS shape: {df_rhis.shape}")
print("RHIS date col:", "SurveyTime" in df_rhis.columns)
print("RHIS min date:", df_rhis['SurveyTime'].min(), "max date:", df_rhis['SurveyTime'].max())

print(f"\neDNA shape: {df_edna.shape}")
print("eDNA min date:", df_edna['Date'].min(), "max date:", df_edna['Date'].max())
print("eDNA unique reefs:", df_edna['ReefName'].nunique())
print("eDNA unique site locations (Lat/Long):", df_edna[['Lat', 'Long']].drop_duplicates().shape[0])
print("eDNA unique years per reef count:")
reef_years = df_edna.groupby('ReefName')['Year'].nunique()
print(reef_years.value_counts())

print("\nSample of eDNA Reefs with >1 year:")
print(reef_years[reef_years > 1].head(10))
