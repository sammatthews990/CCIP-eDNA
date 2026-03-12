import pandas as pd

def inspect_excel(path, name):
    try:
        xl = pd.ExcelFile(path)
        res = f"=== {name} ===\nSheets: {xl.sheet_names}\n"
        for sheet in xl.sheet_names:
            cols = xl.parse(sheet, nrows=0).columns.tolist()
            res += f"  Sheet '{sheet}' cols: {cols}\n"
        return res
    except Exception as e:
        return f"Error with {name}: {e}\n"

with open("sheets_info.txt", "w") as f:
    f.write(inspect_excel("data/260201_COTS-Cull-Data-Ewels.xlsx", "Cull Data"))
    f.write("\n")
    f.write(inspect_excel("data/eDNA data_ALL_20260225.xlsx", "eDNA Data"))
