import sqlite3
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent

INPUT_FILE = SCRIPT_DIR / "input" / "XlinkX_CDI_mainlin_rep1.pdResult"
OUTPUT_DIR = SCRIPT_DIR / "output" / "pdResult_export"

OUTPUT_DIR.mkdir(exist_ok=True)


conn = sqlite3.connect(INPUT_FILE)

# Alle Tabellen finden
tables = pd.read_sql_query(
    """
    SELECT name
    FROM sqlite_master
    WHERE type='table'
    ORDER BY name
    """,
    conn,
)["name"]

print(f"Found {len(tables)} tables.\n")

# Excel-Datei mit allen Tabellen
excel_path = OUTPUT_DIR / "all_tables.xlsx"

with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:

    for table in tables:
        print(f"Exporting {table}...")

        df = pd.read_sql_query(f'SELECT * FROM "{table}"', conn)

        # CSV
        df.to_csv(
            OUTPUT_DIR / f"{table}.csv",
            index=False,
        )

        # Excel Sheet
        sheet = table[:31]  # Excel erlaubt max. 31 Zeichen
        df.to_excel(writer, sheet_name=sheet, index=False)

conn.close()

print("\nDone!")
print(f"CSV files written to: {OUTPUT_DIR}")
print(f"Excel workbook written to: {excel_path}")