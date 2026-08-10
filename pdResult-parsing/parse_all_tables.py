import sqlite3
from pathlib import Path
import pandas as pd
import sys

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))
from helpers.input_file_selection import choose_file

script_dir = Path(__file__).parent
input_dir = script_dir / "input"
output_dir = script_dir / "output"

output_dir.mkdir(exist_ok=True)


def main():
    """ 
    Exports all tables from a selected PDResult SQLite database 
    to individual CSV files and a combined Excel workbook.

    The input PDResult file is selected from the input directory. 
    The exported files are stored in a dedicated output directory 
    named after the selected input file. 

    :return: None. 
    """
    selected_file = choose_file(input_dir, ".pdResult")

    conn = sqlite3.connect(selected_file)

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

    specific_output_dir = output_dir / f"{selected_file.stem}_export"
    specific_output_dir.mkdir(exist_ok=True)

    excel_path = specific_output_dir / "all_tables.xlsx"

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:

        for table in tables:
            print(f"Exporting {table}...")

            df = pd.read_sql_query(f'SELECT * FROM "{table}"', conn)

            df.to_csv(
                specific_output_dir / f"{table}.csv",
                index=False,
            )

            sheet = table[:31]  
            df.to_excel(writer, sheet_name=sheet, index=False)

    conn.close()

    print("\nDone!")
    print(f"CSV files were written to: {specific_output_dir}")
    print(f"Excel workbook was written to: {excel_path}")


if __name__ == "__main__":
    main()