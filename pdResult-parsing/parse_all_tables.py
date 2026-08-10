import sqlite3
from pathlib import Path
import pandas as pd

script_dir = Path(__file__).parent
input_dir = script_dir / "input"
output_dir = script_dir / "output"

output_dir.mkdir(exist_ok=True)


def choose_file() -> Path:
    files = sorted(
        list(input_dir.glob("*.pdResult"))
    )

    if not files:
        raise FileNotFoundError(f"No input file available in {input_dir}")

    print("Available input:")
    for i, file in enumerate(files, start=1):
        print(f"{i}: {file.name}")
    print()

    while True:
        choice = input(f"Select file [1-{len(files)}]: ")

        try:
            choice = int(choice)
        except ValueError:
            print("Please enter a number.")
            continue

        if 1 <= choice <= len(files):
            return files[choice - 1]

        print("Invalid selection.")


def main():
    selected_file = choose_file()

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
    print(f"CSV files written to: {specific_output_dir}")
    print(f"Excel workbook written to: {excel_path}")


if __name__ == "__main__":
    main()