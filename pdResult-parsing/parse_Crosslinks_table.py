import sqlite3
import csv
from pathlib import Path
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
    Exports the CrossLinks table from a selected PDResult SQLite database 
    to a CSV file. 
    
    The input PDResult file is selected from the input directory. 
    The exported CSV file is stored in the output directory 
    and named after the selected input file. 
    
    :return: None. 
    """
    selected_file = choose_file(input_dir, ".pdResult")

    conn = sqlite3.connect(selected_file)
    cursor = conn.cursor()

    result = cursor.execute('SELECT * FROM "CrossLinks"')

    output_file = output_dir / f"{selected_file.stem}_Crosslinks.csv"
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        writer.writerow([column[0] for column in result.description])

        writer.writerows(result)

    conn.close()

    print(f"\nOutput file was written to: {output_file}")


if __name__ == "__main__":
    main()