import sqlite3
import csv
from pathlib import Path

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
    cursor = conn.cursor()

    result = cursor.execute('SELECT * FROM "CrossLinks"')

    output_file = output_dir / f"{selected_file.stem}_Crosslinks.csv"
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        writer.writerow([column[0] for column in result.description])

        writer.writerows(result)

    conn.close()

    print(f"\nOutput files were written to: {output_file}")


if __name__ == "__main__":
    main()