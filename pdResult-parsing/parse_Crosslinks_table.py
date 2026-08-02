import sqlite3
import csv
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

INPUT_FILE = SCRIPT_DIR / "input" / "Annika_CDI_mainlin_rep2.pdResult"
OUTPUT_FILE = SCRIPT_DIR / "output" / "Annika_CDI_mainlin_rep2_CrossLinks.csv"

conn = sqlite3.connect(INPUT_FILE)
cursor = conn.cursor()

# CrossLinks-Tabelle auslesen
result = cursor.execute('SELECT * FROM "CrossLinks"')

with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)

    # Spaltennamen schreiben
    writer.writerow([column[0] for column in result.description])

    # Alle Zeilen schreiben
    writer.writerows(result)

conn.close()

print(f"Export erfolgreich: {OUTPUT_FILE}")