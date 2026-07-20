import sqlite3
import csv
from pathlib import Path

INPUT_FILE = "Cas9Cube_r1.pdResult"
OUTPUT_FILE = "Cas9Cube_r1_CrossLinks.csv"

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