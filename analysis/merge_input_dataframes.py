from pathlib import Path
import pandas as pd
import sys

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))
from helpers.input_file_selection import choose_file


project_dir = Path(__file__).parent
input_dir = project_dir / "input_results"
output_dir = project_dir / "merged_dfs"

output_dir.mkdir(exist_ok=True)

def merge_selected_directories(selected_dirs):
    ca_dfs = []
    exact_dfs = []

    for dir in selected_dirs:
        ca_files = list(dir.glob("*_CA.csv"))
        exact_files = list(dir.glob("*_exact.csv"))

        ca_df = pd.read_csv(ca_files[0])
        exact_df = pd.read_csv(exact_files[0])

        if list(ca_df.columns) != list(exact_df.columns): 
            raise ValueError(
                "Columns do not match."
            )
        
        ca_dfs.append(ca_df)
        exact_dfs.append(exact_df)

    merged_ca = pd.concat(ca_dfs, ignore_index=True)
    merged_exact = pd.concat(exact_dfs, ignore_index=True)

    return merged_ca, merged_exact


def main():
    selected_dirs = []

    while True: 
        folders = sorted(f for f in input_dir.iterdir() if f.is_dir())
        if not folders:
            raise FileNotFoundError(f"No input folders found in {input_dir}")

        print("Available input folders:")
        for i, folder in enumerate(folders, start=1):
            print(f"{i}: {folder.name}")
        print()

        while True:
            choice = input(f"Select folder [1-{len(folders)}]: ")
            if choice == ".":
                selected_folder = None
                break
            try:
                choice = int(choice)
                if 1 <= choice <= len(folders):
                    selected_folder = folders[choice - 1]
                    break
            except ValueError:
                pass
            print("Invalid selection.")

        if selected_folder is None: 
            break 
        selected_dirs.append(selected_folder)

    if not selected_dirs: 
        print("No directories selected.")
        return 
    
    try: 
        merged_ca, merged_exact = merge_selected_directories(selected_dirs)

    except ValueError as error: 
        print(f"\nError: {error}")
        return
    
    ca_output = output_dir / "merged_ca.csv"
    exact_output = output_dir / "merged_exact.csv"

    merged_ca.to_csv(ca_output, index=False)
    merged_exact.to_csv(exact_output, index=False)

    print("\nOutput files were written to:")
    print(ca_output)
    print(exact_output)


if __name__ == "__main__":
    main()