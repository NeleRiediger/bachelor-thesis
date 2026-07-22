import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

VALIDATION_CRITERIA = {
        "manual",
        "maxPAE",
        "minPAE",
        "plDDT",
    }

def main(): 
    project_dir = Path(__file__).parent
    input_dir = project_dir / "input"
    output_dir = project_dir / "output"

    output_dir.mkdir(exist_ok=True)

    print("=== start analysis ===\n")

    """
    print("Available input data:")
    for file in sorted(input_dir.glob("*.csv")):
        print(f" - {file.name}")
    print()

    input_file_ca_binding_site = input("CA binding site data: ")
    file_CA = input_dir / input_file_ca_binding_site
    if not file_CA.exists():
        print(f"Error: {input_file_ca_binding_site} not found.")
        return
    
    input_file_exact_binding_site = input("exact binding site data: ")
    file_exact = input_dir / input_file_exact_binding_site
    if not file_exact.exists():
        print(f"Error: {input_file_exact_binding_site} not found.")
        return
    """

    print("Available input folders:")
    folders = sorted([f for f in input_dir.iterdir() if f.is_dir()])
    for folder in folders:
        print(f" - {folder.name}")
    print()

    folder_name = input("Folder: ")
    selected_folder = input_dir / folder_name
    if not selected_folder.exists() or not selected_folder.is_dir():
        print(f"Error: Folder '{folder_name}' not found.")
        exit()

    expected_files = {
        "manual_CA": None,
        "manual_exact": None,
        "maxPAE_CA": None,
        "maxPAE_exact": None,
        "minPAE_CA": None,
        "minPAE_exact": None,
        "plDDT_CA": None,
        "plDDT_exact": None,
    }

    for file in selected_folder.iterdir():
        if not file.is_file():
            continue

        stem = file.stem  
        
        for key in expected_files:
            if stem.endswith(f"_{key}"):
                expected_files[key] = file
                break

    print("\nAvailable analyses:")
    for name in FUNCTIONS:
        print(f" - {name}")
    print()

    analysis_name = input("Analysis: ")
    if analysis_name not in FUNCTIONS:
        print(f"Error: '{analysis_name}' is not a valid analysis.")
        return
    analysis = FUNCTIONS[analysis_name]

    output_suffix = input("Output suffix: ")
    output_file = output_dir / f"{analysis_name}_{output_suffix}.csv"

    """
    df_CA = pd.read_csv(file_CA)
    df_exact = pd.read_csv(file_exact)

    result_df, result_value = analysis(df_CA, df_exact)
    """

    result_df, result_dict = analysis(expected_files, output_file)

    """
    if result_value: 
        print(f"Result value: {result_value:.2f} Å")

    if result_df is not None and not result_df.empty:
        result_df.to_csv(output_file, index=False)

        print(f"\nDone! Output file was written to: ")
        print(output_file)
    """


"""
def average_cl_length_diff(df_CA: pd.DataFrame, df_exact: pd.DataFrame) -> pd.DataFrame:
    # Alle gemeinsamen Spalten außer alphafold_distance zum Matchen verwenden
    merge_columns = [
        col
        for col in df_CA.columns
            if col not in {
            "alphafold_distance",
            "reactive_atom1",
            "reactive_atom2",
        }
    ]

    merged = df_CA.merge(
        df_exact,
        on=merge_columns,
        suffixes=("_ca", "_exact"),
        how="inner",
    )

    merged["distance_difference"] = (
        merged["alphafold_distance_exact"]
        - merged["alphafold_distance_ca"]
    )

    average_difference = merged["distance_difference"].mean()

    return merged, average_difference
"""

def compare_distances(files: dict, output_file: Path)-> tuple[pd.DataFrame, dict]: 
    validation_criterion = input(
        f"Validation criterion {sorted(VALIDATION_CRITERIA)}: "
    ).strip()

    if validation_criterion not in VALIDATION_CRITERIA:
        print(f"Unknown validation criterion: {validation_criterion}")
        return None

    ca_key = f"{validation_criterion}_CA"
    exact_key = f"{validation_criterion}_exact"

    missing = [
        key for key in (ca_key, exact_key)
        if files.get(key) is None
    ]
    if missing:
        print("Error: Missing required files:")
        for key in missing:
            print(f" - *_{key}")
        return None
    
    df_CA = pd.read_csv(files[ca_key])
    df_exact = pd.read_csv(files[exact_key])

    merge_columns = [
        col
        for col in df_CA.columns
        if col not in {
            "alphafold_distance",
            "reactive_atom1",
            "reactive_atom2",
        }
    ]

    merged = df_CA.merge(
        df_exact,
        on=merge_columns,
        suffixes=("_ca", "_exact"),
        how="inner",
    )

    merged["distance_difference"] = (
        merged["alphafold_distance_exact"]
        - merged["alphafold_distance_ca"]
    )

    merged["absolute_difference"] = merged["distance_difference"].abs()

    stats = {
        "count": len(merged),
        "mean": merged["distance_difference"].mean(),
        "median": merged["distance_difference"].median(),
        "std": merged["distance_difference"].std(),
        "min": merged["distance_difference"].min(),
        "25%": merged["distance_difference"].quantile(0.25),
        "75%": merged["distance_difference"].quantile(0.75),
        "max": merged["distance_difference"].max(),
    }

    print("\nDistance difference statistics")
    print("-" * 40)
    for key, value in stats.items():
        if key == "count":
            print(f"{key:>8}: {value}")
        else:
            print(f"{key:>8}: {value:.2f} Å")

    abs_diff = merged["absolute_difference"]
    thresholds = [0, 0.5, 1, 2, 5]
    print("\nAbsolute distance changes")
    print("-" * 40)
    for t in thresholds:
        if t == 0:
            count = (abs_diff == 0).sum()
            percent = count / len(abs_diff) * 100
            print(f"No change : {count:5} ({percent:5.1f}%)")
        else:
            count = (abs_diff < t).sum()
            percent = count / len(abs_diff) * 100
            print(f"< {t} Å : {count:5} ({percent:5.1f}%)")

    histogram_file = output_file.with_name(
        output_file.stem + "_histogram.png"
    )
    plt.figure(figsize=(8, 5))
    plt.hist(
        merged["distance_difference"],
        bins=40
    )
    plt.xlabel("Distance difference (Å)")
    plt.ylabel("Number of crosslinks")
    plt.title("Distribution of distance differences")
    plt.tight_layout()
    plt.savefig(histogram_file, dpi=300)
    plt.close()

    boxplot_file = output_file.with_name(
        output_file.stem + "_boxplot.png"
    )
    plt.figure(figsize=(4, 6))
    plt.boxplot(
        merged["distance_difference"],
        orientation="vertical"
    )
    plt.ylabel("Distance difference (Å)")
    plt.title("Distance differences")
    plt.tight_layout()
    plt.savefig(boxplot_file, dpi=300)
    plt.close()

    if merged is not None and not merged.empty:
        merged.to_csv(output_file, index=False)
        
    print(f"Output files were written to: ")
    print(output_file)
    print(histogram_file)
    print(boxplot_file)

    return merged, stats


def compare_validation(files: dict) -> pd.DataFrame:
    return files

def compare_all_validation_methods(files: dict) -> pd.DataFrame:
    missing = [key for key, value in expected_files.items() if value is None]

    if missing:
        print("Error: Missing files:")
        for m in missing:
            print(f" - *_{m}")
        exit()

    file_manual_CA = expected_files["manual_CA"]
    file_manual_exact = expected_files["manual_exact"]
    file_maxPAE_CA = expected_files["maxPAE_CA"]
    file_maxPAE_exact = expected_files["maxPAE_exact"]
    file_minPAE_CA = expected_files["minPAE_CA"]
    file_minPAE_exact = expected_files["minPAE_exact"]
    file_plDDT_CA = expected_files["plDDT_CA"]
    file_plDDT_exact = expected_files["plDDT_exact"]
    return files 

FUNCTIONS = {
    #"average_cl_length_diff": average_cl_length_diff,
    "compare_distances": compare_distances,
    "compare_validation": compare_validation,
    "compare_all_validation_methods": compare_all_validation_methods,
}

if __name__ == "__main__":
    main()
