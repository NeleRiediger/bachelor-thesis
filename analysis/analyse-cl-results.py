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

    output_file = output_dir / f"{analysis_name}_{folder_name}.csv"

    analysis(expected_files, output_file)


def merge_dfs(df_ca: pd.DataFrame, df_exact: pd.DataFrame) -> pd.DataFrame:
    merge_columns = [
        col
        for col in df_ca.columns
        if col not in {
            "alphafold_distance",
            "reactive_atom1",
            "reactive_atom2",
            "valid_crosslink",
        }
    ]

    merged = df_ca.merge(
        df_exact,
        on=merge_columns,
        suffixes=("_ca", "_exact"),
        how="inner",
    )

    return merged 


def compare_distances(files: dict, output_file: Path) -> None: 
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
    
    df_ca = pd.read_csv(files[ca_key])
    df_exact = pd.read_csv(files[exact_key])

    merged = merge_dfs(df_ca, df_exact)

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


def compare_validation(files: dict, output_file: Path) -> None:
    """
    Compare CA-based and exact-atom validation results.

    For each selected validation criterion:
    - prints confusion matrix
    - calculates invalid -> valid improvement
    """

    available_criteria = sorted(
        {
            key.replace("_CA", "").replace("_exact", "")
            for key in files.keys()
        }
    )

    validation_criteria = input(
        f"Validation criteria {available_criteria} (comma separated): "
    ).strip()

    selected_criteria = [
        criterion.strip()
        for criterion in validation_criteria.split(",")
    ]

    for criterion in selected_criteria:
        file_ca = files.get(f"{criterion}_CA")
        file_exact = files.get(f"{criterion}_exact")

        if file_ca is None or file_exact is None:
            print(
                f"\nSkipping {criterion}: missing CA or exact file"
            )
            continue

        print("\n" + "=" * 60)
        print(f"Validation criterion: {criterion}")
        print("=" * 60)

        df_ca = pd.read_csv(file_ca)
        df_exact = pd.read_csv(file_exact)

        merged = merge_dfs(df_ca, df_exact)

        ca_valid = (
            merged["valid_crosslink_ca"]
            .astype(str)
            .str.lower()
            .eq("true")
        )

        exact_valid = (
            merged["valid_crosslink_exact"]
            .astype(str)
            .str.lower()
            .eq("true")
        )

        confusion_matrix = pd.crosstab(
            ca_valid,
            exact_valid
        )

        confusion_matrix = confusion_matrix.reindex(
            index=[True, False],
            columns=[True, False],
            fill_value=0
        )

        confusion_matrix.index = ["valid", "invalid"]
        confusion_matrix.columns = ["valid", "invalid"]

        confusion_matrix.index.name = "CA"
        confusion_matrix.columns.name = "Exact"

        print("\nConfusion Matrix:")
        print(confusion_matrix)

        invalid_to_valid = (
            (~ca_valid) & (exact_valid)
        ).sum()

        originally_invalid = (~ca_valid).sum()

        improvement = (
            invalid_to_valid / originally_invalid
            if originally_invalid > 0
            else 0
        )

        print(
            f"\nInvalid → Valid improvement: "
            f"{invalid_to_valid}/{originally_invalid} "
            f"({improvement:.2%})"
        )

        valid_to_invalid = (
            (ca_valid) & (~exact_valid)
        ).sum()

        originally_valid = (ca_valid).sum()

        deterioration = (
            valid_to_invalid / originally_valid
            if originally_valid > 0
            else 0
        )

        print(
            f"\nValid → Invalid deterioration: "
            f"{valid_to_invalid}/{originally_valid} "
            f"({deterioration:.2%})"
        )

def compare_all_validation_methods(files: dict, output_file: Path) -> pd.DataFrame:
    results = []

    for criterion in VALIDATION_CRITERIA:

        file_ca = files.get(f"{criterion}_CA")
        file_exact = files.get(f"{criterion}_exact")

        if file_ca is None or file_exact is None:
            print(f"Skipping {criterion}: missing files.")
            continue

        df_ca = pd.read_csv(file_ca)
        df_exact = pd.read_csv(file_exact)

        merged = merge_dfs(df_ca, df_exact)

        ca_valid = merged["valid_crosslink_ca"]
        exact_valid = merged["valid_crosslink_exact"]

        remained_valid = (ca_valid & exact_valid).sum()
        remained_invalid = ((~ca_valid) & (~exact_valid)).sum()
        became_valid = ((~ca_valid) & exact_valid).sum()
        became_invalid = (ca_valid & (~exact_valid)).sum()

        results.append(
            {
                "Method": criterion,
                "Remained valid": int(remained_valid),
                "Remained invalid": int(remained_invalid),
                "Became valid": int(became_valid),
                "Became invalid": int(became_invalid),
                "Net gain": int(became_valid - became_invalid),
                "Improvement": f"{became_valid / max(1, (~ca_valid).sum()):.1%}",
            }
        )

    result = pd.DataFrame(results)

    print("\nValidation changes by method")
    print("=" * 90)
    print(result.to_string(index=False))

    if result is not None and not result.empty:
        result.to_csv(output_file, index=False)

    print(f"Output file was written to: ")
    print(output_file)
    

FUNCTIONS = {
    "compare_distances": compare_distances,
    "compare_validation": compare_validation,
    "compare_all_validation_methods": compare_all_validation_methods,
}

if __name__ == "__main__":
    main()
