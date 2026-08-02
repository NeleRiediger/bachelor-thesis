import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def compare_distances(merged: pd.DataFrame, output_file: Path) -> None:

    merged["distance_difference"] = (
        merged["alphafold_distance_exact"]
        - merged["alphafold_distance_ca"]
    )

    merged["absolute_difference"] = (
        merged["distance_difference"].abs()
    )

    merged = merged.sort_values(
        "absolute_difference",
        ascending=False
    )

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
    stats_file = output_file.with_name(
    output_file.stem + "_stats.csv"
    )

    pd.DataFrame([stats]).to_csv(
        stats_file,
        index=False
    )


    abs_diff = merged["absolute_difference"]
    thresholds = [0, 0.5, 1, 2, 5]
    absolute_changes = []

    for t in thresholds:
        if t == 0:
            count = (abs_diff == 0).sum()
            label = "No change"
        else:
            count = (abs_diff < t).sum()
            label = f"< {t} Å"

        absolute_changes.append(
            {
                "threshold": label,
                "count": count,
                "percentage": count / len(abs_diff) * 100,
            }
        )

    absolute_changes_df = pd.DataFrame(absolute_changes)
    absolute_changes_file = output_file.with_name(
        output_file.stem + "_abs_stats.csv"
    )

    absolute_changes_df.to_csv(
        absolute_changes_file,
        index=False
    )

    print("\nDistance difference statistics")
    print("-" * 40)

    for key, value in stats.items():
        if key == "count":
            print(f"{key:>10}: {value:.0f}")
        else:
            print(f"{key:>10}: {value:.2f} Å")


    print("\nAbsolute distance changes")
    print("-" * 40)

    for _, row in absolute_changes_df.iterrows():
        print(
            f"{row['threshold']:>10}: "
            f"{row['count']:5.0f} "
            f"({row['percentage']:5.1f}%)"
        )


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


    if not merged.empty:
        merged.to_csv(
            output_file.stem + "_df.csv",
            index=False
        )

    print("\nOutput files were written to:")
    print(output_file)
    print(stats_file)
    print(absolute_changes_file)
    print(histogram_file)
    print(boxplot_file)


def compare_validation(merged: pd.DataFrame, output_file: Path) -> None:
    """
    Compare CA-based and exact-atom validation results.
    """

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

    confusion_matrix_file = output_file.with_name(
        output_file.stem + "_confusion_matrix.csv"
    )
    confusion_matrix.to_csv(confusion_matrix_file)

    unchanged = (
        ((ca_valid) & (exact_valid)) |
        ((~ca_valid) & (~exact_valid))
    ).sum()

    originally_invalid = (~ca_valid).sum()
    originally_valid = (ca_valid).sum() 

    invalid_to_valid = (
        (~ca_valid) & (exact_valid)
    ).sum()

    valid_to_invalid = (
        (ca_valid) & (~exact_valid)
    ).sum()

    validation_changes = [
        {
            "category": "Unchanged",
            "count": unchanged,
            "percentage": unchanged / len(merged) * 100,
        },
        {
            "category": "Originally invalid (CA)",
            "count": originally_invalid,
            "percentage": originally_invalid / len(merged) * 100,
        },
        {
            "category": "Originally valid (CA)",
            "count": originally_valid,
            "percentage": originally_valid / len(merged) * 100,
        },
        {
            "category": "Invalid → Valid",
            "count": invalid_to_valid,
            "percentage": invalid_to_valid / originally_invalid * 100
            if originally_invalid > 0
            else 0,
        },
        {
            "category": "Valid → Invalid",
            "count": valid_to_invalid,
            "percentage": valid_to_invalid / originally_valid * 100
            if originally_valid > 0
            else 0,
        },
    ]

    validation_changes_df = pd.DataFrame(validation_changes)

    validation_changes_file = output_file.with_name(
        output_file.stem + "_stats.csv"
    )

    validation_changes_df.to_csv(
        validation_changes_file,
        index=False
    )

    print("\nValidation changes")
    print("-" * 40)

    for _, row in validation_changes_df.iterrows():
        print(
            f"{row['category']:>25}: "
            f"{row['count']:5.0f} "
            f"({row['percentage']:5.1f}%)"
        )

    print("\nOutput files were written to:")
    print(confusion_matrix_file)
    print(validation_changes_file)
     

"""
def compare_all_validation_methods(files: dict, output_file: Path) -> None:
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
"""  


"""
def analyse_by_binding_site(files: dict, output_file: Path) -> None: 
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
        print("Missing files:")
        for key in missing:
            print(f" - {key}")
        return None

    df_ca = pd.read_csv(files[ca_key])
    df_exact = pd.read_csv(files[exact_key])

    merged = merge_dfs(df_ca, df_exact)

    merged["distance_difference"] = (
        merged["alphafold_distance_exact"]
        - merged["alphafold_distance_ca"]
    )

    merged["validation_change"] = "unchanged"

    merged.loc[
        (~merged["valid_crosslink_ca"]) &
        (merged["valid_crosslink_exact"]),
        "validation_change"
    ] = "became_valid"

    merged.loc[
        (merged["valid_crosslink_ca"]) &
        (~merged["valid_crosslink_exact"]),
        "validation_change"
    ] = "became_invalid"


    merged["atom_pair"] = (
        merged["reactive_atom1_exact"].astype(str)
        + "-"
        + merged["reactive_atom2_exact"].astype(str)
    )


    analysis = (
        merged
        .groupby("atom_pair")
        .agg(
            count=("atom_pair", "size"),
            mean_distance_change=(
                "distance_difference",
                "mean"
            ),
            median_distance_change=(
                "distance_difference",
                "median"
            ),
            became_valid=(
                "validation_change",
                lambda x: (x == "became_valid").sum()
            ),
            became_invalid=(
                "validation_change",
                lambda x: (x == "became_invalid").sum()
            ),
            unchanged=(
                "validation_change",
                lambda x: (x == "unchanged").sum()
            ),
        )
        .reset_index()
    )

    analysis["improvement_rate"] = (
        analysis["became_valid"]
        /
        analysis["count"]
    )

    analysis = analysis.sort_values(
        [
            "became_valid",
            "improvement_rate"
        ],
        ascending=False
    )

    result_file = output_file.with_name(
        output_file.stem + "_binding_atom_analysis.csv"
    )

    analysis.to_csv(
        result_file,
        index=False
    )

    print("\nBinding atom analysis")
    print("=" * 90)
    print(
        analysis.to_string(index=False)
    )

    print(f"\nOutput written to:")
    print(result_file)
"""


def merge_dfs(df_ca: pd.DataFrame, df_exact: pd.DataFrame) -> pd.DataFrame:
    merge_columns = [
        "Protein1",
        "Protein2",
        "Protein_id1",
        "Protein_id2",
        "Is_intra_crosslink",
        "Crosslinker",
        "Peptide1",
        "Peptide2",
        "1_based_CL_position_within_peptide1",
        "1_based_CL_position_within_peptide2",
        "Q_value",
        "Chain_id1",
        "Chain_id2",
        "1_based_crosslinker_position1",
        "1_based_crosslinker_position2",
    ]
    # especially "alphafold_distance", "reactive_atom1", "reactive_atom2", "valid_crosslink" can differ

    merged = df_ca.merge(
        df_exact,
        on=merge_columns,
        suffixes=("_ca", "_exact"),
        how="inner",
    )

    return merged 


def main(): 
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

    validation_criteria = [
        "manual",
        "maxPAE",
        "minPAE",
        "plDDT",
    ]

    project_dir = Path(__file__).parent
    input_dir = project_dir / "input_results"
    output_dir = project_dir / "output_results"

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
        return

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

    print("\nAvailable validation criteria:")
    for criterion in validation_criteria:
        print(f" - {criterion}")
    print()

    validation_criterion = input("Validation criterion: ")

    if validation_criterion not in validation_criteria:
        print(f"Unknown validation criterion: {validation_criterion}")
        return None

    ca_key = f"{validation_criterion}_CA"
    exact_key = f"{validation_criterion}_exact"

    missing = [
        key for key in (ca_key, exact_key)
        if expected_files.get(key) is None
    ]
    if missing:
        print("Error: Missing required files:")
        for key in missing:
            print(f" - *_{key}")
        return None
    
    df_ca = pd.read_csv(expected_files[ca_key])
    df_exact = pd.read_csv(expected_files[exact_key])

    merged = merge_dfs(df_ca, df_exact)

    output_folder = output_dir / folder_name
    output_folder.mkdir(exist_ok=True)
    output_file = output_folder / f"{analysis_name}_{folder_name}_{validation_criterion}.csv"

    analysis(merged, output_file)


FUNCTIONS = {
    "compare_distances": compare_distances,
    "compare_validation": compare_validation,
    #"compare_all_validation_methods": compare_all_validation_methods,
    #"analyse_by_binding_site": analyse_by_binding_site,
}

if __name__ == "__main__":
    main()
