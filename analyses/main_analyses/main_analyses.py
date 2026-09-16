import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math 
from pathlib import Path

def highlight_validity_reclassifications(row):
    if row["valid_crosslink_ca"] == True and row["valid_crosslink_exact"] == False: 
        return ["background-color: orange"] * len(row)

    if row["valid_crosslink_ca"] == False and row["valid_crosslink_exact"] == True: 
        return ["background-color: lightblue"] * len(row)

    return [""] * len(row)

def add_coloring_info_column(df):
    conditions = [
        (df["valid_crosslink_ca"] == True) & (df["valid_crosslink_exact"] == False),
        (df["valid_crosslink_ca"] == False) & (df["valid_crosslink_exact"] == True), 
    ]
    choices = ["orange", "lightblue"]
    df = df.copy()
    df["color"] = np.select(conditions, choices, default="none")
    return df

def get_df_excerpt_middle(df):
    n = 23
    start = 157
    middle = df.iloc[start:start+n]

    omitted=pd.DataFrame(
        [["..."] * len(df.columns)],
        columns = df.columns
    )
    display_df = pd.concat(
        [omitted, middle, omitted],
        ignore_index=True
    )
    return display_df

def get_df_beginning_end(df): 
    if len(df) > 28:
        first = df.head(12)
        last = df.tail(12)
        omitted=pd.DataFrame(
            [["..."] * len(df.columns)],
            columns = df.columns
        )
        display_df = pd.concat(
            [first, omitted, last],
            ignore_index=True
        )
    else:
        display_df = df.copy()
    return display_df


def create_reduced_df(df):
    group_columns = [
        "Protein_id1",
        "reactive_atom1_ca",
        "reactive_atom2_ca",
        "alphafold_distance_ca",
        "valid_crosslink_ca",
        "reactive_atom1_exact",
        "reactive_atom2_exact",
        "alphafold_distance_exact",
        "valid_crosslink_exact",
        "relevant_pae_exact",
        "distance_difference",
        "absolute_difference",
    ]
    df["alphafold_distance_ca"] = df["alphafold_distance_ca"].round(2)
    df["alphafold_distance_exact"] = df["alphafold_distance_exact"].round(2)
    df["relevant_pae_exact"] = df["relevant_pae_exact"].round(2)
    df["distance_difference"] = df["distance_difference"].round(2)
    df["absolute_difference"] = df["absolute_difference"].round(2)
    result = (
        df[group_columns]
        .groupby(group_columns, dropna=False)
        .size()
        .reset_index(name="count")
    )
    result = result.sort_values(
        "absolute_difference",
        ascending=False
    )
    
    df_excerpt_middle = get_df_excerpt_middle(result)
    df_beginning_end = get_df_beginning_end(result)

    color_annotated_df_excerpt_middle = add_coloring_info_column(df_excerpt_middle)
    color_annotated_df_beginning_end = add_coloring_info_column(df_beginning_end)

    return color_annotated_df_excerpt_middle, color_annotated_df_beginning_end 


def compare_distances(merged: pd.DataFrame, output_file: Path) -> None:
    """
    Compares AlphaFold distances calculated using CA and exact binding atoms.

    Calculates distance differences and summary statistics, exports the
    results to CSV and XLSX files, and generates a histogram and ecdf.

    :param merged: DataFrame containing merged CA and exact-atom results.
    :param output_file: Path used as the base name for output files.

    :return: None.
    """
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

    #---------------------------------------------------
    # summary statistics 
    #---------------------------------------------------

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

    print("\nDistance difference statistics")
    print("-" * 40)
    for key, value in stats.items():
        if key == "count":
            print(f"{key:>10}: {value:.0f}")
        else:
            print(f"{key:>10}: {value:.2f} Å")


    #---------------------------------------------------
    # histogram 
    #---------------------------------------------------

    histogram_file = output_file.with_name(
        output_file.stem + "_histogram.png"
    )

    plt.figure(figsize=(8, 5))
    plt.rcParams["font.family"]="serif"
    plt.rcParams["font.size"]=16
    plt.hist(
        merged["distance_difference"],
        bins="auto" 
    )
    plt.xlabel("Distance difference (Å)", fontsize=18)
    plt.ylabel("Number of crosslinks", fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    xmin = math.floor(merged["distance_difference"].min())
    xmax = math.ceil(merged["distance_difference"].max())
    plt.xlim(xmin, xmax)
    
    plt.tight_layout()
    plt.savefig(histogram_file, dpi=300)
    plt.close()


    #---------------------------------------------------
    # ecdf 
    #---------------------------------------------------

    ecdf_file = output_file.with_name(
        output_file.stem + "_ecdf.png"
    )

    x = np.sort(merged["absolute_difference"])
    y = np.arange(1, len(x)+1) / len(x)*100

    plt.figure(figsize=(8, 5))
    plt.rcParams["font.family"]="serif"
    plt.rcParams["font.size"]=14
    plt.step(x, y, where="post")
    plt.xlabel("Absolute distance change (Å)", fontsize=16)
    plt.ylabel("Cumulative percentage of crosslinks (%)", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlim(
        math.floor(merged["absolute_difference"].min()), 
        math.ceil(merged["absolute_difference"].max())
    )
    plt.ylim(0, 100)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(ecdf_file, dpi=300)
    plt.close()


    #---------------------------------------------------
    # full table highlighted .xlsx
    #---------------------------------------------------

    styled = merged.style.apply(highlight_validity_reclassifications, axis=1)
    output_df_file = output_file.with_name(
        output_file.stem + "_df.xlsx"
    )
    if not merged.empty:
        styled.to_excel(
            output_df_file,
            index=False
        )


    #---------------------------------------------------
    # excerpts with color information column .csv
    #---------------------------------------------------

    df_excerpt_middle, df_beginning_end = create_reduced_df(merged)

    output_reduced_df_file = output_file.with_name(
        output_file.stem + "_df_reduced.csv"
    )
    df_beginning_end.to_csv(output_reduced_df_file)

    output_excerpt_df_file = output_file.with_name(
        output_file.stem + "_df_excerpt.csv"
    )
    df_excerpt_middle.to_csv(output_excerpt_df_file)


    print("\nOutput files were written to:")
    print(stats_file)
    print(histogram_file)
    print(ecdf_file)
    print(output_df_file)
    print(output_reduced_df_file)
    print(output_excerpt_df_file)


def compare_validation(merged: pd.DataFrame, output_file: Path) -> None:
    """
    Compares CA-based and exact-atom crosslink validation results.

    Calculates a transition confusion matrix, a stacked bar chart and summarizes changes between the two
    validation results.

    :param merged: DataFrame containing merged CA and exact-atom results.
    :param output_file: Path used as the base name for output files.

    :return: None.
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


    #---------------------------------------------------
    # transition confusion matrix 
    #---------------------------------------------------

    transition_matrix_counts = pd.crosstab(
        pd.Categorical(
            ca_valid, 
            categories=[True, False]
        ),
        pd.Categorical(
            exact_valid, 
            categories=[True, False]
        ),
        margins=True,
        dropna=False,
    )
    transition_matrix_counts.index = [r"$C\alpha$ Valid", r"$C\alpha$ Invalid", "Exact Total"]
    transition_matrix_counts.columns = ["Exact Valid", "Exact Invalid", r"$C\alpha$ Total"]

    transition_matrix_percentages = pd.crosstab(
        pd.Categorical(
            ca_valid, 
            categories=[True, False]
        ),
        pd.Categorical(
            exact_valid, 
            categories=[True, False]
        ),
        normalize="index",
        dropna=False, 
    ) * 100
    transition_matrix_percentages.index = [r"$C\alpha$ Valid", r"$C\alpha$ Invalid"]
    transition_matrix_percentages.columns = ["Exact Valid", "Exact Invalid"]

    total = len(merged)
    ca_total_percentages = transition_matrix_counts.loc[[r"$C\alpha$ Valid", r"$C\alpha$ Invalid"], r"$C\alpha$ Total"] / total * 100
    exact_total_percentages = transition_matrix_counts.loc["Exact Total", ["Exact Valid", "Exact Invalid"]] / total * 100
    
    transition_matrix_whole = transition_matrix_counts.astype(str)
    for row in [r"$C\alpha$ Valid", r"$C\alpha$ Invalid"]:
        for column in ["Exact Valid", "Exact Invalid"]:
            count = transition_matrix_counts.loc[row, column]
            percentage = transition_matrix_percentages.loc[row, column]
            transition_matrix_whole.loc[row, column] = (
                f"{count} ({percentage:.2f}\\%)"
            )

    for row in [r"$C\alpha$ Valid", r"$C\alpha$ Invalid"]:
        count = transition_matrix_counts.loc[row, r"$C\alpha$ Total"]
        percentage = ca_total_percentages.loc[row]
        transition_matrix_whole.loc[row, r"$C\alpha$ Total"] = (
            f"{count} ({percentage:.2f}\\%)"
        )
    
    for column in ["Exact Valid", "Exact Invalid"]:
        count = transition_matrix_counts.loc["Exact Total", column]
        percentage = exact_total_percentages.loc[column]
        transition_matrix_whole.loc["Exact Total", column] = (
            f"{count} ({percentage:.2f}\\%)"
        )

    transition_matrix_whole.loc["Exact Total", r"$C\alpha$ Total"] = f"{total} (100.00\\%)"

    transition_matrix_file = output_file.with_name(
        output_file.stem + "_transition_matrix.csv"
    )
    transition_matrix_whole.to_csv(transition_matrix_file)

    transition_matrix_latex_file = output_file.with_name(
        output_file.stem + "_transition_matrix_latex.tex"
    )
    transition_matrix_latex = transition_matrix_whole.to_latex(
        escape=False
    )
    with open(transition_matrix_latex_file, "w") as f:
        f.write(transition_matrix_latex)

    print("\nTransition Confusion Matrix:")
    print("-" * 40)
    print(transition_matrix_whole)


    #---------------------------------------------------
    # stacked bar chart
    #---------------------------------------------------

    unchanged = (
        ((ca_valid) & (exact_valid)) |
        ((~ca_valid) & (~exact_valid))
    ).sum()

    originally_invalid = (~ca_valid).sum()
    originally_valid = (ca_valid).sum() 

    finally_invalid = (~exact_valid).sum()
    finally_valid = (exact_valid).sum() 

    invalid_to_valid = (
        (~ca_valid) & (exact_valid)
    ).sum()

    valid_to_invalid = (
        (ca_valid) & (~exact_valid)
    ).sum()

    stacked_bar_chart_file = output_file.with_name(
        output_file.stem + "_stacked_bar_chart.png"
    )

    approaches = ["C\N{GREEK SMALL LETTER ALPHA} atoms", "exact reactive atoms"]
    valid = [originally_valid/ len(merged) * 100, finally_valid/ len(merged) * 100]
    invalid = [originally_invalid/ len(merged) * 100, finally_invalid/ len(merged) * 100]

    plt.bar(approaches, valid, label="Valid")
    plt.bar(approaches, invalid, bottom=valid, label="Invalid")

    plt.ylabel("Crosslinks (%)")
    plt.ylim(0, 100)
    plt.legend()
    plt.tight_layout()
    plt.savefig(stacked_bar_chart_file, dpi=300)
    plt.close()


    #---------------------------------------------------
    # summary statistics
    #---------------------------------------------------

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
        {
            "category": "Finally invalid (Exact)",
            "count": finally_invalid,
            "percentage": finally_invalid / len(merged) * 100,
        },
        {
            "category": "Finally valid (Exact)",
            "count": finally_valid,
            "percentage": finally_valid / len(merged) * 100,
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

    print("\nSummary Statistics")
    print("-" * 40)

    for _, row in validation_changes_df.iterrows():
        print(
            f"{row['category']:>25}: "
            f"{row['count']:5.0f} "
            f"({row['percentage']:5.1f}%)"
        )

    
    print("\nOutput files were written to:")
    print(transition_matrix_file)
    print(transition_matrix_latex_file)
    print(stacked_bar_chart_file)
    print(validation_changes_file)


def merge_dfs(df_ca: pd.DataFrame, df_exact: pd.DataFrame) -> pd.DataFrame:
    """
    Merges CA-based and exact-atom DataFrames using shared crosslink data.

    Columns containing distances, reactive atoms, and validation results
    are excluded from the merge keys so they can be compared separately.

    :param df_ca: DataFrame containing CA-based crosslink results.
    :param df_exact: DataFrame containing exact-atom crosslink results.

    :return: Merged DataFrame with CA and exact-atom columns distinguished
        by suffixes.
    """
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
    """
    Runs the crosslink distance and validation analysis.

    Allows the user to select an input folder, analysis, and validation
    criterion before loading the corresponding CA and exact-atom results.

    :return: None.
    """
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

    script_dir = Path(__file__).parent
    input_dir = script_dir / "input_main_analyses"
    output_dir = script_dir / "output_main_analyses"

    output_dir.mkdir(exist_ok=True)

    print("=== start analysis ===\n")

    folders = sorted(f for f in input_dir.iterdir() if f.is_dir())
    if not folders:
        raise FileNotFoundError(f"No input folders found in {input_dir}")

    print("Available input folders:")
    for i, folder in enumerate(folders, start=1):
        print(f"{i}: {folder.name}")
    print()

    while True:
        try:
            choice = int(input(f"Select folder [1-{len(folders)}]: "))
            if 1 <= choice <= len(folders):
                selected_folder = folders[choice - 1]
                break
        except ValueError:
            pass
        print("Invalid selection.")

    folder_name = selected_folder.name

    for file in selected_folder.iterdir():
        if not file.is_file():
            continue

        stem = file.stem  
        
        for key in expected_files:
            if stem.endswith(f"_{key}"):
                expected_files[key] = file
                break

    analysis_names = list(FUNCTIONS.keys())

    print("\nAvailable analyses:")
    for i, name in enumerate(analysis_names, start=1):
        print(f"{i}: {name}")
    print()

    while True:
        try:
            choice = int(input(f"Select analysis [1-{len(analysis_names)}]: "))
            if 1 <= choice <= len(analysis_names):
                analysis_name = analysis_names[choice - 1]
                analysis = FUNCTIONS[analysis_name]
                break
        except ValueError:
            pass

        print("Invalid selection.")

    print("\nAvailable validation criteria:")
    for i, criterion in enumerate(validation_criteria, start=1):
        print(f"{i}: {criterion}")
    print()

    while True:
        try:
            choice = int(input(f"Select validation criterion [1-{len(validation_criteria)}]: "))
            if 1 <= choice <= len(validation_criteria):
                validation_criterion = validation_criteria[choice - 1]
                break
        except ValueError:
            pass

        print("Invalid selection.")

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

    duplicate_columns_ca = df_ca.columns.difference(["id"])
    df_ca = df_ca.drop_duplicates(subset=duplicate_columns_ca)
    duplicate_columns_exact = df_exact.columns.difference(["id"])
    df_exact = df_exact.drop_duplicates(subset=duplicate_columns_exact)

    merged = merge_dfs(df_ca, df_exact)

    output_folder = output_dir / folder_name
    output_folder.mkdir(exist_ok=True)
    output_file = output_folder / f"{analysis_name}_{folder_name}_{validation_criterion}.csv"

    analysis(merged, output_file)


FUNCTIONS = {
    "compare_distances": compare_distances,
    "compare_validation": compare_validation,
}

if __name__ == "__main__":
    main()
