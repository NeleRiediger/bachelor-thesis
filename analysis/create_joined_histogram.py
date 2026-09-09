from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import math 

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))

project_dir = Path(__file__).parent
input_dir = project_dir / "input_results"
output_dir = project_dir / "overall_figures"

output_dir.mkdir(exist_ok=True)

def create_separate_histograms(dfs, output_file):
    labels = [
        "Amine-Reactive",
        "Carboxyl-Reactive",
        "Disulfide Bonds", 
        "Zero-Length"
    ]

    #all_values = pd.concat(
    #    [df["distance_difference"].dropna() for df in dfs]
    #)
    #bins = np.histogram_bin_edges(all_values, bins="auto")

    fig, axes = plt.subplots(
        1, 4,
        figsize=(15,5),
    )

    for ax, df, label in zip(axes, dfs, labels):
        values = df["distance_difference"].dropna()
        bins=np.histogram_bin_edges(
            values,
            bins="auto"
        )
        weights = np.ones(len(values)) / len(values)*100

        ax.hist(
            values, 
            bins=bins,
            weights=weights
        )

        xmin=np.floor(values.min())
        xmax=np.ceil(values.max())
        ax.set_xlim(xmin, xmax)

        ax.set_title(label, fontsize=16, fontfamily="serif")
        ax.set_xlabel(
            "Distance difference (Å)",
            fontsize=16,
            fontfamily="serif"
        )
        ax.set_ylabel(
            "Crosslinks (%)",
            fontsize=16,
            fontfamily="serif"
        )

        ax.tick_params(axis="both", labelsize=14)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontfamily("serif")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def create_histogram(dfs, output_file):
    labels = [
        "Amine-Reactive",
        "Carboxyl-Reactive",
        "Disulfide Bonds",
        "Zero-Length"
    ]

    all_values = pd.concat(
        [df["distance_difference"].dropna() for df in dfs]
    )
    bins = np.histogram_bin_edges(all_values, bins="auto")
        
    plt.figure(figsize=(8, 5))

    for df, label in zip(dfs, labels):
        values = df["distance_difference"].dropna()
        weights = np.ones_like(values) / len(values) * 100
        plt.hist(
            values,
            bins=bins, 
            weights=weights,
            alpha=0.5,
            label=label
        )

    plt.xlabel("Distance difference (Å)")
    plt.ylabel("Crosslinks (%)")
    plt.legend()

    #xmin = math.floor(merged["distance_difference"].min())
    #xmax = math.ceil(merged["distance_difference"].max())
    #plt.xlim(xmin, xmax)
    #plt.xticks(range(xmin, xmax+1, 5))
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def prepare_all_data(selected_dirs):
    dfs = []

    for dir in selected_dirs:
        ca_files = list(dir.glob("*_CA.csv"))
        exact_files = list(dir.glob("*_exact.csv"))

        ca_df = pd.read_csv(ca_files[0])
        exact_df = pd.read_csv(exact_files[0])

        duplicate_columns_ca = ca_df.columns.difference(["id"])
        ca_df = ca_df.drop_duplicates(subset=duplicate_columns_ca)
        duplicate_columns_exact = exact_df.columns.difference(["id"])
        exact_df = exact_df.drop_duplicates(subset=duplicate_columns_exact)

        merged = merge_dfs(ca_df, exact_df)

        merged["distance_difference"] = (
            merged["alphafold_distance_exact"]
            - merged["alphafold_distance_ca"]
        )

        dfs.append(merged)

    return dfs

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

    output_file = output_dir / "joined_histograms.png"
    output_file_separate = output_dir / "joined_histograms_separate.png"

    if not selected_dirs: 
        print("No directories selected.")
        return 
    
    try: 
        dfs = prepare_all_data(selected_dirs)
        create_histogram(dfs, output_file)
        create_separate_histograms(dfs, output_file_separate)

    except ValueError as error: 
        print(f"\nError: {error}")
        return

    print("\nOutput file was written to:")
    print(output_file)
    print(output_file_separate)


if __name__ == "__main__":
    main()