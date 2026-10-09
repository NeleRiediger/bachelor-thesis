from pathlib import Path
import pandas as pd
from scipy import stats 
import numpy as np

script_dir = Path(__file__).parent
project_dir = script_dir.parent
input_dir = project_dir / "main_analyses" / "input_main_analyses"
output_dir = script_dir / "statistical_test_results"

output_dir.mkdir(exist_ok=True)

#def welch_t_test(dfs, output_file):
def welch_t_test(dfs):
    data1 = dfs[0]["distance_difference"].dropna()
    data2 = dfs[1]["distance_difference"].dropna()

    result = stats.ttest_ind(
        data1,
        data2,
        equal_var=False,
        alternative="two-sided"
    )

    n1, n2 = len(data1), len(data2)
    mean_diff = data1.mean()-data2.mean()

    a = data1.var(ddof=1) / n1
    b = data2.var(ddof=1) / n2

    se = np.sqrt(a+b)

    df_welch = (a+b)**2 / (
        a**2 / (n1-1) + b**2 / (n2-1)
    )

    critical_t = stats.t.ppf(0.975, df_welch)

    ci_low = mean_diff - critical_t * se
    ci_high = mean_diff + critical_t * se
    

    print("Welch-t-Test")
    print(f"Stichprobengröße 1: {len(data1)}")
    print(f"Stichprobengröße 2: {len(data2)}")
    print(f"Mittelwert 1: {data1.mean():.2f} Å")
    print(f"Mittelwert 2: {data2.mean():.2f} Å")
    print(f"Mean difference: {mean_diff:.2f} Å")
    print(f"t-Statistik: {result.statistic:.3f}")
    print(f"p-Wert: {result.pvalue:.3f}")
    print(f"95% CI: [{ci_low:.2f}, {ci_high:.2f}] Å")

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


def main():
    selected_dirs = []

    folders = sorted(f for f in input_dir.iterdir() if f.is_dir())
    
    if len(folders)<2:
        raise FileNotFoundError(f"At least two input folders are required in {input_dir}")
    
    for selection_number in range(1, 3):
        print(f"\nSelect dataset {selection_number}:")
        print("Available input folders:")
        for i, folder in enumerate(folders, start=1):
            marker = "(already selected)" if folder in selected_dirs else ""
            print(f"{i}: {folder.name}{marker}")
        print()

        while True:
            choice = input(f"Select folder [1-{len(folders)}]: ")
            
            try:
                choice = int(choice)
                if not 1 <= choice <= len(folders):
                    print("Invalid selection.")
                    continue 

                selected_folder = folders[choice - 1]

                if selected_folder in selected_dirs:
                    print("This folder has already been selected.")
                    continue 
                   
                selected_dirs.append(selected_folder)
                break 
            
            except ValueError:
                print("Invalid selection.")

    output_file = output_dir / "results.png"

    if not selected_dirs: 
        print("No directories selected.")
        return 
    
    try: 
        dfs = prepare_all_data(selected_dirs)
        #welch_t_test(dfs, output_file)
        welch_t_test(dfs)

    except ValueError as error: 
        print(f"\nError: {error}")
        return
    
    print("\nOutput files were written to:")
    print(output_file)
    

if __name__ == "__main__":
    main()