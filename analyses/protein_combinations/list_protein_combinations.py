from pathlib import Path
import pandas as pd
import sys

script_dir = Path(__file__).resolve().parent
project_dir = script_dir.parent.parent
sys.path.insert(0, str(project_dir))
from helpers.input_file_selection import choose_file


input_dir = script_dir / "input_protein_combinations"
output_dir = script_dir / "output_protein_combinations"

output_dir.mkdir(exist_ok=True)


def rank_protein_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """ 
    Ranks protein pairs by their number of occurrences. 
    
    Protein pairs are treated as unordered, 
    so the order of the two protein IDs does not affect their ranking. 
    Intra- and inter-protein crosslinks are ranked separately. 
    
    :param df: DataFrame containing protein pair information. 
    
    :return: DataFrame containing protein pairs, crosslink type, 
        and occurrence counts. 
    """
    required = {
        "Protein_id1", 
        "Protein_id2",
        "Is_intra_crosslink",
    }
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing column: {missing}")

    pairs = (
        df[[ 
            "Protein_id1", 
            "Protein_id2", 
            "Is_intra_crosslink", 
        ]]
        .copy()
    )

    pairs[["Protein_id1", "Protein_id2"]] = ( 
        pairs[["Protein_id1", "Protein_id2"]] 
        .apply( 
            lambda row: pd.Series(sorted(row)), 
            axis=1, 
        ) 
    )

    ranking = ( 
        pairs.value_counts() 
        .reset_index(name="Count") 
        .sort_values("Count", ascending=False) 
        .reset_index(drop=True) 
    )

    return ranking[
        [
            "Protein_id1", 
            "Protein_id2", 
            "Is_intra_crosslink",
            "Count",
        ]
    ]


def main():
    selected_file = choose_file(input_dir, ".csv")
    df = pd.read_csv(selected_file)

    ranking = rank_protein_pairs(df)
    print(ranking)

    output_file = output_dir / f"{selected_file.stem}_protein_pair_ranking.csv"
    ranking.to_csv(output_file, index=False)

    print(f"\nOutput file was written to: {output_file}")


if __name__ == "__main__":
    main()