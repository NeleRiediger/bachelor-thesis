import matplotlib.pyplot as plt
from pathlib import Path

def create_stacked_bar_chart(crosslinkers, valid, invalid, output_file):
    plt.bar(crosslinkers, valid, label="Valid")
    plt.bar(crosslinkers, invalid, bottom=valid, label="Invalid")

    for i, (val, inval) in enumerate(zip(valid, invalid)):
        plt.text(
            i,
            val/2,
            f"{val:.2f}%",
            ha="center",
            va="center"
        )
        if inval < 5:
            plt.text(
                i,
                val+inval/1,
                f"{inval:.2f}%",
                ha="center",
                va="bottom"
            )
        else: 
            plt.text(
                i,
                val+inval/2,
                f"{inval:.2f}%",
                ha="center",
                va="center"
            )

    plt.ylabel("Crosslinks (%)")
    plt.ylim(0, 100)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def main():
    project_dir = Path(__file__).parent
    output_dir = project_dir / "overall_figures"

    output_dir.mkdir(exist_ok=True)

    crosslinkers = ["Amine-Reactive", "Carboxyl-Reactive", "Disulfide Bonds"]

    validation_comparison_ca_file = output_dir / "validity_comparison_ca.png"

    amine_reactive_valid_ca = 28.42
    amine_reactive_invalid_ca = 71.58

    carboxyl_reactive_valid_ca = 41.49
    carboxyl_reactive_invalid_ca = 58.51

    disulfide_bonds_valid_ca = 10.40
    disulfide_bonds_invalid_ca = 89.60

    valid_ca = [amine_reactive_valid_ca, carboxyl_reactive_valid_ca, disulfide_bonds_valid_ca]
    invalid_ca = [amine_reactive_invalid_ca, carboxyl_reactive_invalid_ca, disulfide_bonds_invalid_ca]

    create_stacked_bar_chart(
        crosslinkers=crosslinkers,
        valid=valid_ca,
        invalid=invalid_ca,
        output_file=validation_comparison_ca_file,
    )

    #crosslinkers = ["C\N{GREEK SMALL LETTER ALPHA} atoms", "exact reactive atoms"]
    #valid = [originally_valid/ len(merged) * 100, finally_valid/ len(merged) * 100]
    #invalid = [originally_invalid/ len(merged) * 100, finally_invalid/ len(merged) * 100]

    validation_comparison_exact_file = output_dir / "validity_comparison_exact.png"

    amine_reactive_valid_exact = 23.30
    amine_reactive_invalid_exact = 76.70

    carboxyl_reactive_valid_exact = 37.23
    carboxyl_reactive_invalid_exact = 62.77

    disulfide_bonds_valid_exact = 97.45
    disulfide_bonds_invalid_exact = 2.55

    valid_exact = [amine_reactive_valid_exact, carboxyl_reactive_valid_exact, disulfide_bonds_valid_exact ]
    invalid_exact = [amine_reactive_invalid_exact, carboxyl_reactive_invalid_exact, disulfide_bonds_invalid_exact]

    create_stacked_bar_chart(
        crosslinkers=crosslinkers,
        valid=valid_exact,
        invalid=invalid_exact,
        output_file=validation_comparison_exact_file,
    )

    print("\nOutput files were written to:")
    print(validation_comparison_ca_file)
    print(validation_comparison_exact_file)


if __name__ == "__main__":
    main()