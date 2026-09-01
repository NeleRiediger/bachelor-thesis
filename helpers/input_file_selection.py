from pathlib import Path

def choose_file(input_dir: Path, file_type: str) -> Path:
    if not file_type.startswith("."):
        file_type = f".{file_type}"

    files = sorted(
        list(input_dir.glob(f"*{file_type}"))
    )

    if not files:
        raise FileNotFoundError(f"No input file available in {input_dir}")

    print("Available input:")
    for i, file in enumerate(files, start=1):
        print(f"{i}: {file.name}")
    print()

    while True:
        choice = input(f"Select file [1-{len(files)}]: ")

        if choice == ".":
            return None

        try:
            choice = int(choice)
        except ValueError:
            print("Please enter a number.")
            continue

        if 1 <= choice <= len(files):
            return files[choice - 1]

        print("Invalid selection.")