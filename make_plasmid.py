from plasmid_builder.pipeline import run_pipeline


def main():
    run_pipeline(
        input_fasta="Input.Fa",
        design_file="Design.txt",
        output_fasta="Output.Fa",
        output_genbank="Output.gb"
    )


if __name__ == "__main__":
    main()
