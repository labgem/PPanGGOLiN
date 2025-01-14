
from pathlib import Path


def write_fasta_for_intergenomic_sequences(intergenomic, output_dir: Path):
    """Write the intergenomic sequences to FASTA files for each organism."""
    # Ensure the directory exists
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    # Iterate over the sequences for each organism
    for organism_name, sequences in intergenomic.extract_sequences().items():
        # Create a new FASTA file for each organism in the output directory
        fasta_file_path = output_dir / f"{organism_name}_intergenomic_sequences.fasta"

        with open(fasta_file_path, 'w') as fasta_file:
            for seq_data in sequences:
                gene_pair = seq_data['gene_pair']
                seq_id = seq_data['sequence_id']
                sequence = seq_data['sequence']
                overlap_type = seq_data['overlap_type']

                # Write each sequence in FASTA format
                header = f">{organism_name}_gene_pair_{gene_pair[0]}_{gene_pair[1]}_seqid_{seq_id}_overlap_{overlap_type}"
                fasta_file.write(f"{header}\n")
                fasta_file.write(f"{sequence}\n")

        print(f"FASTA file for {organism_name} written to {fasta_file_path}")
