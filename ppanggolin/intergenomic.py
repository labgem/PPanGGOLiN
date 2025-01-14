from ppanggolin.genome import Feature, Gene, Organism, Contig
from ppanggolin.edge import Edge


class Intergenomic(Feature):
    def __init__(self, edge: Edge):
        """Initializes the Intergenomic object by associating it with an edge."""

        if not isinstance(edge, Edge):
            raise TypeError("Input must be an Edge object.")

        super().__init__(edge.id)
        self.edge = edge
        self.gene_pairs = edge.gene_pairs
        self.organisms_genes = edge.get_organisms_dict()  # Get organisms and associated gene pairs

        # Dictionaries to manage sequences and gene pair associations
        self.sequence_dict = {}  # Maps sequence IDs to actual sequences
        self.gene_pair_to_seqid = {}  # Maps gene pairs (ID tuple) to sequenc

    def extract_sequences(self):
        """Extract the intergenomic sequences for each organism's gene pair."""
        sequences = {}

        for organism, gene_pairs in self.organisms_genes.items():
            sequences_for_organism = []
            for gene_pair in gene_pairs:
                # Directly use the gene pair (ID tuple) as the key
                gene_pair_key = (gene_pair[0].ID, gene_pair[1].ID)

                if gene_pair_key not in self.gene_pair_to_seqid:
                    # Extract the sequence if not found in dictionary
                    sequence, overlap = self.extract_sequence_between_genes(gene_pair[0], gene_pair[1])

                    # Create a new sequence ID for this sequence
                    seq_id = len(self.sequence_dict) + 1  # Sequence ID is just an increment
                    self.sequence_dict[seq_id] = sequence  # Store sequence with its ID

                    # Mark overlap if applicable
                    if overlap:
                        sequence_type = 'overlap'
                    else:
                        sequence_type = 'non-overlap'

                    # Map gene pair to the new sequence ID
                    self.gene_pair_to_seqid[gene_pair_key] = (seq_id, sequence_type)
                else:
                    # Use existing sequence ID for this gene pair
                    seq_id, sequence_type = self.gene_pair_to_seqid[gene_pair_key]

                # Store the gene pair's sequence ID and the sequence
                sequences_for_organism.append({
                    "gene_pair": (gene_pair[0].ID, gene_pair[1].ID),
                    "sequence_id": seq_id,
                    "sequence": self.sequence_dict[seq_id],
                    "sequence_type": sequence_type
                })
            sequences[organism.name] = sequences_for_organism
        return sequences

    def extract_sequence_between_genes(self, gene1: Gene, gene2: Gene):
        """Extract the sequence between two genes, handling overlap."""
        contig = gene1.contig
        start_pos = gene1.stop
        end_pos = gene2.start

        # Check if there's an overlap
        overlap = False
        if start_pos > end_pos:
            # Genes overlap, adjust the extraction range
            overlap = True
            sequence = contig.dna[start_pos:end_pos]
        else:
            sequence = contig.dna[start_pos:end_pos]  # Extract the non-overlapping sequence

        return sequence, overlap

    def get_gene_pairs_with_sequences(self):
        """Returns a dictionary of organisms with their corresponding gene pairs and sequences."""
        gene_pair_dict = {}

        for organism, gene_pairs in self.organisms_genes.items():
            gene_pair_dict[organism.name] = []
            for gene_pair in gene_pairs:
                # Get the sequence ID directly for this gene pair
                gene_pair_key = (gene_pair[0].ID, gene_pair[1].ID)
                seq_id, overlap_type = self.gene_pair_to_seqid[gene_pair_key]

                gene_pair_dict[organism.name].append({
                    "gene_pair": (gene_pair[0].ID, gene_pair[1].ID),
                    "sequence_id": seq_id,
                    "sequence": self.sequence_dict[seq_id],
                    "overlap_type": overlap_type
                })

        return gene_pair_dict

    def write_fasta_for_intergenomic_sequences(self, output_dir: str):
        """
        Write the intergenomic sequences to separate FASTA files for each organism.

        :param output_dir: Directory where the FASTA files will be saved.
        """
        # Ensure output directory exists
        import os
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Iterate over the sequences for each organism
        for organism_name, sequences in self.extract_sequences().items():
            # Create a new FASTA file for this organism
            fasta_file_path = os.path.join(output_dir, f"{organism_name}_intergenomic_sequences.fasta")
            with open(fasta_file_path, 'w') as fasta_file:
                for seq_data in sequences:
                    gene_pair = seq_data['gene_pair']
                    seq_id = seq_data['sequence_id']
                    sequence = seq_data['sequence']
                    overlap_type = seq_data['overlap_type']

                    # Create a header for each gene pair
                    header = f">{organism_name}_gene_pair_{gene_pair[0]}_{gene_pair[1]}_seqid_{seq_id}_overlap_{overlap_type}"
                    fasta_file.write(f"{header}\n")
                    fasta_file.write(f"{sequence}\n")

            print(f"FASTA file for {organism_name} written to {fasta_file_path}")
