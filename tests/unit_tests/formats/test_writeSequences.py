from unittest.mock import patch

from ppanggolin.formats.writeSequences import read_genome_file, write_regions_sequences
from ppanggolin.genome import Contig, Gene, Organism
from ppanggolin.pangenome import Pangenome
from ppanggolin.region import Region


def _parse_fasta(fasta_path):
    records = {}
    header = None
    chunks = []
    for line in fasta_path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header = line[1:].strip()
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def test_write_regions_sequences_one_based_and_circular(tmp_path):
    # 36-nt circular contig: 1..6 = ATGCCC, 10..21 = GGGTTTAAACCC, 31..36 = CCCTAA
    seq = "ATGCCC" + "AAA" + "GGGTTTAAACCC" + "TTTTTTTTT" + "CCCTAA"
    assert len(seq) == 36

    genome_fasta = tmp_path / "org1.fna"
    genome_fasta.write_text(f">ctg1\n{seq}\n")
    fasta_list = tmp_path / "genomes.tsv"
    fasta_list.write_text(f"org1\t{genome_fasta}\n")

    pan = Pangenome()
    org = Organism("org1")
    pan.add_organism(org)
    ctg = Contig(1, "ctg1", is_circular=True)
    ctg.length = len(seq)
    org.add(ctg)

    g0 = Gene("g0")
    g0.fill_annotations(start=1, stop=6, strand="+", position=0)
    g0.fill_parents(org, ctg)
    ctg.add(g0)

    g1 = Gene("g1")
    g1.fill_annotations(start=10, stop=21, strand="+", position=1)
    g1.fill_parents(org, ctg)
    ctg.add(g1)

    g2 = Gene("g2")
    g2.fill_annotations(start=31, stop=36, strand="+", position=2)
    g2.fill_parents(org, ctg)
    ctg.add(g2)

    r_start = Region("r_start")
    r_start.add(g0)
    r_start.add(g1)
    pan.add_region(r_start)

    r_mid = Region("r_mid")
    r_mid.add(g1)
    pan.add_region(r_mid)

    r_circ = Region("r_circ")
    r_circ.add(g2)
    r_circ.add(g0)
    pan.add_region(r_circ)

    with patch(
        "ppanggolin.formats.writeSequences.read_genome_file",
        wraps=read_genome_file,
    ) as mock_read:
        write_regions_sequences(
            pan, tmp_path, regions="all", fasta=fasta_list, disable_bar=True
        )
        assert mock_read.call_count == 1

    records = _parse_fasta(tmp_path / "all_rgp_genomic_sequences.fasta")
    assert records["r_start"] == seq[0:21]
    assert len(records["r_start"]) == r_start.length
    assert records["r_mid"] == seq[9:21]
    assert len(records["r_mid"]) == r_mid.length
    assert records["r_circ"] == seq[30:36] + seq[0:6]
    assert len(records["r_circ"]) == r_circ.length
