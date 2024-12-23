from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import subprocess
import pandas as pd


def filter_sequences_by_length(input_file, output_file, max_lines=None):
    sequences = []
    with open(input_file, "r") as file:
        current_line_count = 0
        while True:
            header = file.readline().strip()
            sequence = file.readline().strip()
            plus = file.readline().strip()
            quality = file.readline().strip()
            current_line_count += 4

            if not header or (max_lines and current_line_count > max_lines):
                break

            sequences.append(SeqRecord(Seq(sequence), id=header[1:], description="",
                                       letter_annotations={"phred_quality": [ord(c) - 33 for c in quality]}))

    sequence_lengths = [len(seq.seq) for seq in sequences]

    max_length = max(sequence_lengths)
    threshold_length = int(0.8 * max_length)

    filtered_sequences = [seq for seq in sequences if len(seq.seq) >= threshold_length]

    SeqIO.write(filtered_sequences, output_file, "fastq")

    print(f"Total sequences: {len(sequences)}")
    print(f"Sequences retained: {len(filtered_sequences)}")


def convert_fastq_to_fasta(input_file, output_file):
    SeqIO.convert(input_file, "fastq", output_file, "fasta")


def run_blast(input_file, output_file):
    try:
        subprocess.run(["C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\makeblastdb.exe", "-in", input_file, "-dbtype", "nucl"], check=True)

        result = subprocess.run(
            ["C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\blastn.exe", "-query", input_file, "-db", input_file, "-out", output_file, "-outfmt", "6"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print(f"run_blast() complete. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"BLAST failed: {e.stderr}")


def process_blast_output(blast_file, probe_min_len=35, probe_max_len=40, max_mismatches=2):
    columns = ["query_id", "subject_id", "%identity", "alignment_length", "mismatches",
               "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    blast_df = pd.read_csv(blast_file, sep="\t", names=columns)
    initial_sequences = blast_df["query_id"].nunique()
    print(f"Initial number of sequences in the dataset: {initial_sequences}")

    high_identity_alignments = blast_df[blast_df["%identity"] >= 98]

    probe_coverage = defaultdict(list)
    sequence_dict = {row["subject_id"]: fetch_sequence("./filtered_sequences.fasta", row["subject_id"])
                     for _, row in high_identity_alignments.iterrows()}

    for query_id, group in high_identity_alignments.groupby("query_id"):
        sequence_probes = []

        for _, row in group.iterrows():
            for start in range(row["q_start"], row["q_end"] - probe_min_len + 1):
                for probe_len in range(probe_min_len, probe_max_len + 1):
                    if start + probe_len > row["q_end"]:
                        continue
                    sequence_probes.append((query_id, start, start + probe_len))

        for probe in sequence_probes:
            matches = []
            probe_start, probe_end = probe[1], probe[2]
            for _, row in group.iterrows():
                target_start, target_end = row["s_start"], row["s_end"]
                target_seq = sequence_dict.get(row["subject_id"])
                if target_seq is None:
                    continue

                if target_start <= probe_start and target_end >= probe_end:
                    probe_seq = target_seq[probe_start - target_start:probe_end - target_start]
                    mismatches = sum(1 for a, b in zip(probe_seq, target_seq[probe_start - target_start:probe_end - target_start]) if a != b)
                    if mismatches <= max_mismatches:
                        matches.append(row["subject_id"])

            probe_coverage[probe] = matches

    probe_targets = {}
    selected_probes = []
    covered_sequences = set()
    while len(covered_sequences) < len(high_identity_alignments["subject_id"].unique()):
        best_probe = max(probe_coverage, key=lambda p: len(set(probe_coverage[p]) - covered_sequences))
        selected_probes.append(best_probe)
        probe_targets[best_probe] = list(set(probe_coverage[best_probe]))
        covered_sequences.update(probe_coverage[best_probe])
        del probe_coverage[best_probe]

    return selected_probes, probe_targets


def fetch_sequence(db, sequence_id):
    try:
        result = subprocess.run(
            ["C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\blastdbcmd.exe", "-db", db, "-entry", sequence_id],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        sequence = "".join(result.stdout.splitlines()[1:])  # Skip the FASTA header
        return sequence
    except subprocess.CalledProcessError as e:
        print(f"Error fetching sequence for {sequence_id}: {e.stderr}")
        return None


if __name__ == '__main__':
    input_file = "./input.fastq"
    filtered_fastq_file = "./filtered_sequences.fastq"
    filtered_fasta_file = "./filtered_sequences.fasta"
    aligned_file = "./aligned_sequences.fasta"
    aligned_file_small = "./aligned_sequences_small.fasta"

    # filter_sequences_by_length(input_file, filtered_fastq_file, 100)
    # convert_fastq_to_fasta(filtered_fastq_file, filtered_fasta_file)
    # run_blast(filtered_fasta_file, aligned_file)

    selected_probes, probe_targets = process_blast_output(aligned_file)
    print(f"Total probes selected: {len(selected_probes)}")
    for probe, targets in probe_targets.items():
        print(f"Probe: {probe}\nMatches sequences: {targets}")
