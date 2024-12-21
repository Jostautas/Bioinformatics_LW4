from Bio import SeqIO
import subprocess

def filter_sequences_by_length(input_file, output_file):
    # Parse the sequences from the input FASTQ file
    sequences = list(SeqIO.parse(input_file, "fastq"))

    # Calculate the length of each sequence
    sequence_lengths = [len(seq.seq) for seq in sequences]

    # Determine the threshold length (80% of the longest sequence)
    max_length = max(sequence_lengths)
    threshold_length = int(0.8 * max_length)

    # Filter sequences
    filtered_sequences = [seq for seq in sequences if len(seq.seq) >= threshold_length]

    # Write the filtered sequences to the output FASTQ file
    SeqIO.write(filtered_sequences, output_file, "fastq")

    print(f"Total sequences: {len(sequences)}")
    print(f"Sequences retained: {len(filtered_sequences)}")


def convert_fastq_to_fasta(input_file, output_file):
    SeqIO.convert(input_file, "fastq", output_file, "fasta")


def run_blast(input_file, output_file):
    """
    Performs sequence alignment using BLAST.

    Parameters:
    - input_file (str): Path to the input FASTA file.
    - output_file (str): Path to the BLAST results file.
    """
    try:
        # Create a BLAST database
        subprocess.run(["C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\makeblastdb.exe", "-in", input_file, "-dbtype", "nucl"], check=True)

        # Run BLAST
        result = subprocess.run(
            ["C:\\Program Files\\NCBI\\blast-2.16.0+\\bin\\blastn.exe", "-query", input_file, "-db", input_file, "-out", output_file, "-outfmt", "6"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print(f"BLAST alignment complete. Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"BLAST failed: {e.stderr}")


if __name__ == '__main__':
    input_file = "./input.fastq"
    filtered_fastq_file = "./filtered_sequences.fastq"
    filtered_fasta_file = "./filtered_sequences.fasta"
    aligned_file = "./aligned_sequences.fasta"

    filter_sequences_by_length(input_file, filtered_fastq_file)
    convert_fastq_to_fasta(filtered_fastq_file, filtered_fasta_file)
    run_blast(filtered_fasta_file, aligned_file)
