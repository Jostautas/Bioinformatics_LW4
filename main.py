from Bio import SeqIO


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


if __name__ == '__main__':
    input_file = "./input.fastq"
    output_file = "./filtered_sequences.fastq"

    filter_sequences_by_length(input_file, output_file)
