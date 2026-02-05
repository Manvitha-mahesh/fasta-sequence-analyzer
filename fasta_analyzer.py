def read_fasta(file_path):
    header = ""
    sequence = ""

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                header = line
            else:
                sequence += line.upper()

    return header, sequence


def identify_sequence_type(sequence):
    dna_letters = set("ATGC")
    if set(sequence).issubset(dna_letters):
        return "DNA"
    else:
        return "Protein"


def gc_content(sequence):
    g = sequence.count("G")
    c = sequence.count("C")
    return ((g + c) / len(sequence)) * 100


def main():
    fasta_file = "data/sample.fasta"
    header, sequence = read_fasta(fasta_file)

    seq_type = identify_sequence_type(sequence)

    print("Header:", header)
    print("Sequence length:", len(sequence))
    print("Sequence type:", seq_type)

    if seq_type == "DNA":
        gc = gc_content(sequence)
        print("GC content: {:.2f}%".format(gc))


if __name__ == "__main__":
    main()
