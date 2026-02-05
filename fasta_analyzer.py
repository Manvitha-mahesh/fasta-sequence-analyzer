def validate_sequence(sequence):
    """
    Validate that a sequence contains only valid nucleotide characters.
    
    Args:
        sequence (str): The nucleotide sequence to validate
        
    Returns:
        bool: True if sequence is valid
        
    Raises:
        ValueError: If sequence is empty or contains invalid characters
    """
    # Check if sequence is empty
    if not sequence or len(sequence.strip()) == 0:
        raise ValueError("Error: Sequence cannot be empty. Please check your FASTA file.")
    
    # Valid nucleotide characters (DNA/RNA)
    # A, T, G, C for DNA
    # A, U, G, C for RNA
    # N for unknown nucleotide
    valid_bases = set("ATGCUN")
    
    # Find invalid characters
    invalid_chars = set(sequence.upper()) - valid_bases
    
    if invalid_chars:
        raise ValueError(
            f"Error: Invalid characters found in sequence: {', '.join(sorted(invalid_chars))}. "
            f"Valid bases are: A, T, G, C, U, N (case insensitive)"
        )
    
    return True


def read_fasta(file_path):
    """
    Read a FASTA file and extract header and sequence with validation.
    
    Args:
        file_path (str): Path to the FASTA file
        
    Returns:
        tuple: (header, sequence) where both are strings
        
    Raises:
        FileNotFoundError: If file does not exist
        IOError: If there are issues reading the file
        ValueError: If sequence is empty or contains invalid characters
    """
    header = ""
    sequence = ""

    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                # Skip empty lines
                if not line:
                    continue
                # Process header
                if line.startswith(">"):
                    header = line
                # Process sequence
                else:
                    sequence += line.upper()
    
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_path}' not found. Please check the file path.")
    except IOError as e:
        raise IOError(f"Error reading file '{file_path}': {str(e)}")
    
    # Validate the sequence after reading
    if sequence:  # Only validate if sequence exists
        validate_sequence(sequence)
    else:
        raise ValueError("Error: No sequence found in the FASTA file. Please check the file format.")
    
    return header, sequence


def identify_sequence_type(sequence):
    """
    Identify whether a sequence is DNA, RNA, or Protein.
    
    Args:
        sequence (str): The nucleotide sequence (should be uppercase)
        
    Returns:
        str: "DNA", "RNA", or "Protein"
    """
    dna_letters = set("ATGC")
    rna_letters = set("AUGC")
    
    seq_set = set(sequence)
    
    # Check for RNA (contains U but no T)
    if "U" in seq_set and "T" not in seq_set:
        return "RNA"
    
    # Check for DNA (contains T but no U)
    elif "T" in seq_set and "U" not in seq_set:
        return "DNA"
    
    # If both T and U, or neither, check against standard bases
    elif seq_set.issubset(dna_letters):
        return "DNA"
    elif seq_set.issubset(rna_letters):
        return "RNA"
    else:
        return "Protein"


def gc_content(sequence):
    """
    Calculate GC content (percentage of G and C nucleotides).
    
    Args:
        sequence (str): The nucleotide sequence
        
    Returns:
        float: GC content as a percentage
        
    Raises:
        ValueError: If sequence is empty
    """
    if not sequence:
        raise ValueError("Error: Cannot calculate GC content for an empty sequence.")
    
    g = sequence.count("G")
    c = sequence.count("C")
    gc_percent = ((g + c) / len(sequence)) * 100
    
    return gc_percent


def at_content(sequence):
    """
    Calculate AT content (percentage of A and T nucleotides).
    
    Args:
        sequence (str): The nucleotide sequence
        
    Returns:
        float: AT content as a percentage
    """
    if not sequence:
        raise ValueError("Error: Cannot calculate AT content for an empty sequence.")
    
    a = sequence.count("A")
    t = sequence.count("T")
    at_percent = ((a + t) / len(sequence)) * 100
    
    return at_percent


def main():
    """
    Main function to analyze a FASTA file with error handling.
    """
    fasta_file = "data/sample.fasta"
    
    try:
        # Read and validate the FASTA file
        header, sequence = read_fasta(fasta_file)
        
        # Identify sequence type
        seq_type = identify_sequence_type(sequence)
        
        # Print results
        print("\n" + "="*50)
        print("FASTA Sequence Analysis")
        print("="*50)
        print("Header:", header)
        print("Sequence length:", len(sequence))
        print("Sequence type:", seq_type)
        
        # Calculate and display content percentages
        if seq_type in ["DNA", "RNA"]:
            gc = gc_content(sequence)
            print("GC content: {:.2f}%".format(gc))
            
            # Show AT content if DNA
            if seq_type == "DNA":
                at = at_content(sequence)
                print("AT content: {:.2f}%".format(at))
        
        print("="*50 + "\n")
    
    except FileNotFoundError as e:
        print(f"\n❌ {e}\n")
    except IOError as e:
        print(f"\n❌ {e}\n")
    except ValueError as e:
        print(f"\n❌ {e}\n")
    except Exception as e:
        print(f"\n❌ Unexpected error: {str(e)}\n")


if __name__ == "__main__":
    main()
