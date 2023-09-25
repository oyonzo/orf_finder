
START = "ATG"

def is_stop(codon):
    """
    Determines if a string is a stop codon
    
    Args:
        codon: A string representing the codon to be tested
    
    Returns:
        A boolean signifying if codon is or is not a stop codon.
    """
    codon = codon.upper()
    if codon == "TGA" or codon == "TAG" or codon == "TAA":
        return True
    else:
        return False

def orf_sequence(dna_sequence):
    """
    Finds and returns the first orf sequence in a dna seqeunce
    
    Args:
        dna_sequence: A string representing the dna sequence to be sequenced
        
    Return:
        A string representing the orf
    """
    orf = ""
    dna_sequence = dna_sequence.upper()
    for i in range(0,len(dna_sequence)):
        if dna_sequence[i:i+3] == START:
            while is_stop(dna_sequence[i:i+3]) == False:
                orf += dna_sequence[i:i+3]
                i += 3
                if is_stop(dna_sequence[i:i+3]) == True or i > len(dna_sequence):
                    return orf
        return orf
            
def find_orfs(dna_sequence):
    """
    Finds returns all the orf sequences in a dna sequence
    
    Args:
        dna_sequence: A string representing the dna sequence to be sequenced
        
    Return:
        A list of strings representing all of the orfs
    """
    orfs = []
    dna_sequence = dna_sequence.upper()
    i = 0
    while i <= len(dna_sequence):
        if dna_sequence[i:i+3] == START:
            orfs.append(orf_sequence(dna_sequence[i:]))
            i += len(orf_sequence(dna_sequence[i:]))
        i += 3
    return orfs

def reverse_complement(dna_sequence):
    """
    Reverses a DNA sequence and replaces each base with its complement
    
    Args:
        dna_sequence: A string representing the DNA sequence to reversed and complemented
        
    Return:
        A string representing the reverse complement of the given DNA sequence
    """
    dna_sequence = dna_sequence.upper()
    reverse = ""
    for i in range(len(dna_sequence)):
        reverse += dna_sequence[-(i+1)]
    reverse = reverse.replace("A", "t")
    reverse = reverse.replace("T", "a")
    reverse = reverse.replace("C", "g")
    reverse = reverse.replace("G", "c")
    
    return reverse.upper()

def gene_finder(dna_sequence):
    """
    Finds all the genes for a given dna_sequence for all three reading frames on both the given strand and its reverse complement
    
    Args:
        dna_sequence: A string that represents the DNA sequence in which genes will looked for
        
    Return:
        A list of strings representing all the genes found
    """
    genes = []
    genes += find_orfs(dna_sequence)
    genes += find_orfs(dna_sequence[1:])
    genes += find_orfs(dna_sequence[2:])
    genes += find_orfs(reverse_complement(dna_sequence))
    genes += find_orfs(reverse_complement(dna_sequence)[1:])
    genes += find_orfs(reverse_complement(dna_sequence)[2:])
    
    return genes
    
def read_fasta(filename):
    """
    Read a single DNA sequence from a FASTA-formatted file
    
    For example, to read the sequence from a file named "X73525.fasta.txt"
    >>> sequence = read_fasta("X73525.fasta.txt")
    
    Args:
        filename: Filename as a string
        
    Returns: Upper case DNA sequence as a string
    """
    with open(filename, "r") as file:
        # Read (and ignore) header line
        header = file.readline()
        # Read sequence lines
        sequence = ""
        for line in file:
            sequence += line.strip()
        return sequence.upper()

def filter_orfs(orfs, min_length):
    """
    Filter ORFs to have a minimum length
    
    Args:
        orfs: List of candidate ORF sequences
        min_length: Minimum length for an ORF
    
    Returns:
        A new list containing only ORF strings longer than min_length bases
    """
    filtered_orfs = []
    for orf in orfs:
        if len(orf) > min_length:
            filtered_orfs.append(orf)
    return filtered_orfs


def write_fasta(filename, orfs):
    """
    Write list of ORFs to a FASTA formatted text file.
    
    For example, to save a list of orfs assigned to the variable my_orfs to a
    file named "genes.txt"
    >>> write_fasta("genes.txt", my_orfs)
    
    Args:
        filename: Filename as a string. Note that any existing file with this name
            will be overwritten.
        orfs: List of ORF sequences to write to the file
    """
    with open(filename, "w") as file:
        for i in range(len(orfs)):
            # A FASTA entry is a header line that begins with a '>', and then the sequence on the next line(s)
            print(">seq" + str(i), file=file)
            print(orfs[i], file=file)

           
if __name__ == "__main__":
    # If you want perform the steps for creating genes.txt as part of your
    # program, implement that code here. It will run every time you click the
    # green arrow, but will be ignored by Gradescope.
    pass
