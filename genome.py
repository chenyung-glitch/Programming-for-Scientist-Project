from random import randint

def generate_random_genome(g_length = int) -> str:
    """
    takes as input an integer length and returns a randomly generated genome string with size of g_length
    
    input:
    g_length (int): length of genome
    
    output:
    str: genome string with length of g_length
    """

    genome_map = {1 : "A",2 : "T",3 : "C",4 : "G"}
    genome = ""
    for _ in range(g_length):
        genome += genome_map[randint(1,4)]
    return genome

def generate_random_short_reads(genome = str, r_length = int, r_number = int) -> list[str]:
    """
    takes as input a string and two ints as the genome string, length of the read, and number of reads and returns
    a list of size r_number with reads of size r_length
    
    input:
    genome (str): genome string
    r_length (int): length of the read
    r_number (int): the number of reads
    
    output:
    list[str]: list of size r_number with reads with length r_length
    """

    genome_upper_bound = len(genome) - r_length
    read_list = []
    for _ in range(r_number):
        read_position = randint(0,genome_upper_bound)
        read_list.append(genome[read_position : read_position + r_length])
    return read_list

def generate_kmers_for_read_list(read_list = list[str], k_length = int) -> list[str]:
    """
    takes as input a list of strings read_list and an int k_length and returns a list of all substrings for each read
    in read_list with k-mer length of k_length
    
    input:
    read_list (list[str]): list of read string
    k_length (int): length of the k-mer
    
    output:
    list[str]: list of all substrings for each read in read_list with k-mers of length k_length
    """

    kmer_list = []
    for read in read_list:
        kmer_list += generate_kmer(read, k_length)
    return kmer_list
    
def generate_kmer(read = str, k_length = int) -> list[str]:
    """
    takes as input a string read and an int k_length and returns a list of all substrings of read with k-mer length of k_length
    
    input:
    read (str): read string
    k_length (int): length of the k-mer
    
    output:
    list[str]: list of all substrings of string read with k-mers of length k_length
    """

    kmer_list = []
    for pos in range(0, len(read) - k_length + 1):
        kmer_list.append(read[pos : pos + k_length])
    return kmer_list