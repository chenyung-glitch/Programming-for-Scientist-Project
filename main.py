from genome import generate_random_genome, generate_random_short_reads, generate_kmers_for_read_list, generate_kmer

def main():
    print("greedy_genome_gobblers")
    genome = generate_random_genome(10)
    read_list = generate_random_short_reads(genome, 3, 5)
    kmer_list = generate_kmers_for_read_list(read_list, 2)
    print(genome)
    print(read_list)
    print(kmer_list)

if __name__ == "__main__":
    main()