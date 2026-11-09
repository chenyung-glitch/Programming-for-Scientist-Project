"""
Greedy Genome Gobblers - milestone demonstration script.

This script connects the pieces that already existed in the repository:

    genome.generate_kmers_for_read_list()  ->  list of k-mers
    de_brujin.DeBruijnGraph.add_kmer()     ->  graph nodes and edges
    de_brujin.eulerian_path()              ->  ordered node path
    de_brujin.path_to_sequence()           ->  reconstructed DNA

It does NOT add a new assembly algorithm and it does NOT evaluate accuracy on
realistic read sets. The goal is to show that the existing implementation can
reconstruct a small synthetic genome whose answer we already know.

Usage:
    python main.py             # deterministic ATGCGT example (primary demo)
    python main.py --random    # optional seeded random example, same pipeline
"""

import random
import sys

from data.small_dataset import SMALL_EXAMPLE
from de_brujin import DeBruijnGraph, eulerian_path, path_to_sequence
from genome import (
    generate_kmer,
    generate_kmers_for_read_list,
    generate_random_genome,
    generate_random_short_reads,
)


def assemble(reads, k):
    """
    Build a De Bruijn graph from reads and return one reconstruction.

    input:
    reads (list[str]): short read strings
    k (int): k-mer length

    output:
    dict: graph, k-mers, Eulerian node path (None if there is no valid path)
          and the reconstructed sequence (None if reconstruction failed)
    """

    kmer_list = generate_kmers_for_read_list(reads, k)

    graph = DeBruijnGraph()
    for kmer in kmer_list:
        graph.add_kmer(kmer)

    path = eulerian_path(graph)

    # Never invent a sequence when the graph does not admit an Eulerian path.
    sequence = path_to_sequence(path) if path is not None else None

    return {
        "kmers": kmer_list,
        "graph": graph,
        "path": path,
        "sequence": sequence,
    }


def report(reference_genome, reads, k, result):
    """
    Print one assembly attempt in a fixed, readable format.

    input:
    reference_genome (str): ground truth used by the dataset
    reads (list[str]): the reads that were assembled
    k (int): k-mer length
    result (dict): output of assemble()
    """

    graph = result["graph"]
    sequence = result["sequence"]

    print("Reference genome: {}".format(reference_genome))
    print("Reads: {}".format(reads))
    print("K-mer length: {}".format(k))
    print()

    print("Generated k-mers:")
    print(", ".join(result["kmers"]))
    print()

    print("Graph nodes: {}".format(len(graph.nodes)))
    print("Graph edges: {}".format(len(graph.edges)))
    print()

    if sequence is None:
        print("Reconstructed sequence: <none>")
        print("Reconstruction successful: False")
        print(
            "Reason: the graph does not admit a single Eulerian path over all "
            "edges, so no sequence was produced."
        )
        return False

    print("Reconstructed sequence: {}".format(sequence))
    successful = sequence == reference_genome
    print("Reconstruction successful: {}".format(successful))

    if not successful:
        print(
            "Note: the assembled sequence differs from the reference genome. "
            "A valid Eulerian path over the observed k-mers is not by itself "
            "proof of a correct biological assembly."
        )

    return successful


def run_small_example():
    """Run the deterministic ATGCGT example and return whether it matched."""

    dataset = SMALL_EXAMPLE
    k = dataset["k"]

    # Confirm the fixture's reads really produce the documented k-mers, so any
    # future edit to the data or to generate_kmer() fails loudly.
    expected_kmers = dataset["expected_kmers"]
    actual_kmers = generate_kmers_for_read_list(dataset["reads"], k)
    if actual_kmers != expected_kmers:
        print("Dataset error: k-mers {} do not match expected {}".format(
            actual_kmers, expected_kmers
        ))
        return False

    print("Greedy Genome Gobblers")
    print("Milestone: De Bruijn graph assembly of a known synthetic genome")
    print("Dataset: {}".format(dataset["name"]))
    print()

    result = assemble(dataset["reads"], k)
    return report(dataset["reference_genome"], dataset["reads"], k, result)


def run_random_example(seed=601, genome_length=40, read_length=8,
                       read_number=12, k=5):
    """
    Optional seeded random example.

    The seed makes the run reproducible: the same seed always produces the same
    genome, the same reads, and therefore the same graph. Random sampling of
    reads may leave parts of the genome uncovered, so a failed reconstruction
    here is expected behaviour, not a bug in the graph traversal.
    """

    random.seed(seed)

    reference_genome = generate_random_genome(genome_length)
    reads = generate_random_short_reads(reference_genome, read_length, read_number)

    print("Greedy Genome Gobblers")
    print("Optional random example (seed = {})".format(seed))
    print("Genome length: {}  read length: {}  reads: {}  k: {}".format(
        genome_length, read_length, read_number, k
    ))
    print()

    result = assemble(reads, k)
    report(reference_genome, reads, k, result)

    # Show the k-mers of the first read so the overlap structure stays visible.
    print()
    print("Example k-mers of the first read: {}".format(
        generate_kmer(reads[0], k)
    ))
    return result


def main():
    run_small_example()

    if "--random" in sys.argv[1:]:
        print()
        print("=" * 60)
        print()
        run_random_example()


if __name__ == "__main__":
    main()
