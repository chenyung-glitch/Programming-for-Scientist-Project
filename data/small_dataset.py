"""
Deterministic synthetic datasets for the Greedy Genome Gobblers assembler.

The project proposal asks us to study genome assembly on *simulated* DNA rather
than downloaded biological data, so the fixtures in this module are small,
hand-checked, and fully reproducible: every value below is a literal, nothing is
random, and no file or network access is involved.

Why keep the examples so small?
-------------------------------
A De Bruijn graph built from real short reads is only an approximation of the
original genome. Repeats, uneven coverage, and missing overlaps all create
branching and disconnected pieces. To test the *implementation* (nodes, edges,
Eulerian traversal, sequence reconstruction) we first need a case where the
correct answer is known by hand, and where the graph is guaranteed to hold a
single Eulerian path. Accuracy on imperfect read sets is a later milestone.

Dataset 1 - SMALL_EXAMPLE (primary, used by `python main.py`)
-------------------------------------------------------------
    reference genome : "ATGCGT"
    reads            : ["ATGC", "GCGT"]
    k                : 3

Each read is broken into 3-mers:

    "ATGC" -> "ATG", "TGC"
    "GCGT" -> "GCG", "CGT"

so the k-mer list is ["ATG", "TGC", "GCG", "CGT"].

Every 3-mer becomes a directed edge between two 2-mers (k-1 = 2):

    AT -ATG-> TG -TGC-> GC -GCG-> CG -CGT-> GT

That chain is a single Eulerian path, so reading the node labels in order and
appending the last base of each new node gives back "ATGCGT".

Dataset 2 - LARGER_EXAMPLE (secondary, exercised by the tests)
-------------------------------------------------------------
A slightly longer, still deterministic example. It is built the same way and
is included to show the pipeline is not special-cased to six bases.
"""

# Smallest interesting case: four k-mers, five nodes, one Eulerian path.
SMALL_EXAMPLE = {
    "name": "small_6bp_2reads",
    "reference_genome": "ATGCGT",
    "reads": ["ATGC", "GCGT"],
    "k": 3,
    "expected_kmers": ["ATG", "TGC", "GCG", "CGT"],
    "expected_assembly": "ATGCGT",
    "note": (
        "Hand-checkable fixture. The four k-mers form one Eulerian path "
        "AT -> TG -> GC -> CG -> GT."
    ),
}

# A longer deterministic case with overlapping reads, still a single path.
LARGER_EXAMPLE = {
    "name": "larger_9bp_2reads",
    "reference_genome": "ACGTTGCA",
    "reads": ["ACGTT", "GTTGCA"],
    "k": 4,
    "expected_assembly": "ACGTTGCA",
    "note": (
        "Two overlapping reads sharing the 3-mer GTT, so the graph stays "
        "connected and Eulerian."
    ),
}

# Convenience lookup so callers can select a fixture by name.
DATASETS = {
    SMALL_EXAMPLE["name"]: SMALL_EXAMPLE,
    LARGER_EXAMPLE["name"]: LARGER_EXAMPLE,
}


def get_dataset(name):
    """
    Return one deterministic fixture by name.

    input:
    name (str): key of the requested dataset, see DATASETS

    output:
    dict: the dataset dictionary (reference genome, reads, k, expected assembly)
    """

    return DATASETS[name]
