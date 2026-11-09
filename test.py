"""
Unit tests for the Greedy Genome Gobblers milestone.

Run with:
    python -m unittest test -v

The tests use only the standard library (unittest + random) so that no extra
dependency is needed beyond requirements.txt.
"""

import random
import os
import tempfile
import unittest

import matplotlib.image as mpimg

from data.small_dataset import LARGER_EXAMPLE, SMALL_EXAMPLE
from de_brujin import DeBruijnGraph, Edge, eulerian_path, path_to_sequence
from genome import (
    generate_kmer,
    generate_kmers_for_read_list,
    generate_random_genome,
    generate_random_short_reads,
)
import graph_visualization as gv

DNA_BASES = set("ATCG")


def build_graph(kmers):
    """Small helper: build a DeBruijnGraph from a list of k-mers."""

    graph = DeBruijnGraph()
    for kmer in kmers:
        graph.add_kmer(kmer)
    return graph


class TestGenomeGeneration(unittest.TestCase):
    """genome.generate_random_genome"""

    def test_length_and_alphabet(self):
        random.seed(1234)
        genome = generate_random_genome(50)
        self.assertEqual(len(genome), 50)
        self.assertTrue(set(genome).issubset(DNA_BASES))

    def test_seed_is_reproducible(self):
        random.seed(7)
        first = generate_random_genome(20)
        random.seed(7)
        second = generate_random_genome(20)
        self.assertEqual(first, second)


class TestReadGeneration(unittest.TestCase):
    """genome.generate_random_short_reads"""

    def test_read_count_and_length(self):
        random.seed(42)
        genome = generate_random_genome(30)
        reads = generate_random_short_reads(genome, 6, 8)

        self.assertEqual(len(reads), 8)
        for read in reads:
            self.assertEqual(len(read), 6)
            # Every read must be a substring of the genome it came from.
            self.assertIn(read, genome)

    def test_seed_is_reproducible(self):
        random.seed(11)
        genome = generate_random_genome(25)
        first = generate_random_short_reads(genome, 5, 6)

        random.seed(11)
        genome_again = generate_random_genome(25)
        second = generate_random_short_reads(genome_again, 5, 6)

        self.assertEqual(genome, genome_again)
        self.assertEqual(first, second)


class TestKmerGeneration(unittest.TestCase):
    """genome.generate_kmer / generate_kmers_for_read_list"""

    def test_single_read_overlapping_substrings(self):
        # A read of length 4 contains 4 - 3 + 1 = 2 overlapping 3-mers.
        self.assertEqual(generate_kmer("ATGC", 3), ["ATG", "TGC"])

    def test_read_list_is_concatenation_of_reads(self):
        reads = ["ATGC", "GCGT"]
        self.assertEqual(
            generate_kmers_for_read_list(reads, 3),
            ["ATG", "TGC", "GCG", "CGT"],
        )

    def test_kmer_longer_than_read_gives_no_kmers(self):
        # A read shorter than k has no valid k-mer window; the function should
        # return an empty list instead of a malformed fragment.
        self.assertEqual(generate_kmer("AT", 3), [])

    def test_kmer_length_equals_read_length(self):
        self.assertEqual(generate_kmer("ATGC", 4), ["ATGC"])

    def test_small_example_kmers_match_fixture(self):
        dataset = SMALL_EXAMPLE
        self.assertEqual(
            generate_kmers_for_read_list(dataset["reads"], dataset["k"]),
            dataset["expected_kmers"],
        )


class TestGraphConstruction(unittest.TestCase):
    """de_brujin.DeBruijnGraph.add_kmer"""

    def test_nodes_are_k_minus_1_mers(self):
        graph = build_graph(SMALL_EXAMPLE["expected_kmers"])
        self.assertEqual(
            set(graph.nodes.keys()), {"AT", "TG", "GC", "CG", "GT"}
        )
        for node_key, node in graph.nodes.items():
            self.assertEqual(node.node_key, node_key)

    def test_edges_match_kmers(self):
        graph = build_graph(SMALL_EXAMPLE["expected_kmers"])
        self.assertEqual(
            [edge.kmer for edge in graph.edges], SMALL_EXAMPLE["expected_kmers"]
        )

    def test_adjacency_links_source_and_target(self):
        graph = build_graph(["ATG"])
        at_node = graph.nodes["AT"]
        tg_node = graph.nodes["TG"]

        self.assertEqual(len(at_node.out_edges), 1)
        self.assertEqual(len(at_node.in_edges), 0)
        self.assertEqual(len(tg_node.in_edges), 1)
        self.assertEqual(len(tg_node.out_edges), 0)

        edge = at_node.out_edges[0]
        self.assertIs(edge.source, at_node)
        self.assertIs(edge.target, tg_node)
        self.assertIs(tg_node.in_edges[0], edge)

    def test_shared_node_connects_two_reads(self):
        graph = build_graph(["ATG", "TGC"])
        self.assertEqual(len(graph.nodes["TG"].in_edges), 1)
        self.assertEqual(len(graph.nodes["TG"].out_edges), 1)


class TestReconstruction(unittest.TestCase):
    """eulerian_path + path_to_sequence on the deterministic fixtures"""

    def test_small_example_reconstructs_atgcgt(self):
        graph = build_graph(SMALL_EXAMPLE["expected_kmers"])
        path = eulerian_path(graph)

        self.assertIsNotNone(path)
        self.assertEqual(len(path), len(graph.edges) + 1)  # all edges used
        self.assertEqual(path_to_sequence(path), "ATGCGT")

    def test_small_example_via_main_pipeline(self):
        # Same calls main.py uses, to make sure the integrated path works.
        reads = SMALL_EXAMPLE["reads"]
        kmers = generate_kmers_for_read_list(reads, SMALL_EXAMPLE["k"])
        graph = build_graph(kmers)
        sequence = path_to_sequence(eulerian_path(graph))

        self.assertEqual(sequence, SMALL_EXAMPLE["expected_assembly"])

    def test_larger_example_reconstructs_reference(self):
        dataset = LARGER_EXAMPLE
        kmers = generate_kmers_for_read_list(dataset["reads"], dataset["k"])
        graph = build_graph(kmers)
        sequence = path_to_sequence(eulerian_path(graph))

        self.assertEqual(sequence, dataset["expected_assembly"])

    def test_path_to_sequence_empty_path(self):
        self.assertEqual(path_to_sequence([]), "")

    def test_cycle_starts_anywhere_and_still_reconstructs(self):
        # "ATAT" with k = 3 produces ATA, TAT, ATA: a graph with one unused
        # Eulerian trail, used here only to confirm the cycle fallback works.
        graph = build_graph(["ATA", "TAT", "ATA"])
        path = eulerian_path(graph)

        self.assertIsNotNone(path)
        self.assertEqual(len(path), len(graph.edges) + 1)
        # All three observed k-mers are used, so the repeated ATA appears twice.
        self.assertEqual(path_to_sequence(path), "ATATA")


class TestInvalidConfigurations(unittest.TestCase):
    """The traversal must refuse graphs without a valid Eulerian path."""

    def test_disconnected_graph_returns_none(self):
        # ATG (AT -> TG) and CGT (CG -> GT) share no node.
        graph = build_graph(["ATG", "CGT"])
        self.assertIsNone(eulerian_path(graph))

    def test_unbalanced_degrees_return_none(self):
        # Two parallel ATG edges: AT has out-degree 2, TG in-degree 2, which
        # violates the +-1 rule, so no Eulerian path exists.
        graph = build_graph(["ATG", "ATG"])
        self.assertIsNone(graph.find_start_node())
        self.assertIsNone(eulerian_path(graph))

    def test_multiple_start_nodes_return_none(self):
        # Two separate start nodes, each with out - in == 1.
        graph = build_graph(["ATG", "CGA"])
        degree_one_starts = [
            key for key, node in graph.nodes.items()
            if len(node.out_edges) - len(node.in_edges) == 1
        ]
        self.assertEqual(len(degree_one_starts), 2)
        self.assertIsNone(graph.find_start_node())
        self.assertIsNone(eulerian_path(graph))

    def test_empty_graph_returns_none(self):
        graph = DeBruijnGraph()
        self.assertIsNone(graph.find_start_node())
        self.assertIsNone(eulerian_path(graph))

    def test_legal_degrees_but_two_components_return_none(self):
        # A balanced 4-edge cycle plus one separate edge: every degree is legal
        # and there is only one start (CG) and one end (GA), yet the edges live
        # in two components, so the traversal cannot consume all of them.
        graph = build_graph(["ATG", "TGC", "GCA", "CAT", "CGA"])
        self.assertIsNotNone(graph.find_start_node())
        self.assertIsNone(eulerian_path(graph))

    def test_no_fabricated_sequence_when_traversal_fails(self):
        graph = build_graph(["ATG", "CGT"])
        path = eulerian_path(graph)
        sequence = path_to_sequence(path) if path is not None else None
        self.assertIsNone(sequence)


class TestDuplicateKmers(unittest.TestCase):
    """Repeated k-mers must stay separate edge objects."""

    def test_duplicate_kmers_create_separate_edges(self):
        graph = build_graph(["ATG", "ATG"])

        self.assertEqual(len(graph.edges), 2)
        self.assertIsNot(graph.edges[0], graph.edges[1])

        at_node = graph.nodes["AT"]
        self.assertEqual(len(at_node.out_edges), 2)
        self.assertIsNot(at_node.out_edges[0], at_node.out_edges[1])
        self.assertIsInstance(at_node.out_edges[0], Edge)

    def test_duplicate_kmers_are_not_collapsed(self):
        graph = build_graph(["ATA", "TAT", "ATA"])
        self.assertEqual(len(graph.edges), 3)
        self.assertEqual(len(graph.nodes["AT"].out_edges), 2)

    def test_duplicate_kmer_edge_degrees_are_counted(self):
        graph = build_graph(["ATG", "ATG"])
        # Degrees are counted per edge object, not per distinct k-mer string.
        self.assertEqual(len(graph.nodes["AT"].out_edges), 2)
        self.assertEqual(len(graph.nodes["TG"].in_edges), 2)


class TestGraphVisualization(unittest.TestCase):
    """graph_visualization reads the graph, it never changes it."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.small_graph = build_graph(SMALL_EXAMPLE["expected_kmers"])
        self.small_path = eulerian_path(self.small_graph)

    def tearDown(self):
        self.temp_dir.cleanup()

    def output(self, name):
        return os.path.join(self.temp_dir.name, name)

    # --- element extraction ---------------------------------------------

    def test_node_labels_match_graph_nodes(self):
        elements = gv.graph_elements(self.small_graph)
        labels = sorted(node["label"] for node in elements["nodes"])
        self.assertEqual(labels, sorted(self.small_graph.nodes.keys()))

    def test_edge_labels_match_graph_kmers(self):
        elements = gv.graph_elements(self.small_graph)
        labels = [edge["kmer"] for edge in elements["edges"]]
        self.assertEqual(labels, [edge.kmer for edge in self.small_graph.edges])

    def test_start_and_end_nodes_are_highlighted(self):
        elements = gv.graph_elements(self.small_graph)
        self.assertEqual(elements["start_keys"], ["AT"])
        self.assertEqual(elements["end_keys"], ["GT"])

    def test_layout_is_deterministic(self):
        first = gv.graph_elements(self.small_graph)
        second = gv.graph_elements(self.small_graph)
        self.assertEqual(
            [node["position"] for node in first["nodes"]],
            [node["position"] for node in second["nodes"]],
        )

    # --- path to edge mapping -------------------------------------------

    def test_ordered_edges_follow_the_node_path(self):
        from graph_visualization import ordered_edges_for_path

        ordered = ordered_edges_for_path(self.small_graph, self.small_path)
        self.assertIsNotNone(ordered)
        self.assertEqual(len(ordered), len(self.small_graph.edges))
        self.assertEqual(
            [edge.kmer for edge in ordered],
            SMALL_EXAMPLE["expected_kmers"],
        )

    def test_duplicate_edges_are_not_reused(self):
        from graph_visualization import ordered_edges_for_path

        graph = build_graph(["ATA", "TAT", "ATA"])
        path = eulerian_path(graph)
        ordered = ordered_edges_for_path(graph, path)

        self.assertIsNotNone(ordered)
        self.assertEqual(len(ordered), 3)
        self.assertEqual(len({id(edge) for edge in ordered}), 3)

    def test_ordered_edges_reject_unknown_path(self):
        from graph_visualization import ordered_edges_for_path

        # AT -> TG exists, but TG -> GG does not, so the step cannot be resolved.
        graph = build_graph(["ATG"])
        fake_path = [graph.nodes["TG"], graph.nodes["AA"] if "AA" in graph.nodes
                     else graph.nodes["AT"]]
        self.assertIsNone(ordered_edges_for_path(graph, fake_path))
        self.assertEqual(ordered_edges_for_path(graph, None), [])

    # --- graph must not be mutated --------------------------------------

    def test_rendering_does_not_mutate_the_graph(self):
        before = (
            len(self.small_graph.nodes),
            len(self.small_graph.edges),
            {key: [id(edge) for edge in node.out_edges]
             for key, node in self.small_graph.nodes.items()},
        )

        gv.visualize_debruijn_graph(self.small_graph, self.output("g.png"))
        gv.visualize_eulerian_path(
            self.small_graph, self.small_path, self.output("t.png")
        )

        after = (
            len(self.small_graph.nodes),
            len(self.small_graph.edges),
            {key: [id(edge) for edge in node.out_edges]
             for key, node in self.small_graph.nodes.items()},
        )
        self.assertEqual(before, after)

    # --- rendering success ----------------------------------------------

    def test_graph_figure_is_written_and_nonempty(self):
        path = gv.visualize_debruijn_graph(self.small_graph, self.output("g.png"))
        self.assertEqual(path, self.output("g.png"))
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 1000)

    def test_eulerian_figure_is_written_and_nonempty(self):
        path = gv.visualize_eulerian_path(
            self.small_graph, self.small_path, self.output("t.png")
        )
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 1000)

    def test_sequence_figure_is_written_and_nonempty(self):
        path = gv.visualize_sequence_reconstruction(
            SMALL_EXAMPLE["reference_genome"],
            SMALL_EXAMPLE["reads"],
            SMALL_EXAMPLE["expected_kmers"],
            "ATGCGT",
            k=SMALL_EXAMPLE["k"],
            output_path=self.output("s.png"),
        )
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 1000)

    def test_missing_path_is_handled_safely(self):
        # A disconnected graph yields no Eulerian path; the figure must still be
        # produced and must not raise.
        graph = build_graph(["ATG", "CGT"])
        path = gv.visualize_eulerian_path(
            graph, eulerian_path(graph), self.output("broken.png")
        )
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 1000)

    def test_duplicate_edges_and_self_loops_render(self):
        graph = build_graph(["ATG", "ATG", "AAA", "TGA"])
        path = gv.visualize_debruijn_graph(graph, self.output("dup.png"))
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 1000)

        # The self-loop k-mer AAA maps to the node AA twice.
        self.assertEqual(graph.nodes["AA"].node_key, "AA")
        self.assertEqual(len(graph.nodes["AA"].out_edges), 1)

    def test_all_figures_are_created_together(self):
        written = gv.create_all_figures(
            graph=self.small_graph,
            path=self.small_path,
            reference_genome=SMALL_EXAMPLE["reference_genome"],
            reads=SMALL_EXAMPLE["reads"],
            kmers=SMALL_EXAMPLE["expected_kmers"],
            assembled_sequence="ATGCGT",
            output_directory=os.path.join(self.temp_dir.name, "figures"),
            k=SMALL_EXAMPLE["k"],
        )
        self.assertEqual(len(written), 3)
        for path in written:
            self.assertTrue(os.path.exists(path))
            self.assertGreater(os.path.getsize(path), 1000)

    def test_rendering_is_reproducible(self):
        first = gv.visualize_debruijn_graph(
            self.small_graph, self.output("first.png")
        )
        second = gv.visualize_debruijn_graph(
            self.small_graph, self.output("second.png")
        )

        first_pixels = mpimg.imread(first)
        second_pixels = mpimg.imread(second)
        self.assertEqual(first_pixels.shape, second_pixels.shape)
        self.assertTrue((first_pixels == second_pixels).all())


if __name__ == "__main__":
    unittest.main()
