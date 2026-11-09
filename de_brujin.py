"""
De Brujin Graph Structural Def:
- Edges = K-mers
- Nodes = (K-1)-mers

Each k-mer gives us two (k-1)-mers
EX: ATG --> AT, TG
We can think of these as:
	Prefix = everything minus the last letter
	Suffix = everything minus the first letter

So every k-mer is immediately reconstructed using a graph
EX: AT --ATG--> TG --TGC--> GC --GCA--> CA --CAT--> AT

NOTE: Each neighboring k-mer overlaps by k-1characters

Neither path can use an edge more than once!
* Eulerian path ==> visits every edge in a graph exactly once
	- CAN pass through nodes more than once to do this
* Hamiltonian path ==> visits every node in a graph exactly once
	- DOESN't have to use every edge in the graph to do this

Because of De Brujin's structure where edges represent K-mers, we want to use each k-mer only once! 
Therefore, we use a Eulerian path traversal

Note how each read string broken up into k-mers starts with in-degree == 0 and out-degree == 0

Above, we assessed the k-mer and (k-1)-mer of a single read sequence and how it is represented in the De Brujin Graph

Now, we must consider how to incorporate multiple reads.
Identical (k-1)-mers automatically become the same node, which is what connects the reads.
The de Bruijn graph captures overlap through shared nodes

Sometimes, when doing short read sequencing, you don't get any overlap (happens with low coverage). No overlap can result in different connected components within the deBrujin Graph.

NOTE: condensing overlapped structures into shared nodes is very helpful for creating links between reads. However a few problems arise by doing this:

- repeat-induced ambiguity ==> repeated sequence can cause multiple possible paths through the graph
- local connectivity can be ambiguous, representing either:
	- multiple reads covering the same region
	- a repeated sequence
	- sequencing variation
	- or sequencing errors
    """


class Node:
    def __init__(self, node_key):
        self.node_key = node_key # each node is the (k-1)-mer
        self.out_edges = [] # store in + out edges to create adjacency list
        self.in_edges = []


class Edge:
    def __init__(self, kmer, source, target):
        self.kmer = kmer # NOTE: edges represent k-mers!!
        self.source = source # store prefix and suffix of k-mer to define node
        self.target = target


class DeBruijnGraph:
    def __init__(self):
        self.nodes = {} # string mapped to node object, EX: "AT": Node("AT")
        self.edges = [] # list of k-mer edges as strings EX: ATG, TGC, etc.

    def add_kmer(self, kmer):
        prefix_key = kmer[:-1]
        suffix_key = kmer[1:]

        if prefix_key not in self.nodes: # create nodes only if they don't already exist
            self.nodes[prefix_key] = Node(prefix_key)

        if suffix_key not in self.nodes:
            self.nodes[suffix_key] = Node(suffix_key)

        source = self.nodes[prefix_key] # call the node object
        target = self.nodes[suffix_key]

        
        edge = Edge(kmer, source, target) # pass created nodes to Edge()

        source.out_edges.append(edge)
        target.in_edges.append(edge) # add to both in an out edges since undirected
        self.edges.append(edge) # add edge to graph object to store (master edge list)

    def find_start_node(self):
        """
        Return the node an Eulerian path has to start at, or None when no
        Eulerian path exists.

        Degree conditions for a directed graph whose edges must each be used
        exactly once:

          * every node must satisfy in_degree == out_degree, except
          * exactly one node may have out_degree - in_degree == 1 (the start),
          * exactly one node may have in_degree - out_degree == 1 (the end).

        When every node is balanced the graph holds an Eulerian cycle, so any
        node that owns an outgoing edge can serve as the start.

        A graph whose degrees are balanced but whose edges fall into several
        connected components still has no single Eulerian path. That case is
        caught in eulerian_path(), because the traversal cannot consume every
        edge. Nodes left with no edges at all are ignored, they do not affect
        edge traversal.
        """

        start_candidates = [] # nodes with out_degree - in_degree == 1
        end_candidates = [] # nodes with in_degree - out_degree == 1

        for node in self.nodes.values(): # iterate over nodes
            out_degree = len(node.out_edges) # out degree
            in_degree = len(node.in_edges) # in degree
            degree_difference = out_degree - in_degree

            if degree_difference == 1:
                start_candidates.append(node)
            elif degree_difference == -1:
                end_candidates.append(node)
            elif degree_difference != 0: # difference is neither +1 nor -1
                return None # invalid configuration, no Eulerian path

        # More than one start or more than one end means the edges cannot be
        # traversed as a single path, even if every local degree looks legal.
        # (The number of start candidates always equals the number of end
        # candidates, because out_degree - in_degree sums to zero over all
        # nodes, so checking one side would already be enough.)
        if len(start_candidates) > 1 or len(end_candidates) > 1:
            return None

        if start_candidates: # exactly one start and one matching end
            return start_candidates[0]

        for node in self.nodes.values(): # balanced graph: Eulerian cycle
            if node.out_edges: # pick any node with an outgoing edge as start
                return node

        return None # graph has no edges at all


def eulerian_path(graph):
    start_node = graph.find_start_node() # return the start node
    path = [] # initialize path var
    
    if start_node is None:
        return None

    stack = [start_node] # create a stack (LIFO)

    unused = {} # track unused edges separate of graph traversed

    for node_key, node in graph.nodes.items():  # for every node, make a copy in the unused
        unused[node_key] = list(node.out_edges)

    while stack:
        current_node = stack[-1] # current node is last in stack

        if unused[current_node.node_key]: # if it exists in unused, use the edge
            edge = unused[current_node.node_key].pop() # pop edge using current node from unused
            stack.append(edge.target) # append target node reached

        else:
            path.append(stack.pop()) 
            # NOTE: ONLY append nodes to the stack when a deadend is reached. 
            # This ensures that the eulerian traversal will finish completely on a given path before backtracking 
            # and adding to the output path. By this method, we need to reverse the returned path.

    if len(path) != len(graph.edges) + 1: # eulerian path contains V-1 = E. IF more, return none
        return None # not every edge was traversed

    path.reverse()

    return path


def path_to_sequence(path):
    if not path: # no input
        return ""

    sequence = path[0].node_key # use label of first node in path

    for node in path[1:]: # iterate over path
        sequence += node.node_key[-1] # extend sequence

    return sequence
