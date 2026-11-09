"""
De Bruijn graph visualization for the Greedy Genome Gobblers project.

This module only *draws* what the assembly pipeline already produced. It does
not build graphs, does not traverse them, and does not assemble anything: every
node, edge, and path shown here is read out of the existing `DeBruijnGraph`
object, so the figures cannot drift away from the real data structure.

Two conventions match the rest of the project:

  * a node is a (k-1)-mer,
  * an edge is a k-mer, drawn as an arrow from its prefix to its suffix.

The Matplotlib backend is forced to "Agg" (non-interactive) so the figures can
be generated on a headless machine or in a test run. Nothing here opens a GUI
window; figures are written to disk when `output_path` is given.
"""

import heapq
import os

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, must be set before pyplot

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle

from genome import generate_kmer

# One colour per nucleotide, shared by every figure so that the same base always
# means the same colour within the presentation.
BASE_COLORS = {
    "A": "#4C78A8",  # blue
    "T": "#E45756",  # red
    "C": "#54A24B",  # green
    "G": "#EECA3B",  # yellow/gold
}

# Node fill colours by role in the Eulerian traversal.
NODE_FILL = {
    "start": "#D7F0DA",   # green tint  : valid Eulerian start
    "end": "#FBE0D1",     # orange tint : valid Eulerian end
    "inner": "#E8ECF4",   # neutral
    "isolated": "#F2F2F2",  # node without any edge
}

NODE_EDGE_COLOR = "#41506B"
UNUSED_EDGE_COLOR = "#AAB3C2"
NODE_WIDTH = 1.0
NODE_HEIGHT = 0.6
LEVEL_SPACING = 1.75   # horizontal distance between topological levels
STACK_SPACING = 1.55   # vertical distance between nodes in the same level


# ---------------------------------------------------------------------------
# Layout and element extraction (pure functions, no drawing)
# ---------------------------------------------------------------------------

def _compute_layout(graph):
    """
    Deterministic left-to-right layout for a De Bruijn graph.

    Nodes are placed by topological level (the longest path from a source), so
    a linear chain such as AT -> TG -> GC -> CG -> GT comes out as one straight
    row. Nodes sharing a level are stacked vertically in sorted order. Cycles
    cannot be topologically ordered, so the remaining nodes are appended after
    the resolved ones instead of being dropped.

    input:
    graph (DeBruijnGraph): the graph to lay out

    output:
    dict: {node_key: (x, y)}
    """

    node_keys = list(graph.nodes.keys())

    # Self-loops say nothing about ordering, so they are ignored here.
    successors = {
        key: sorted({
            edge.target.node_key
            for edge in graph.nodes[key].out_edges
            if edge.target.node_key != key
        })
        for key in node_keys
    }

    in_degree = {key: 0 for key in node_keys}
    for key in node_keys:
        for target_key in successors[key]:
            in_degree[target_key] += 1

    levels = {key: 0 for key in node_keys}
    ready = [key for key in node_keys if in_degree[key] == 0]
    heapq.heapify(ready)  # sorted processing keeps the layout reproducible

    resolved = set()
    while ready:
        key = heapq.heappop(ready)
        resolved.add(key)
        for target_key in successors[key]:
            levels[target_key] = max(levels[target_key], levels[key] + 1)
            in_degree[target_key] -= 1
            if in_degree[target_key] == 0:
                heapq.heappush(ready, target_key)

    unresolved = sorted(key for key in node_keys if key not in resolved)
    if unresolved:
        next_level = max(levels.values()) + 1 if levels else 0
        for key in unresolved:  # cyclic remainder, still placed deterministically
            levels[key] = next_level

    by_level = {}
    for key in node_keys:
        by_level.setdefault(levels[key], []).append(key)

    positions = {}
    for level, members in sorted(by_level.items()):
        members = sorted(members)
        middle = (len(members) - 1) / 2.0
        for index, key in enumerate(members):
            # Highest key on top, centred around y = 0.
            positions[key] = (
                level * LEVEL_SPACING,
                (middle - index) * STACK_SPACING,
            )

    return positions


def _classify_nodes(graph):
    """
    Label every node by its role in an Eulerian traversal, using only degrees.

    A node with out_degree - in_degree == 1 can only start a path, a node with
    in_degree - out_degree == 1 can only end one. Several candidates of the same
    role means the graph has no single Eulerian path, but the figure still shows
    every candidate so the problem is visible instead of hidden.

    output:
    dict: {node_key: "start" | "end" | "inner" | "isolated"}
    """

    roles = {}
    for key, node in graph.nodes.items():
        out_degree = len(node.out_edges)
        in_degree = len(node.in_edges)

        if out_degree - in_degree == 1:
            roles[key] = "start"
        elif in_degree - out_degree == 1:
            roles[key] = "end"
        elif out_degree == 0 and in_degree == 0:
            roles[key] = "isolated"
        else:
            roles[key] = "inner"
    return roles


def _edge_curvatures(graph):
    """
    Assign one distinct curvature to every edge so parallel edges stay visible.

    Repeated k-mers produce several Edge objects with identical endpoints. If
    they were drawn on top of each other only one would be visible, so each
    parallel edge gets its own bow. Self-loops are flagged for special drawing.

    output:
    dict: {id(edge): {"rad": float, "self_loop": bool}}
    """

    groups = {}
    for edge in graph.edges:
        groups.setdefault((edge.source.node_key, edge.target.node_key), []).append(edge)

    styles = {}
    for (source_key, target_key), group in groups.items():
        is_self_loop = source_key == target_key
        count = len(group)
        for index, edge in enumerate(group):
            if is_self_loop:
                # A loop from a node to itself never overlaps another loop of the
                # same node unless the k-mer repeats, so grow the radius instead.
                styles[id(edge)] = {"rad": 1.4 + 0.35 * index, "self_loop": True}
            elif count == 1:
                styles[id(edge)] = {"rad": 0.0, "self_loop": False}
            else:
                # Symmetric fan: -0.35, 0, +0.35 for three parallel edges, etc.
                offset = index - (count - 1) / 2.0
                styles[id(edge)] = {"rad": offset * 0.45, "self_loop": False}

    return styles


def graph_elements(graph):
    """
    Extract everything the figures need, straight from the graph object.

    input:
    graph (DeBruijnGraph): graph to describe

    output:
    dict with keys:
        "nodes": list of {key, label, position, role}
        "edges": list of {edge, kmer, source, target, rad, self_loop}
        "node_count": int
        "edge_count": int
        "start_keys": list[str]
        "end_keys": list[str]
    """

    positions = _compute_layout(graph)
    roles = _classify_nodes(graph)
    curvatures = _edge_curvatures(graph)

    nodes = []
    for key, node in graph.nodes.items():
        nodes.append({
            "key": key,
            "label": node.node_key,
            "position": positions[key],
            "role": roles[key],
        })

    edges = []
    for edge in graph.edges:
        style = curvatures[id(edge)]
        edges.append({
            "edge": edge,
            "kmer": edge.kmer,
            "source": edge.source.node_key,
            "target": edge.target.node_key,
            "rad": style["rad"],
            "self_loop": style["self_loop"],
        })

    return {
        "nodes": nodes,
        "edges": edges,
        "node_count": len(nodes),
        "edge_count": len(edges),
        "start_keys": [n["key"] for n in nodes if n["role"] == "start"],
        "end_keys": [n["key"] for n in nodes if n["role"] == "end"],
    }


def ordered_edges_for_path(graph, path):
    """
    Map an Eulerian node path onto the actual Edge objects it walks over.

    `eulerian_path()` returns Node objects, not edges, so each step is resolved
    by looking for an unused outgoing edge of the current node whose target is
    the next node. Edges are consumed from a private copy of each node's
    out_edges list, which means:

      * repeated k-mers are never reused (the second identical step picks the
        second identical Edge instance),
      * self-loops are resolved normally (source and target are the same node),
      * the graph object itself is never modified.

    input:
    graph (DeBruijnGraph): the graph the path was computed on
    path (list[Node] | None): node sequence returned by eulerian_path()

    output:
    list[Edge] | None: edge instances in traversal order, or None when the path
    is missing or a step has no matching unused edge
    """

    if not path or len(path) < 2:
        return []

    available = {
        key: list(node.out_edges) for key, node in graph.nodes.items()
    }

    ordered = []
    for current_node, next_node in zip(path, path[1:]):
        bucket = available.get(current_node.node_key, [])
        chosen_index = None
        for index, edge in enumerate(bucket):
            if edge.target.node_key == next_node.node_key:
                chosen_index = index
                break

        if chosen_index is None:
            return None  # path does not correspond to unused edges of this graph

        ordered.append(bucket.pop(chosen_index))

    return ordered


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _text_color(base):
    """Pick a readable label colour for a nucleotide block."""

    return "#FFFFFF" if base in ("A", "T", "C") else "#202020"


def _draw_colored_sequence(ax, x, y, sequence, fontsize=11, weight="bold"):
    """
    Draw a short string one nucleotide at a time, each base in its own colour.

    Characters are spaced with display-space offsets (points) rather than data
    units, so the label stays centred and legible whatever the axis scale is.
    """

    advance_points = 0.62 * fontsize
    centre = (len(sequence) - 1) / 2.0
    for index, base in enumerate(sequence):
        ax.annotate(
            base,
            xy=(x, y),
            xytext=((index - centre) * advance_points, 0),
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=fontsize,
            fontweight=weight,
            family="monospace",
            color=BASE_COLORS.get(base, "#202020"),
        )


def _perpendicular(source, target):
    """Unit perpendicular vector for an edge, used to offset labels."""

    dx = target[0] - source[0]
    dy = target[1] - source[1]
    length = (dx * dx + dy * dy) ** 0.5
    if length == 0:
        return (0.0, 1.0)
    return (-dy / length, dx / length)


def _edge_arrow(ax, start, end, rad, color, linewidth=1.6, alpha=1.0,
                arrow_size=13, zorder=2):
    """Draw one directed edge as a curved arrow."""

    arrow = FancyArrowPatch(
        start,
        end,
        connectionstyle="arc3,rad={}".format(rad),
        arrowstyle="-|>",
        mutation_scale=arrow_size,
        linewidth=linewidth,
        color=color,
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(arrow)


def _shorten(start, end, offset):
    """Move both endpoints `offset` data units along the edge direction."""

    dx = end[0] - start[0]
    dy = end[1] - start[1]
    length = (dx * dx + dy * dy) ** 0.5
    if length <= 2 * offset or length == 0:
        return start, end
    ux, uy = dx / length, dy / length
    return (
        (start[0] + ux * offset, start[1] + uy * offset),
        (end[0] - ux * offset, end[1] - uy * offset),
    )


def _draw_node(ax, x, y, label, role, alpha=1.0):
    """Draw one (k-1)-mer as a rounded box with per-base coloured text."""

    box = FancyBboxPatch(
        (x - NODE_WIDTH / 2.0, y - NODE_HEIGHT / 2.0),
        NODE_WIDTH,
        NODE_HEIGHT,
        boxstyle="round,pad=0.02,rounding_size=0.14",
        facecolor=NODE_FILL.get(role, NODE_FILL["inner"]),
        edgecolor=NODE_EDGE_COLOR,
        linewidth=1.4,
        alpha=alpha,
        zorder=3,
    )
    ax.add_patch(box)
    _draw_colored_sequence(ax, x, y, label, fontsize=11)


def _draw_graph_base(ax, elements, alpha=1.0, draw_nodes=True):
    """Draw every node and every edge of the graph, edges underneath nodes."""

    for item in elements["edges"]:
        source = elements_position(elements, item["source"])
        target = elements_position(elements, item["target"])

        if item["self_loop"]:
            start = (source[0] - 0.22, source[1] + NODE_HEIGHT / 2.0)
            end = (source[0] + 0.22, source[1] + NODE_HEIGHT / 2.0)
            _edge_arrow(ax, start, end, item["rad"], UNUSED_EDGE_COLOR,
                        linewidth=1.4, alpha=alpha)
        else:
            start, end = _shorten(source, target, 0.64)
            _edge_arrow(ax, start, end, item["rad"], UNUSED_EDGE_COLOR,
                        linewidth=1.4, alpha=alpha)

    if not draw_nodes:
        return

    for node in elements["nodes"]:
        x, y = node["position"]
        _draw_node(ax, x, y, node["label"], node["role"], alpha=alpha)


def elements_position(elements, key):
    """Look up a node position inside an element dictionary."""

    for node in elements["nodes"]:
        if node["key"] == key:
            return node["position"]
    raise KeyError(key)


def _edge_label_point(source, target, rad, self_loop, lift=0.0):
    """Where to place an edge label so it does not sit on the arrow."""

    if self_loop:
        return (source[0], source[1] + NODE_HEIGHT / 2.0 + 0.85 + lift)

    midpoint = ((source[0] + target[0]) / 2.0, (source[1] + target[1]) / 2.0)
    perpendicular = _perpendicular(source, target)
    distance = (
        (target[0] - source[0]) ** 2 + (target[1] - source[1]) ** 2
    ) ** 0.5
    curve_lift = rad * distance / 2.0
    offset = 0.26 + curve_lift + lift
    return (
        midpoint[0] + perpendicular[0] * offset,
        midpoint[1] + perpendicular[1] * offset,
    )


def _figure_size(elements):
    """Figure size that grows with the graph instead of cramping it."""

    xs = [node["position"][0] for node in elements["nodes"]]
    ys = [node["position"][1] for node in elements["nodes"]]
    width = max(6.0, (max(xs) - min(xs)) * 0.95 + 3.0) if xs else 6.0
    height = max(3.4, (max(ys) - min(ys)) * 0.85 + 2.6) if ys else 3.4
    return width, height


def _finish(ax, elements, title, subtitle=None):
    """Apply common axis styling and the title block."""

    xs = [node["position"][0] for node in elements["nodes"]]
    ys = [node["position"][1] for node in elements["nodes"]]

    margin_x = 1.6
    margin_y = 1.4
    if xs and ys:
        ax.set_xlim(min(xs) - margin_x, max(xs) + margin_x)
        ax.set_ylim(min(ys) - margin_y, max(ys) + margin_y)
    ax.set_aspect("equal")
    ax.axis("off")

    # Generous title padding leaves room for the subtitle line underneath it.
    ax.set_title(title, fontsize=13, fontweight="bold", pad=26)
    if subtitle:
        ax.text(
            0.5, 1.004, subtitle,
            transform=ax.transAxes,
            ha="center", va="bottom",
            fontsize=9.5, color="#444444",
        )


def _save(fig, output_path):
    """Save and close a figure, returning the path or the figure itself."""

    if output_path is None:
        return fig

    directory = os.path.dirname(os.path.abspath(output_path))
    if directory:
        os.makedirs(directory, exist_ok=True)

    fig.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output_path


# ---------------------------------------------------------------------------
# Public figures
# ---------------------------------------------------------------------------

def visualize_debruijn_graph(graph, output_path=None):
    """
    Draw the De Bruijn graph: one box per (k-1)-mer, one arrow per k-mer.

    Parallel k-mers get separate bows and self-loops are drawn as loops, so no
    edge is hidden. Valid Eulerian start and end nodes are highlighted and
    labelled in the legend.

    input:
    graph (DeBruijnGraph): graph to draw
    output_path (str | None): PNG path; when omitted the figure is returned

    output:
    str | Figure: saved path, or the Matplotlib figure when nothing is saved
    """

    elements = graph_elements(graph)
    fig, ax = plt.subplots(figsize=_figure_size(elements))

    _draw_graph_base(ax, elements)

    for item in elements["edges"]:
        source = elements_position(elements, item["source"])
        target = elements_position(elements, item["target"])
        label_x, label_y = _edge_label_point(
            source, target, item["rad"], item["self_loop"]
        )
        ax.text(
            label_x, label_y, item["kmer"],
            ha="center", va="center",
            fontsize=9, color="#1F2A44", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor="#D3D9E2", linewidth=0.7),
            zorder=5,
        )

    # Node keys are already visible inside the boxes; annotate the roles only in
    # the legend so the figure stays readable.
    handles = [
        Patch(facecolor=NODE_FILL["start"], edgecolor=NODE_EDGE_COLOR,
              label="Eulerian start node"),
        Patch(facecolor=NODE_FILL["end"], edgecolor=NODE_EDGE_COLOR,
              label="Eulerian end node"),
        Patch(facecolor=NODE_FILL["inner"], edgecolor=NODE_EDGE_COLOR,
              label="internal node"),
    ]
    ax.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
              fontsize=8.5, bbox_to_anchor=(0.5, -0.12))

    subtitle = "nodes = (k-1)-mers, edges = k-mers ({})".format(
        "start {} / end {}".format(
            ", ".join(elements["start_keys"]) or "none",
            ", ".join(elements["end_keys"]) or "none",
        )
    )
    _finish(
        ax,
        elements,
        "De Bruijn graph ({} nodes, {} edges)".format(
            elements["node_count"], elements["edge_count"]
        ),
        subtitle,
    )

    return _save(fig, output_path)


def visualize_eulerian_path(graph, path, output_path=None):
    """
    Draw the graph together with its Eulerian traversal.

    Traversed edges are drawn on top of the grey background graph in traversal
    order, each numbered, and each in its own shade so the direction of travel
    is obvious. The reconstructed sequence is put in the title block.

    Nothing is recomputed: the ordered edges come from `ordered_edges_for_path`,
    which only reads the graph. A missing or invalid path is reported in the
    figure instead of raising.

    input:
    graph (DeBruijnGraph): graph the path belongs to
    path (list[Node] | None): output of eulerian_path()
    output_path (str | None): PNG path; when omitted the figure is returned

    output:
    str | Figure: saved path, or the Matplotlib figure when nothing is saved
    """

    from de_brujin import path_to_sequence  # local import avoids a cycle

    elements = graph_elements(graph)
    fig, ax = plt.subplots(figsize=_figure_size(elements))

    # Nodes are drawn twice on purpose: a faded copy belongs to the background
    # graph, the solid copy is painted last so arrows never cover the labels.
    # Skipping the faded node labels avoids duplicate text at the same spot.
    _draw_graph_base(ax, elements, alpha=0.55, draw_nodes=False)

    ordered_edges = ordered_edges_for_path(graph, path)

    if not path or len(path) < 2:
        ax.text(0.5, -0.16, "No Eulerian path: the graph cannot be traversed as "
                            "a single trail over all edges.",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=10, color="#B03A2E", fontweight="bold")
        _finish(ax, elements, "Eulerian traversal: no valid path")
        return _save(fig, output_path)

    sequence = path_to_sequence(path)

    if ordered_edges is None:
        ax.text(0.5, -0.16, "The node path does not match the edges of this "
                            "graph; traversal could not be resolved.",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=10, color="#B03A2E", fontweight="bold")
        _finish(ax, elements, "Eulerian traversal: unresolved path")
        return _save(fig, output_path)

    positions = {node["key"]: node["position"] for node in elements["nodes"]}
    styles = _edge_curvatures(graph)
    step_count = len(ordered_edges)
    colormap = plt.get_cmap("viridis")

    for step, edge in enumerate(ordered_edges, start=1):
        source = positions[edge.source.node_key]
        target = positions[edge.target.node_key]
        style = styles[id(edge)]
        color = colormap(step / max(step_count, 1))

        if style["self_loop"]:
            start = (source[0] - 0.22, source[1] + NODE_HEIGHT / 2.0)
            end = (source[0] + 0.22, source[1] + NODE_HEIGHT / 2.0)
        else:
            start, end = _shorten(source, target, 0.64)

        _edge_arrow(ax, start, end, style["rad"], color,
                    linewidth=2.6, arrow_size=15, zorder=4)

        label_x, label_y = _edge_label_point(
            source, target, style["rad"], style["self_loop"], lift=0.30
        )
        ax.text(
            label_x, label_y, "{}. {}".format(step, edge.kmer),
            ha="center", va="center",
            fontsize=9, color="#FFFFFF", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.22", facecolor=color,
                      edgecolor="#FFFFFF", linewidth=0.8),
            zorder=6,
        )

    # Node labels are redrawn on top so the traversal arrows do not cover them.
    for node in elements["nodes"]:
        x, y = node["position"]
        _draw_node(ax, x, y, node["label"], node["role"], alpha=1.0)

    ax.text(
        0.5, -0.14,
        "reconstructed sequence: {}   ({} of {} edges traversed)".format(
            sequence, step_count, elements["edge_count"]
        ),
        transform=ax.transAxes, ha="center", va="top",
        fontsize=10.5, color="#1F2A44", fontweight="bold",
    )

    _finish(
        ax,
        elements,
        "Eulerian traversal reconstructs the DNA sequence",
        "edge numbers follow the traversal order; grey edges are not part of "
        "this path",
    )

    return _save(fig, output_path)


def visualize_sequence_reconstruction(reference_genome, reads, kmers,
                                      assembled_sequence, k=None,
                                      output_path=None):
    """
    Educational figure of the whole pipeline: reference, reads, k-mers, assembly.

    Reads are placed at their position on the reference only when that position
    is unambiguous. A read that occurs several times, or not at all, is drawn in
    a separate panel and labelled as ambiguous rather than assigned a made-up
    origin.

    input:
    reference_genome (str): ground-truth sequence
    reads (list[str]): simulated reads
    kmers (list[str]): flat k-mer list produced from the reads
    assembled_sequence (str | None): sequence rebuilt from the graph
    k (int | None): k-mer length, inferred from the first k-mer when omitted
    output_path (str | None): PNG path; when omitted the figure is returned

    output:
    str | Figure: saved path, or the Matplotlib figure when nothing is saved
    """

    reads = list(reads)
    kmers = list(kmers)

    if k is None:
        k = len(kmers[0]) if kmers else 1
    if reads and k > len(reads[0]):
        raise ValueError("k is larger than the shortest read")

    n_bases = len(reference_genome)
    ambiguous_x = n_bases + 1.6   # panel for reads with no unique position
    block_height = 0.62

    # --- decide whether each read has one unambiguous position --------------
    placements = []
    for read in reads:
        start = reference_genome.find(read)
        if start == -1:
            placements.append({"read": read, "start": None,
                               "status": "not found in reference"})
            continue

        second = reference_genome.find(read, start + 1)
        if second != -1:
            placements.append({"read": read, "start": None,
                               "status": "ambiguous, appears several times"})
        else:
            placements.append({"read": read, "start": start, "status": "placed"})

    # --- build rows ---------------------------------------------------------
    rows = []  # (label, y, kind)
    y_cursor = 0.0

    def add_row(label, kind, height=1.0):
        nonlocal y_cursor
        y_cursor -= height
        rows.append({"label": label, "y": y_cursor, "kind": kind})
        return y_cursor

    reference_y = add_row("Reference genome", "reference")

    read_rows = []
    for index, placement in enumerate(placements, start=1):
        if placement["status"] == "placed":
            label = "Read {} (position {})".format(index, placement["start"])
        else:
            label = "Read {} ({})".format(index, placement["status"])
        y = add_row(label, "read")
        read_rows.append({"placement": placement, "y": y})

    kmer_rows = []
    for index, entry in enumerate(read_rows, start=1):
        y = add_row("k-mers of read {}".format(index), "kmer", height=1.15)
        kmer_rows.append({"placement": entry["placement"], "y": y})

    assembly_y = add_row("Reconstructed genome", "assembly")

    # --- draw ---------------------------------------------------------------
    fig_height = max(4.4, 0.95 * len(rows) + 1.9)
    fig, ax = plt.subplots(figsize=(max(7.6, 0.52 * n_bases + 4.2), fig_height))

    def draw_blocks(x_start, y, sequence, fontsize=10):
        for offset, base in enumerate(sequence):
            ax.add_patch(Rectangle(
                (x_start + offset, y - block_height / 2.0), 1.0, block_height,
                facecolor=BASE_COLORS.get(base, "#CCCCCC"),
                edgecolor="white", linewidth=1.0, zorder=2,
            ))
            ax.text(
                x_start + offset + 0.5, y, base,
                ha="center", va="center", fontsize=fontsize,
                color=_text_color(base), family="monospace",
                fontweight="bold", zorder=3,
            )

    # 1) reference genome, with base indices
    draw_blocks(0, reference_y, reference_genome, fontsize=10.5)
    for index in range(n_bases + 1):
        ax.text(index, reference_y - 0.62, str(index), ha="center", va="top",
                fontsize=6.5, color="#8A93A3")

    # 2) reads, placed only when unambiguous
    ambiguous_stack = 0
    for entry in read_rows:
        placement = entry["placement"]
        y = entry["y"]
        if placement["status"] == "placed":
            draw_blocks(placement["start"], y, placement["read"])
        else:
            draw_blocks(ambiguous_x, y, placement["read"])
            ax.text(
                ambiguous_x + len(placement["read"]) + 0.25, y,
                placement["status"],
                ha="left", va="center", fontsize=8.5,
                color="#B03A2E", style="italic",
            )
            ambiguous_stack += 1

    # 3) k-mers, drawn as a staircase to show the (k-1) overlap
    for entry in kmer_rows:
        placement = entry["placement"]
        base_x = ambiguous_x if placement["status"] != "placed" else placement["start"]
        read_kmers = generate_kmer(placement["read"], k)
        for index, kmer in enumerate(read_kmers):
            y = entry["y"] - index * 0.16
            ax.add_patch(Rectangle(
                (base_x + index, y - block_height / 2.0), k, block_height,
                facecolor="#F4F6FA", edgecolor=BASE_COLORS.get(kmer[0], "#8899AA"),
                linewidth=1.0, zorder=2,
            ))
            for offset, base in enumerate(kmer):
                ax.text(
                    base_x + index + offset + 0.5, y, base,
                    ha="center", va="center", fontsize=7.5,
                    color=BASE_COLORS.get(base, "#202020"),
                    family="monospace", fontweight="bold", zorder=3,
                )

    # 4) reconstructed genome
    if assembled_sequence:
        draw_blocks(0, assembly_y, assembled_sequence, fontsize=10.5)
    else:
        ax.text(0.05, assembly_y, "no sequence produced (graph has no valid "
                                  "Eulerian path)",
                ha="left", va="center", fontsize=9.5,
                color="#B03A2E", style="italic")

    matches = assembled_sequence == reference_genome

    # k-mer list summary under the figure
    summary = "k = {}   k-mers: {}".format(k, ", ".join(kmers))
    ax.text(0, y_cursor - 0.95, summary, ha="left", va="center",
            fontsize=8.5, color="#444444")

    verdict = ("reconstructed sequence equals the reference genome" if matches
               else "reconstructed sequence differs from the reference genome")
    ax.text(0, y_cursor - 1.35, verdict, ha="left", va="center",
            fontsize=9.5, color="#1F6F3F" if matches else "#B03A2E",
            fontweight="bold")

    # Row names are drawn as y tick labels, which Matplotlib lays out outside the
    # axes and keeps inside the saved image, so no manual column offset is used.
    ax.set_yticks([row["y"] for row in rows])
    ax.set_yticklabels([row["label"] for row in rows], fontsize=9.5)
    for tick_label, row in zip(ax.get_yticklabels(), rows):
        if row["kind"] in ("reference", "assembly"):
            tick_label.set_fontweight("bold")
        tick_label.set_color("#1F2A44")

    ax.tick_params(axis="y", left=False, length=0)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlim(-0.6, max(n_bases, ambiguous_x + 4.0) + 1.0)
    ax.set_ylim(y_cursor - 1.7, 0.9)

    ax.legend(
        handles=[Patch(facecolor=BASE_COLORS[base], label=base) for base in "ATCG"],
        loc="lower right", ncol=4, frameon=False, fontsize=8.5,
        bbox_to_anchor=(1.0, -0.04),
    )

    ax.set_title(
        "From short reads to reconstructed genome",
        fontsize=13, fontweight="bold", pad=26,
    )
    ax.text(
        0.0, 1.012,
        "reference {}  |  {} reads  |  k = {}".format(
            reference_genome, len(reads), k
        ),
        transform=ax.transAxes, ha="left", va="bottom",
        fontsize=9, color="#444444",
    )

    return _save(fig, output_path)


def create_all_figures(graph, path, reference_genome, reads, kmers,
                       assembled_sequence, output_directory="figures",
                       k=None):
    """
    Write the three milestone figures into one directory.

    input:
    graph (DeBruijnGraph): assembled graph
    path (list[Node] | None): Eulerian path used for the assembly
    reference_genome (str): ground truth
    reads (list[str]): reads that produced the graph
    kmers (list[str]): flat k-mer list
    assembled_sequence (str | None): reconstructed sequence
    output_directory (str): folder that receives the PNG files
    k (int | None): k-mer length, inferred when omitted

    output:
    list[str]: absolute paths of the generated PNG files, in figure order
    """

    if k is None:
        k = len(kmers[0]) if kmers else 1

    written = []
    written.append(visualize_debruijn_graph(
        graph, os.path.join(output_directory, "debruijn_graph.png")
    ))
    written.append(visualize_eulerian_path(
        graph, path, os.path.join(output_directory, "eulerian_traversal.png")
    ))
    written.append(visualize_sequence_reconstruction(
        reference_genome, reads, kmers, assembled_sequence, k=k,
        output_path=os.path.join(output_directory, "sequence_reconstruction.png"),
    ))
    return written
