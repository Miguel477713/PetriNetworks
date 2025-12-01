import os, numpy, math
from graphviz import Digraph
from typing import Optional
from collections import deque, defaultdict

if os.name == 'nt':
    os.system('cls')  # Clear the console for Windows
elif os.name == 'posix':
    os.system('clear')  # Clear the console for Linux/Mac
    
# New: build reachability graph and draw with graphviz
def build_reachability_graph(pre: numpy.ndarray, post: numpy.ndarray, start_marking: numpy.ndarray):
    """
    Explore reachable markings using Karp-Miller coverability tree algorithm.
    Returns (nodes_map, edges, depths)
    """
    adjacency = post - pre
    num_transitions = pre.shape[1]
    # Convert start to tuple because numpy arrays are not hashable
    start = tuple(int(x) for x in start_marking)

    # nodes_map: Maps a Marking Tuple -> Unique ID or same Tuple
    nodes = {start: start}
    # Queue items: (current_marking, list_of_ancestors_on_this_path)
    queue = deque([(start, [])]) # BFS queue
    edges = []
    depths = {start: 0}

    # Helper to check if edge exists
    def edge_exists(edgelist, a, b, t):
        return (a, b, t) in edgelist

    # Helper to add values dealing with infinity
    def add_markings(m, change):
        res = [] # result marking
        for i in range(len(m)):
            if m[i] is math.inf: # we use math.inf to represent ω since numpy.inf can cause issues in comparisons
                res.append(math.inf) # ω + x = ω
            else:
                res.append(m[i] + change[i]) # normal addition
        return tuple(res)

    while queue: # while there are markings to explore
        marking, ancestors = queue.popleft() # we use popleft for BFS since we want shortest paths, equivalent to pop(0)

        for t in range(num_transitions):
            # 1. Check enabled: marking >= pre[:, t]
            enabled = True
            for i, req in enumerate(pre[:, t].astype(int)): # for each place check if enough tokens
                mi = marking[i] # current tokens in place i
                if mi is math.inf: # omega can always satisfy requirements
                    continue
                if mi < req: # not enough tokens
                    enabled = False
                    break # no need to check further
            if not enabled:
                continue

            change_vec = adjacency[:, t] # change vector for transition t
            new_m = add_markings(marking, change_vec) # compute new marking with addition

            # Check against ALL ancestors on this specific path
            current_path_ancestors = ancestors + [marking]

            # We might need to apply omega multiple times if it dominates multiple ancestors
            temp_m_list = list(new_m)
            for ancestor in current_path_ancestors:
                # Check if new_m >= ancestor (component-wise)
                is_ge = True # is greater or equal
                strictly_greater = False # is strictly greater in at least one component

                for k in range(len(temp_m_list)): # for each place
                    val_new = temp_m_list[k] # new marking value
                    # Assign the value from the ancestor marking to a variable for clarity
                    val_anc = ancestor[k]

                    if val_anc is math.inf and val_new is not math.inf: # we compare if ancestor has ω and new does not
                        is_ge = False; # new cannot be >= ω if it is not ω
                        break
                    if val_new is not math.inf and val_anc is not math.inf and val_new < val_anc: # both finite and new < ancestor
                        is_ge = False;# we set is_ge to false
                        break
                    if val_new != val_anc: # if they are different and new is not ω
                        strictly_greater = True # we have at least one strictly greater component

                if is_ge and strictly_greater:
                    # Replace growing components with infinity
                    for k in range(len(temp_m_list)):
                        if temp_m_list[k] is not math.inf and ancestor[k] is not math.inf:
                            if temp_m_list[k] > ancestor[k]: # strictly greater component
                                temp_m_list[k] = math.inf # set to ω

            new_m = tuple(temp_m_list)

            if not edge_exists(edges, marking, new_m, t + 1):
                edges.append((marking, new_m, t + 1)) # Add edge (1-based transition index)

            if new_m not in nodes: # here we will add new nodes only if they are unique
                nodes[new_m] = new_m
                depths[new_m] = depths[marking] + 1
                queue.append((new_m, current_path_ancestors))
            else:
                # It's a duplicate state (or previously visited), so we linked it (Edge added above)
                # but we do NOT add it to queue to avoid infinite loops
                pass

    return nodes, edges, depths


# Tarjan's algorithm to find strongly connected components
def tarjan_scc(adj_map: dict):
    """Return list of SCCs where each SCC is a list of nodes (node keys).
    adj_map: dict[node] -> list of neighbor nodes
    """
    index = {}
    lowlink = {}
    onstack = set()
    stack = []
    sccs = []
    next_index = 0

    def strongconnect(v):
        nonlocal next_index
        index[v] = next_index
        lowlink[v] = next_index
        next_index += 1
        stack.append(v)
        onstack.add(v)

        for w in adj_map.get(v, []):
            if w not in index:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in onstack:
                lowlink[v] = min(lowlink[v], index[w])

        # If v is a root node, pop the stack and generate an SCC
        if lowlink[v] == index[v]:
            scc = []
            while True:
                w = stack.pop()
                onstack.remove(w)
                scc.append(w)
                if w == v:
                    break
            sccs.append(scc)

    for v in adj_map.keys():
        if v not in index:
            strongconnect(v)

    return sccs


def analyze_properties(nodes: dict, edges: list, start_marking: tuple, num_transitions: int):
    """
    Analyze Boundedness, Liveness, and Reversibility based on the generated graph.
    References definitions from '11-Redes de Petri.pdf'.
    """
    print("\n" + "=" * 60)
    print("QUALITATIVE ANALYSIS (Based on course notes: 11-Redes de Petri)")
    print("=" * 60)

    # Helper: format marking elements to Python ints or 'ω'
    def format_marking(m):
        return tuple('ω' if x is math.inf else int(x) for x in m)

    # ---------------------------------------------------------
    # 1. BOUNDEDNESS ANALYSIS
    # ---------------------------------------------------------
    is_bounded = True
    omega_node = None

    for m in nodes.keys():
        if math.inf in m:
            is_bounded = False
            omega_node = format_marking(m)
            break

    print(f"\n1. BOUNDEDNESS: {'BOUNDED' if is_bounded else 'UNBOUNDED'}")
    if is_bounded:
        print(f"   The symbol ω does not appear in the generated graph.")
    else:
        print(f"   Reason: UNBOUNDED because ω appears in marking {omega_node}.")

    # Build adjacency list for graph traversal (needed for Liveness/Reversibility)
    adj = defaultdict(list)
    for from_marking, to_marking, t_idx in edges:
        adj[from_marking].append((to_marking, t_idx))

    # Also build a neighbor-only map for SCC computation
    neighbor_map = {n: [nbr for nbr, _ in adj[n]] for n in adj}
    # ensure isolated nodes are present
    for n in nodes:
        neighbor_map.setdefault(n, [])

    # ---------------------------------------------------------
    # 1.5 STRONGLY CONNECTED COMPONENTS (Tarjan)
    # ---------------------------------------------------------
    sccs = tarjan_scc(neighbor_map)
    print(f"\n1.5 STRONGLY CONNECTED COMPONENTS (Tarjan) - total: {len(sccs)}")
    for i, scc in enumerate(sccs, start=1):
        formatted = [format_marking(x) for x in scc]
        print(f"   SCC {i}: {formatted}")

    # ---------------------------------------------------------
    # Terminal SCC computation (for liveness and reversibility reasoning)
    # ---------------------------------------------------------
    # Map each node to its strongly connected component index
    strongly_connected_component_index_by_node = {}
    for component_index, component_nodes in enumerate(sccs):
        for node in component_nodes:
            strongly_connected_component_index_by_node[node] = component_index

    # Build condensation graph of strongly connected components
    strongly_connected_component_successors = {
        component_index: set() for component_index in range(len(sccs))
    }

    for from_marking, to_marking, t_idx in edges:
        component_index_from = strongly_connected_component_index_by_node[from_marking]
        component_index_to = strongly_connected_component_index_by_node[to_marking]
        if component_index_from != component_index_to:
            strongly_connected_component_successors[component_index_from].add(component_index_to)

    # Terminal strongly connected components: no outgoing edges in condensation graph
    terminal_strongly_connected_component_indices = [
        component_index
        for component_index, successor_set in strongly_connected_component_successors.items()
        if len(successor_set) == 0
    ]

    print("\nTerminal strongly connected components (no outgoing edges in condensation graph):")
    if not terminal_strongly_connected_component_indices:
        print("   None (every SCC has a successor SCC).")
    else:
        for component_index in terminal_strongly_connected_component_indices:
            component_nodes = sccs[component_index]
            formatted = [format_marking(x) for x in component_nodes]
            # +1 so the numbering matches the earlier SCC print (SCC 1, 2, 3, ...)
            print(f"   SCC {component_index + 1}: {formatted}")

    # ---------------------------------------------------------
    # 2. LIVENESS ANALYSIS (via terminal strongly connected components)
    # ---------------------------------------------------------
    deadlock_nodes = []
    for node in nodes:
        if not adj[node]:  # No outgoing edges means no transition is enabled
            deadlock_nodes.append(format_marking(node))

    print(f"\n2. LIVENESS: ", end="")

    if deadlock_nodes:
        # Deadlocks are terminal strongly connected components with no internal transitions.
        print("NOT LIVE (Deadlock Detected)")
        print(f"   Reason: The system contains deadlocks at marking(s): {deadlock_nodes}.")
        is_live = False
    else:
        # Check strict liveness using terminal SCCs:
        # For a live net, for every transition t and for every terminal SCC,
        # t must be enabled (there must be an edge labeled t) in at least one marking of that SCC.
        all_trans_live = True
        failed_trans = None
        failed_component_example_marking = None
        failed_component_index_for_display = None

        for t_idx in range(1, num_transitions + 1):
            transition_live_in_all_terminal_components = True

            for component_index in terminal_strongly_connected_component_indices:
                component_nodes = sccs[component_index]
                transition_occurs_inside_component = False

                # Look for an edge labeled t_idx that stays inside this SCC
                for from_marking, to_marking, edge_t in edges:
                    if edge_t == t_idx and from_marking in component_nodes and to_marking in component_nodes:
                        transition_occurs_inside_component = True
                        break

                if not transition_occurs_inside_component:
                    transition_live_in_all_terminal_components = False
                    failed_trans = t_idx
                    failed_component_index_for_display = component_index + 1  # match printed SCC numbering
                    break

            if not transition_live_in_all_terminal_components:
                all_trans_live = False
                break

        if all_trans_live:
            print("LIVE")
            print("   Reason: The net is deadlock-free AND every transition appears in all terminal strongly connected components.")
        else:
            print("NOT LIVE (Partial Functioning)")
            print(f"   Reason: Transition t{failed_trans} is not live according to the terminal strongly connected components.")
            print("   Explanation:")
            print("         Once the reachability graph enters terminal SCC {failed_component_index_for_display},")
            print(f"       transition t{failed_trans} can never fire again from any marking in that region.")
            is_live = False

    # ---------------------------------------------------------
    # 3. REVERSIBILITY ANALYSIS (via strongly connected components)
    # ---------------------------------------------------------
    # Reversibility: for every reachable marking M, M0 must be reachable from M.
    # On the reachability graph, this is equivalent to:
    #   "All reachable markings lie in the same strongly connected component as M0."
    print(f"\n3. REVERSIBILITY: ", end="")

    # Identify the SCC that contains the initial marking
    component_index_initial = strongly_connected_component_index_by_node[start_marking]
    strongly_connected_component_of_initial_marking = sccs[component_index_initial]

    if len(strongly_connected_component_of_initial_marking) == len(nodes):
        is_reversible = True
        irreversible_node = None
    else:
        is_reversible = False
        # Find an example marking that is not in the SCC of the initial marking
        example_irreversible_marking = None
        for node in nodes:
            if node not in strongly_connected_component_of_initial_marking:
                example_irreversible_marking = node
                break
        irreversible_node = format_marking(example_irreversible_marking)

    if is_reversible:
        print("REVERSIBLE")
        print("   Reason: All reachable markings belong to the same strongly connected component as the initial marking M0.")
    else:
        print("NOT REVERSIBLE")
        print("   Reason: The coverability graph is not strongly connected.")
        print(f"   The initial marking M0 belongs to SCC {component_index_initial + 1},")
        print("   but there exist other strongly connected components that do not contain M0.")
        print(f"   Example: From marking {irreversible_node} it is not possible to reach the initial marking M0.")

    print("=" * 60 + "\n")


def draw_reachability_graph(nodes_map, edges, depths, filename='reachability_graph'):
    """
    Draw the reachability graph using graphviz. Shows ω for unbounded places.
    """
    dot = Digraph(name='ReachabilityGraph', format='png')
    dot.attr(rankdir='TB', imagescale='true', splines='true')# Top to Bottom, splines set to true for better edge routing
    dot.attr(graph='true', ranksep='1.2', nodesep='0.5') # increase rank and node separation
    dot.attr('node', shape='box')
    dot.attr('edge', color='darkblue')

    node_ids = {}
    depth_groups = {}
    for idx, (mkey, marking) in enumerate(nodes_map.items()):
        nid = f'M{idx}'
        node_ids[mkey] = nid
        # render math.inf as ω
        label_elems = []
        for v in marking:
            if v is math.inf:
                label_elems.append('ω') # represent infinity as ω
            else:
                label_elems.append(int(v))# convert to plain int for display
        label = str(label_elems)
        dot.node(nid, label=label)

        d = depths.get(mkey, 0)
        depth_groups.setdefault(d, []).append(nid)  # group nodes by depth for ranking because depths may not be continuous

    for d in sorted(depth_groups.keys()):
        with dot.subgraph() as s: # create subgraph for same rank
            s.attr(rank='same') # set same rank for nodes at same depth
            for nid in depth_groups[d]: # add nodes at this depth
                s.node(nid) # add node to subgraph

    for from_key, to_key, tindex in edges:
        # if keys were stored as tuples, map to ids
        if from_key in node_ids and to_key in node_ids:
            # use proper Graphviz edge label so transitions appear as t1..tn
            dot.edge(node_ids[from_key], node_ids[to_key], label=f't{tindex}')
    output_path = dot.render(filename, cleanup=True)

def main() -> None:

    pre = numpy.array([
        [1, 1, 0, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1]
    ])
    post = numpy.array([
        [0, 0, 0, 0, 0],
        [1, 0, 0, 1, 0],
        [0, 1, 0, 0, 0],
        [1, 0, 1, 0, 0],
        [0, 1, 0, 1, 1]
    ])
    startMarking = numpy.array([0, 1, 1, 0, 0])

    # Generate and display the full reachability graph first
    print("Generating reachability graph...")
    nodes, edges, depths = build_reachability_graph(pre, post, startMarking)

    # Analyze properties (Liveness, Boundedness, Reversibility)
    analyze_properties(nodes, edges, tuple(startMarking), pre.shape[1])

    filename = input('Enter filename: ')
    draw_reachability_graph(nodes, edges, depths, filename=filename)
    print(f"Reachability graph has been generated as {filename}.png")

if __name__ == "__main__":
    main()
