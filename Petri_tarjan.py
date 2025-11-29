import os, numpy, math
from graphviz import Digraph
from typing import Optional
from collections import deque, defaultdict

if os.name == 'nt':
    os.system('cls')
elif os.name == 'posix':
    os.system('clear')


def MarkingKey(marking):
    if marking is not None:
        return tuple(int(x) for x in marking)


class PetriNet:
    def __init__(self, pre, post, startMarking) -> None:
        if (len(pre) != len(post) or pre.shape[1] != post.shape[1]):
            raise ValueError("Pre and Post sizes are different")

        self.pre = pre
        self.post = post
        self.adjacencyMatrix = self.post - self.pre
        self.startMarking = startMarking
        self.availableTransitions = self.GetAvailableTransistions()

    def MarkingKey(self):
        return tuple(int(x) for x in self.startMarking)

    def GetAvailableTransistions(self) -> numpy.ndarray:
        availableTransistions = numpy.empty(self.pre.shape[1])
        for transition in range(self.pre.shape[1]):
            if numpy.all(self.startMarking >= self.pre[:, transition]):
                availableTransistions[transition] = 1
            else:
                availableTransistions[transition] = 0
        return availableTransistions

    def GetNewMarking(self, firingTransition) -> Optional[numpy.ndarray]:
        if (firingTransition < 1 or firingTransition > self.pre.shape[1]):
            print(f"Transition t{firingTransition} is not defined")
            return None

        newMarking = numpy.empty(len(self.pre))
        firingTransition -= 1
        if self.availableTransitions[firingTransition] == 1:
            newMarking = self.startMarking + self.adjacencyMatrix[:, firingTransition]
        else:
            print(f"Transition t{firingTransition + 1} is not available")
            return None
        return newMarking

    def PrintAvailableTransitions(self) -> None:
        print(f"Available transitions from marking {self.startMarking} are: ")
        for i in range(len(self.availableTransitions)):
            if self.availableTransitions[i] == 1:
                print(f"t{i + 1}")


def build_reachability_graph(pre: numpy.ndarray, post: numpy.ndarray, start_marking: numpy.ndarray):
    adjacency = post - pre
    num_transitions = pre.shape[1]
    start = tuple(int(x) for x in start_marking)

    nodes = {start: start}
    queue = deque([(start, [])])
    edges = []
    depths = {start: 0}

    def edge_exists(edgelist, a, b, t):
        return (a, b, t) in edgelist

    def add_markings(m, change):
        res = []
        for i in range(len(m)):
            if m[i] is math.inf:
                res.append(math.inf)
            else:
                res.append(m[i] + change[i])
        return tuple(res)

    while queue:
        marking, ancestors = queue.popleft()

        for t in range(num_transitions):
            enabled = True
            for i, req in enumerate(pre[:, t].astype(int)):
                mi = marking[i]
                if mi is math.inf: continue
                if mi < req:
                    enabled = False;
                    break
            if not enabled: continue

            change_vec = adjacency[:, t]
            new_m = add_markings(marking, change_vec)

            # Omega logic
            current_path_ancestors = ancestors + [marking]
            temp_m_list = list(new_m)
            for ancestor in current_path_ancestors:
                is_ge = True;
                strictly_greater = False
                for k in range(len(temp_m_list)):
                    val_new = temp_m_list[k];
                    val_anc = ancestor[k]
                    if val_anc is math.inf and val_new is not math.inf: is_ge = False; break
                    if val_new is not math.inf and val_anc is not math.inf and val_new < val_anc: is_ge = False; break
                    if val_new != val_anc: strictly_greater = True

                if is_ge and strictly_greater:
                    for k in range(len(temp_m_list)):
                        if temp_m_list[k] is not math.inf and ancestor[k] is not math.inf:
                            if temp_m_list[k] > ancestor[k]: temp_m_list[k] = math.inf

            new_m = tuple(temp_m_list)

            if not edge_exists(edges, marking, new_m, t + 1):
                edges.append((marking, new_m, t + 1))

            if new_m not in nodes:
                nodes[new_m] = new_m
                depths[new_m] = depths[marking] + 1
                queue.append((new_m, current_path_ancestors))

    return nodes, edges, depths


def tarjan_scc(adj_map: dict):
    index = {};
    lowlink = {};
    onstack = set();
    stack = [];
    sccs = [];
    next_index = 0

    def strongconnect(v):
        nonlocal next_index
        index[v] = next_index;
        lowlink[v] = next_index;
        next_index += 1
        stack.append(v);
        onstack.add(v)

        for w in adj_map.get(v, []):
            if w not in index:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in onstack:
                lowlink[v] = min(lowlink[v], index[w])

        if lowlink[v] == index[v]:
            scc = []
            while True:
                w = stack.pop();
                onstack.remove(w);
                scc.append(w)
                if w == v: break
            sccs.append(scc)

    for v in adj_map.keys():
        if v not in index: strongconnect(v)
    return sccs


def find_elementary_paths(adj, start_node, scc_set):
    """Finds paths but returns a dictionary grouping them by the transitions used."""
    paths = []
    LIMIT = 50

    def dfs(current, current_path, visited):
        if len(paths) >= LIMIT: return

        if current == start_node and len(current_path) > 0:
            paths.append(list(current_path))
            return

        neighbors = sorted(adj[current], key=lambda x: x[1])
        for next_node, t_idx in neighbors:
            if next_node in scc_set:
                if next_node not in visited or next_node == start_node:
                    dfs(next_node, current_path + [(current, t_idx, next_node)], visited | {next_node})

    dfs(start_node, [], {start_node})

    # Logic: Group by signature
    grouped_paths = {}
    for p in paths:
        t_indices = tuple(sorted([step[1] for step in p]))
        if t_indices not in grouped_paths:
            grouped_paths[t_indices] = []
        grouped_paths[t_indices].append(p)

    return grouped_paths


def analyze_properties(nodes: dict, edges: list, start_marking: tuple, num_transitions: int):
    print("\n" + "=" * 60)
    print("QUALITATIVE ANALYSIS (Based on course notes: 11-Redes de Petri)")
    print("=" * 60)

    def pretty_marking(m):
        parts = []
        for i, val in enumerate(m):
            place_name = f"P{i + 1}"
            if val is math.inf:
                parts.append(f"{place_name}(ω)")
            elif val > 0:
                if val == 1:
                    parts.append(place_name)
                else:
                    parts.append(f"{int(val)}{place_name}")
        return "{" + ", ".join(parts) + "}" if parts else "{}"

    # 1. BOUNDEDNESS
    is_bounded = True
    omega_node = None
    for m in nodes.keys():
        if math.inf in m:
            is_bounded = False;
            omega_node = pretty_marking(m);
            break

    print(f"\n1. BOUNDEDNESS: {'BOUNDED' if is_bounded else 'UNBOUNDED'}")
    if is_bounded:
        print(f"   Reason: No place count exceeded a finite limit k.")
        print(f"   The symbol ω (infinity) does not appear in the generated graph.")
    else:
        print(f"   Reason: UNBOUNDED because ω appears in marking {omega_node}.")
        print(f"   Reference: 'If ω appears -> PN unbounded' (Slide 305).")

    # Build Adjacency
    adj = defaultdict(list)
    for from_m, to_m, t_idx in edges:
        adj[from_m].append((to_m, t_idx))
    neighbor_map = {n: [nbr for nbr, _ in adj[n]] for n in adj}
    for n in nodes: neighbor_map.setdefault(n, [])

    # 1.5 CYCLIC BEHAVIOR
    sccs = tarjan_scc(neighbor_map)
    cyclic_sccs = [scc for scc in sccs if len(scc) > 1 or (len(scc) == 1 and scc[0] in neighbor_map.get(scc[0], []))]

    print(f"\n1.5 STRONGLY CONNECTED COMPONENTS (Tarjan) - total: {len(sccs)}")

    if not cyclic_sccs:
        print("   No cycles detected.")
    else:
        for i, scc in enumerate(cyclic_sccs, start=1):
            print(f"   SCC #{i} (Cyclic, Contains {len(scc)} states):")

            root = start_marking if start_marking in scc else sorted(scc)[0]

            grouped_cycles = find_elementary_paths(adj, root, set(scc))

            if not grouped_cycles:
                print("      (Complex loop structure, displaying basic connections)")
                for u in scc:
                    for v, t in adj[u]:
                        if v in scc: print(f"      {pretty_marking(u)} -- t{t} --> {pretty_marking(v)}")
            else:
                for t_sig, path_list in grouped_cycles.items():
                    print(f"      LOGICAL CYCLE DETECTED:")
                    t_str = ", ".join([f"t{t}" for t in t_sig])
                    print(f"        Transitions used: {{{t_str}}}")
                    if len(path_list) > 1:
                        print(f"        (Merged {len(path_list)} redundant sequences caused by parallel execution)")

                    print(f"        Representative Sequence:")
                    print(f"          Start: {pretty_marking(path_list[0][0][0])}")
                    for (u, t, v) in path_list[0]:
                        print(f"             |-- t{t} --> {pretty_marking(v)}")
                    print("          (Return to Start)\n")

    # 2. LIVENESS
    deadlock_nodes = [pretty_marking(n) for n in nodes if not adj[n]]
    print(f"\n2. LIVENESS: ", end="")
    if deadlock_nodes:
        print("NOT LIVE (Deadlock Detected)")
        print(f"   Reason: The system contains deadlocks at marking(s): {', '.join(deadlock_nodes)}.")
        print(f"   Reference: 'Deadlocks: for a given marking, any transition cannot be fired' (Slide 260).")
    else:
        all_live = True
        failed_trans = -1
        failed_node = None
        for t_idx in range(1, num_transitions + 1):
            for start_node in nodes:
                queue = [start_node];
                visited = {start_node};
                found = False
                while queue:
                    curr = queue.pop(0)
                    for nbr, edge_t in adj[curr]:
                        if edge_t == t_idx: found = True; break
                        if nbr not in visited: visited.add(nbr); queue.append(nbr)
                    if found: break
                if not found: all_live = False; failed_trans = t_idx; failed_node = pretty_marking(start_node); break
            if not all_live: break

        if all_live:
            print("LIVE")
            print("   Reason: The net is Deadlock-free AND every transition is live.")
            print("   (Every transition is reachable from every reachable marking).")
            print("   Reference: 'A PN is live... iff t is live for all t' (pp.6/9).")
        else:
            print("NOT LIVE (Partial Functioning)")
            print(f"   Reason: Transition t{failed_trans} is not live.")
            print(f"   Once the system reaches {failed_node}, t{failed_trans} can never fire again.")

    # 3. REVERSIBILITY
    reversible = True
    irreversible_node = None
    for node in nodes:
        if node == start_marking: continue
        q = [node];
        vis = {node};
        found = False
        while q:
            c = q.pop(0)
            if c == start_marking: found = True; break
            for nbr, _ in adj[c]:
                if nbr not in vis: vis.add(nbr); q.append(nbr)
        if not found: reversible = False; irreversible_node = pretty_marking(node); break

    print(f"\n3. REVERSIBILITY: {'REVERSIBLE' if reversible else 'NOT REVERSIBLE'}")
    if reversible:
        print(f"   Reason: It is possible to return to the initial marking M_0 from ALL reachable states.")
        print(f"   Reference: 'A PN is reversible... iff for all M, M_0 is in R(NS, M)' (Slide 7/9).")
    else:
        print(f"   Reason: It is NOT possible to return to M0 from marking {irreversible_node}.")
        print(
            f"   Reference: 'Related with the existence of sequences leading the PN to the initial marking' (Slide 285).")
    print("=" * 60 + "\n")


# --- RESTORED ORIGINAL GRAPHVIZ FUNCTION ---
def draw_reachability_graph(nodes_map, edges, depths, filename='reachability_graph'):
    """
    Draw the reachability graph using graphviz. Shows ω for unbounded places.
    """
    dot = Digraph(name='ReachabilityGraph', format='png')
    dot.attr(rankdir='TB', imagescale='true',
             splines='true')  # Top to Bottom, splines set to true for better edge routing
    dot.attr(graph='true', ranksep='1.2', nodesep='0.5')  # increase rank and node separation
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
                label_elems.append('ω')  # represent infinity as ω
            else:
                label_elems.append(int(v))  # convert to plain int for display
        label = str(label_elems)
        dot.node(nid, label=label)

        d = depths.get(mkey, 0)
        depth_groups.setdefault(d, []).append(
            nid)  # group nodes by depth for ranking because depths may not be continuous

    for d in sorted(depth_groups.keys()):
        with dot.subgraph() as s:  # create subgraph for same rank
            s.attr(rank='same')  # set same rank for nodes at same depth
            for nid in depth_groups[d]:  # add nodes at this depth
                s.node(nid)  # add node to subgraph

    for from_key, to_key, tindex in edges:
        # if keys were stored as tuples, map to ids
        if from_key in node_ids and to_key in node_ids:
            # use proper Graphviz edge label so transitions appear as t1..tn
            dot.edge(node_ids[from_key], node_ids[to_key], label=f't{tindex}')

    try:
        output_path = dot.render(filename, cleanup=True)
        print(f"Graph generated: {output_path}")
    except:
        print("GraphViz error: Ensure GraphViz is installed and in your PATH.")


def main() -> None:
    # --- MATRICES FROM FIGURA 6.8 (Pre-configured) ---
    print("Using Matrix definition from Figura 6.8...")

    # pre = numpy.array([
    #     [1, 0, 0, 0, 0],  # P1 -> t1
    #     [0, 1, 0, 0, 0],  # P2 -> t2
    #     [0, 0, 1, 0, 0],  # P3 -> t3
    #     [1, 0, 0, 0, 0],  # P4 -> t1
    #     [0, 0, 0, 1, 0],  # P5 -> t4
    #     [0, 0, 0, 0, 1]  # P6 -> t5
    # ])
    #
    # post = numpy.array([
    #     [0, 0, 1, 0, 0],  # t3 -> P1
    #     [1, 0, 0, 0, 0],  # t1 -> P2
    #     [0, 1, 0, 0, 0],  # t2 -> P3
    #     [0, 0, 0, 0, 1],  # t5 -> P4
    #     [1, 0, 0, 0, 0],  # t1 -> P5
    #     [0, 0, 0, 1, 0]  # t4 -> P6
    # ])

    startMarking = numpy.array([1, 0, 0, 1, 0, 0])

    # APUNTES UANL pp. 65
    # pre = numpy.array([
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 1, 0, 0],
    #     [0, 0, 0, 1, 1],
    #     [0, 0, 0, 0, 1],
    # ])
    # post = numpy.array([
    #     [0, 0, 0, 0, 1],
    #     [1, 0, 0, 1, 0],
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 1, 0, 0],
    # ])
    # currentMarking = startMarking = numpy.array([1, 0, 0, 0, 0])

    # Ejemplo de diapositivas pp. 7, red no acotada:
    # pre = numpy.array(
    #     [[1, 1, 0],
    #     [0, 0, 1],
    #     [0, 0, 1],
    #     ])
    # post = numpy.array(
    #     [[1, 0, 0],
    #      [1, 1, 0],
    #      [0, 1, 1],
    #      ])
    # currentMarking = startMarking = numpy.array([1, 0, 0])

    # Apuntes Diapositivas
    # pre = numpy.array([[1, 0, 0, 0, 0],
    #                 [0, 1, 0, 0, 0],
    #                 [0, 0, 1, 0, 0],
    #                 [1, 0, 0, 0, 0],
    #                 [0, 0, 0, 1, 0],
    #                 [0, 0, 0, 0, 1],
    #                 ])
    #
    # post = numpy.array([[0, 0, 1, 0, 0],
    #                     [1, 0, 0, 0, 0],
    #                     [0, 1, 0, 0, 0],
    #                     [0, 0, 0, 0, 1],
    #                     [1, 0, 0, 0, 0],
    #                     [0, 0, 0, 1, 0],
    #                 ])
    # currentMarking = startMarking = numpy.array([1, 0, 0, 1, 0, 0])

    # Ejemplo presentación viernes
    # pre = numpy.array([
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 1, 0, 0],
    #     [0, 0, 0, 1, 1],
    #     [0, 0, 0, 1, 0],
    # ])
    # post = numpy.array([
    #     [0, 0, 0, 1, 0],
    #     [1, 0, 0, 0, 1],
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 1, 0, 0],
    # ])
    # currentMarking = startMarking = numpy.array([1, 0, 0, 0, 0])

    # Ejemplo pagina 66 libro UANL
    pre = numpy.array([
        [1, 0, 0],
        [0, 1, 1],
        [0, 0, 0],
        [0, 1, 0]
    ])
    post = numpy.array([
        [0, 0, 1],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    startMarking = numpy.array([1, 0, 0, 0])

    # Ejemplo pagina 66 libro UANL
    # pre = numpy.array([
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 1, 0, 0],
    #     [1, 0, 0, 0, 0],
    #     [0, 0, 0, 1, 0],
    #     [0, 0, 0, 0, 1]
    # ])
    # post = numpy.array([
    #     [0, 0, 1, 0, 0],
    #     [1, 0, 0, 0, 0],
    #     [0, 1, 0, 0, 0],
    #     [0, 0, 0, 0, 1],
    #     [1, 0, 0, 0, 0],
    #     [0, 0, 0, 1, 0]
    # ])
    # startMarking = numpy.array([1, 0, 0, 1, 0, 0])

    print("Generating reachability graph...")
    nodes, edges, depths = build_reachability_graph(pre, post, startMarking)

    analyze_properties(nodes, edges, tuple(startMarking), pre.shape[1])

    filename = input('Enter filename for graph image (default: petri_graph): ')
    if not filename: filename = "petri_graph"
    draw_reachability_graph(nodes, edges, depths, filename=filename)

    print("\nInteracting now with the Petri Net (Manual Mode):")
    net = PetriNet(pre=pre, post=post, startMarking=startMarking)
    petriNetworks = {net.MarkingKey(): net}

    while (True):
        net.PrintAvailableTransitions()
        userInput = input(f"Enter a transition [1..{pre.shape[1]}] or enter * to stop: ")

        if userInput == '*': break

        try:
            number = int(userInput)
        except ValueError:
            print("Invalid input.")
            continue

        currentMarking = net.GetNewMarking(number)
        if currentMarking is None: continue

        key = MarkingKey(currentMarking)
        if key in petriNetworks:
            net = petriNetworks[key]
        else:
            net = PetriNet(pre=pre, post=post, startMarking=currentMarking)
            petriNetworks[net.MarkingKey()] = net


if __name__ == "__main__":
    main()