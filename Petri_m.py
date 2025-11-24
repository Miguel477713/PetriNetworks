import os, numpy, math
from graphviz import Digraph
from typing import Optional
from collections import deque, defaultdict

if os.name == 'nt':
    os.system('cls')  # Clear the console for Windows
elif os.name == 'posix':
    os.system('clear')  # Clear the console for Linux/Mac

def MarkingKey(marking):
    """Return a hashable key representing the start marking."""
    # Convert numpy array to tuple of ints
    if marking is not None:
        return tuple(int(x) for x in marking) # Convert numpy array to tuple of ints since numpy arrays are not hashable

class PetriNet:
    def __init__(self, pre, post, startMarking) -> None:

        if(len(pre) != len(post) or pre.shape[1] != post.shape[1]): # Check dimensions
            raise ValueError("Pre and Post sizes are different")

        self.pre = pre
        self.post = post

        self.adjacencyMatrix = self.post - self.pre # the adjacency matrix is simply post - pre

        self.startMarking = startMarking
        self.availableTransitions = self.GetAvailableTransistions()

    def MarkingKey(self):
        """Return a hashable key representing the start marking."""
        # Convert numpy array to tuple of ints in order to be hashable
        return tuple(int(x) for x in self.startMarking)

    def GetAvailableTransistions(self) -> numpy.ndarray:
        availableTransistions = numpy.empty(self.pre.shape[1]) # we use empty since we will fill all values

        for transition in range(self.pre.shape[1]): # for each transition in the net
            # Check if startMarking >= the specific column for this transition
            if numpy.all(self.startMarking >= self.pre[:, transition]): # if all places have enough tokens
                availableTransistions[transition] = 1 # transition is available
            else:
                availableTransistions[transition] = 0 # transition is not available
        return availableTransistions

    def GetNewMarking(self, firingTransition) -> Optional[numpy.ndarray]:
        if(firingTransition < 1 or firingTransition > self.pre.shape[1]): # Check valid transition number
            print(f"Transition t{firingTransition} is not defined")
            return None

        newMarking = numpy.empty(len(self.pre)) # prepare new marking array
        firingTransition -= 1 # adjust for 0-based index
        if self.availableTransitions[firingTransition] == 1: # if the transition is available
            newMarking = self.startMarking + self.adjacencyMatrix[:, firingTransition] # calculate new marking by adding the corresponding column of the adjacency matrix
        else:
            print(f"Transition t{firingTransition + 1} is not available")
            return None
        return newMarking

    def PrintAvailableTransitions(self) -> None:
        print(f"Available transitions from marking {self.startMarking} are: ")

        for i in range(len(self.availableTransitions)):
            if self.availableTransitions[i] == 1: # if transition is available
                print(f"t{i + 1}") # print transition number (1-based)

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

def analyze_properties(nodes: dict, edges: list, start_marking: tuple, num_transitions: int):
    """
    Analyze Boundedness, Liveness, and Reversibility based on the generated graph.
    """
    print("\nQualitative Analysis of the Petri Net:")
    # Helper: format marking elements to Python ints or 'ω'
    def format_marking(m):
        return tuple('ω' if x is math.inf else int(x) for x in m)

    #BOUNDEDNESS
    is_bounded = True
    for m in nodes.keys():
        if math.inf in m:
            is_bounded = False
            break

    print(f"Boundedness: {'BOUNDED' if is_bounded else 'UNBOUNDED (ω appears)'}")
    # Build adjacency list for graph traversal
    adj = defaultdict(list) # using defaultdict for convenience since we may have sparse connections
    for from_marking, to_marking, t_idx in edges:  # build adjacency list for liveness and reversibility checks
        adj[from_marking].append((to_marking, t_idx))  # append (neighbor, transition_index)
    # LIVENESS & DEADLOCKS
    deadlock_nodes = []
    for node in nodes:
        if not adj[node]:  # No outgoing edges
            # store formatted marking with plain ints (or 'ω')
            deadlock_nodes.append(format_marking(node))

    if deadlock_nodes:
        print(f"Deadlock: DETECTED at markings: {deadlock_nodes}")
        is_live = False  # If deadlock exists, it is not Live
    else:
        print("Deadlock: None (Deadlock-free)")

        # Check Strict Liveness:
        # A PN is live if from EVERY node, EVERY transition must be reachable.
        all_trans_live = True
        for t_idx in range(1, num_transitions + 1):
            # Check if transition t_idx is reachable from EVERY node
            for start_node in nodes:
                # BFS to find t_idx
                queue = [start_node]
                visited = {start_node}
                found_t = False
                while queue: # while there are nodes to explore
                    curr = queue.pop(0) # dequeue the first node
                    for neighbor, edge_t in adj[curr]: # for each neighbor and the transition that leads to it
                        if edge_t == t_idx: # if this edge corresponds to the transition we are looking for
                            found_t = True
                            break
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                    if found_t:
                        break

                if not found_t: # transition t_idx is not reachable from start_node
                    all_trans_live = False # this transition is not live
                    break # no need to check further for this transition
            if not all_trans_live:
                break
        is_live = all_trans_live # if all transitions are live

    print(f"Liveness: {'LIVE' if is_live else 'NOT LIVE'}")

    #REVERSIBILITY
    # A PN is reversible if and only if for all M in R(M0), M0 in R(M)
    # this means: Can we get back to start_marking from EVERY reachable node?
    is_reversible = True

    for node in nodes:
        if node == start_marking:
            continue

        # BFS to find start_marking
        queue = [node]
        visited = {node}
        found_start = False

        while queue:
            curr = queue.pop(0)
            if curr == start_marking:# found a path back to start_marking
                found_start = True
                break # no need to explore further

            for neighbor, _ in adj[curr]: # for each neighbor (ignore transition index)
                if neighbor not in visited:
                    visited.add(neighbor) # mark as visited
                    queue.append(neighbor) # enqueue neighbor

        if not found_start: # could not find path back to start_marking
            is_reversible = False # not reversible
            break # no need to check further

    print(f"Reversibility: {'REVERSIBLE' if is_reversible else 'NOT REVERSIBLE'}")
    print("-"*20)
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
    # output_path = dot.render(filename, cleanup=True)

def main() -> None:
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
    pre = numpy.array([
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 0, 1, 0],
    ])
    post = numpy.array([
        [0, 0, 0, 1, 0],
        [1, 0, 0, 0, 1],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
    ])
    currentMarking = startMarking = numpy.array([1, 0, 0, 0, 0])

    # Ejemplo pagina 66 libro UANL
    # pre = numpy.array([
    #     [1, 0, 0],
    #     [0, 1, 1],
    #     [0, 0, 0],
    #     [0, 1, 0]
    # ])
    # post = numpy.array([
    #     [0, 0, 1],
    #     [1, 1, 0],
    #     [0, 1, 0],
    #     [0, 0, 1]
    # ])
    # startMarking = numpy.array([1, 0, 0, 0])
    # Generate and display the full reachability graph first
    print("Generating reachability graph...")
    nodes, edges, depths = build_reachability_graph(pre, post, startMarking)

    # Analyze properties (Liveness, Boundedness, Reversibility)
    analyze_properties(nodes, edges, tuple(startMarking), pre.shape[1])

    filename = input('Enter filename: ')
    draw_reachability_graph(nodes, edges, depths, filename=filename)
    print(f"Reachability graph has been generated as {filename}.png")
    print("\nInteracting now with the Petri Net:")

    net = PetriNet(pre=pre, post=post, startMarking=startMarking)
    petriNetworks = {net.MarkingKey(): net}  # Initialize with start marking


    while(True):
        net.PrintAvailableTransitions()
        userInput = input(f"Enter a transition [1..{pre.shape[1]}] or enter * to stop: ")

        if userInput == '*':
            break

        # try to parse transition number
        try:
            number = int(userInput) # transition numbers start at 1
        except ValueError: # invalid input
            print("Invalid input, enter a number or * to stop.")
            continue # ask again

        currentMarking = net.GetNewMarking(number) # get new marking after firing transition
        if currentMarking is None:
            continue

        key = MarkingKey(currentMarking) # get marking key for lookup
        # check if we already have a PetriNet for this marking
        if key in petriNetworks:
            net = petriNetworks[key]
        else:
            # otherwise create and store a new PetriNet for this marking
            net = PetriNet(pre=pre, post=post, startMarking=currentMarking)
            petriNetworks[net.MarkingKey()] = net

if __name__ == "__main__":
    main()
