## Una RP es viva ssi en el grafo de cobertura, todas las componentes fuertemente conexas contienen todas las transiciones.

# Petri Net Reachability Information

This project implements a Petri Net analyzer using the **Karp-Miller algorithm**, which generates a reachability (coverability) graph and performs qualitative analysis for Boundedness, Liveness, Reversibility, and computes Strongly Connected Components (SCCs) using Tarjan's algorithm.

---

## What's new (changelog)

- Added Tarjan's algorithm to compute Strongly Connected Components (SCCs) in the reachability/coverability graph.
- The script `Petri_m.py` now prints SCCs found and includes a clearer, refactored analysis flow.
- Graph drawing with Graphviz remains available (output PNG) and nodes containing ω are rendered as 'ω'.

---

## Generalities of the Karp-Miller Algorithm

The **Karp-Miller algorithm** constructs the **Coverability Graph** for a Petri Net. Its core purpose is to transform a potentially infinite state space (found in nets that can accumulate tokens endlessly) into a finite, analyzable graph representation by introducing the **omega (ω) symbol**.

The algorithm operates as a **Breadth-First Search (BFS)** across the markings. It systematically explores all reachable states by firing every enabled transition from a given state. The process halts when all possible transitions have been fired, and no new, unique markings (including those containing ω) can be generated.

### The ω Abstraction Rule

This rule ensures the process ends by identifying growing components.

- When a new state is generated, it is compared against **all ancestors** on the path from the starting state.
- If the new state is found to be greater than or equal to an ancestor in all token positions, and **strictly greater in at least one position**, the rule is triggered.
- The positions showing strict growth are replaced by ω, which signifies the place is **unbounded** (can accumulate an arbitrarily large number of tokens).

---

## Reachability & Coverability Details (example)

### Initialization
The algorithm begins at the initial state. Example: `[1, 0, 0, 1, 0, 0]`.

### Enabling check
A transition is allowed to fire if and only if the current state has enough tokens for every input required by that transition.

### Firing and State Evolution
If a transition is enabled, tokens are removed from input places and added to output places, generating the next state.

### The ω Rule Example
The ω appears when the system finds a recurring growth pattern relative to an older state (a "Grandparent" scenario). See the code comments and the Karp-Miller description above.

---

## Qualitative Analysis Performed by `Petri_m.py`

The script prints the following information after generating the reachability / coverability graph:

1. BOUNDEDNESS: Whether the net is BOUNDED or UNBOUNDED (ω presence means UNBOUNDED).
2. STRONGLY CONNECTED COMPONENTS (SCCs) computed with Tarjan's algorithm. Each SCC is listed as a set of markings (ω rendered as 'ω').
3. LIVENESS: Reports deadlocks and attempts to determine strict liveness (whether every transition is live from every marking).
4. REVERSIBILITY: Whether the initial marking M0 is reachable from every reachable marking (graph strongly connected to M0).

### SCCs (Tarjan)

- SCCs are useful to reason about recurrent behavior, livelocks, and closed sets of reachable states.
- The new output prints all SCCs found and their members (markings). This helps locate cycles or closed components where certain transitions may never fire.

---

## How to run

1. Install dependencies (Graphviz required and Python packages):

```powershell
pip install -r requirements.txt
```

2. Run the script from the project root:

```powershell
python Petri_m.py
```

3. The program prints the reachability analysis and asks for a filename to render the graph as PNG.

---

## Common questions / Troubleshooting

- Q: The PNG is not being generated, why?
  - A: Make sure `graphviz` (the system package) is installed and `dot` is on your PATH. The Python `graphviz` package is only a wrapper: you still need the Graphviz binaries. On Windows, add Graphviz `bin` directory to PATH and reopen the terminal/IDE.

- Q: I see `np.int64` in outputs (like deadlock nodes), how to avoid it?
  - A: The code was refactored to convert numeric marking elements to plain Python `int` for printing and to render ω as a string. This prevents `np.int64(...)` in printed output.

- Q: Why are there no transitions labeled `t3` or others in the generated graph PNG?
  - A: The script only draws edges for transitions that are actually fired between markings while building the graph. If a transition is never enabled in any reachable marking (or its edges were filtered), it won't appear in the diagram. Check the printed list of edges or analyze the available transitions at specific markings.

- Q: How to check why a line was commented out in Git?
  - A: Use `git blame <file>` or inspect the history with `git log -p -- <file>` to see when and why changes were made. (Run these commands in your shell/IDE.)

---

## Suggested small next steps (if you want me to continue)

- Add unit tests for the Karp-Miller and Tarjan parts (simple nets to validate SCC counts and ω appearance).
- Make the graph drawing optional (CLI flag) and write the output path to the console.
- Add more descriptive node labels (index + marking) so the PNG is easier to cross-reference with printed SCCs.

---

## References

- Notes and slides from your course (referenced in code comments).
- Reisig, W. (2013). Understanding Petri Nets.


---

(End of README)
