# Petri Net Reachability Information

This project implements a Petri Net analyzer using the **Karp-Miller algorithm**, which generates a reachability graph and performs qualitative analysis for Boundedness, Liveness, and Reversibility.

---

## Generalities of the Karp-Miller Algorithm

The **Karp-Miller algorithm** constructs the **Coverability Graph** for a Petri Net. Its core purpose is to transform a potentially infinite state space (found in nets that can accumulate tokens endlessly) into a finite, analyzable graph representation by introducing the **omega ($\omega$) symbol**.
Reference: Understanding Petri Nets by Wolfgang Reisig,2013, chapter 14

The algorithm operates as a **Breadth-First Search (BFS)** across the markings. It systematically explores all reachable states by firing every enabled transition from a given state. The process halts when all possible transitions have been fired, and no new, unique markings (including those containing $\omega$) can be generated.

### The $\omega$ Abstraction Rule

This rule ensures the process ends by identifying growing components. 

* When a new state is generated, it is compared against **all ancestors** on the path from the starting state.
* If the new state is found to be greater than or equal to an ancestor in all token positions, and **strictly greater in at least one position**, the rule is triggered.
* The positions showing strict growth are replaced by $\omega$, which signifies the place is **unbounded** and can accumulate an arbitrarily large number of tokens.

---

## Reachability & Coverability Details

### Initialization
The algorithm begins at the initial state.
* **Example:** `[1, 0, 0, 1, 0, 0]`

### The Enabling Check
A transition is allowed to fire if and only if the current state has enough tokens for every input required by that transition.

**Example at State `[1, 0, 0, 1, 0, 0]`:**
* **Transition $t_1$ (Requires 1 from position $p_1$, 1 from $p_4$):**
    * $M(p_1) \ge 1$ and $M(p_4) \ge 1$.
    * **Result:** Enabled.
* **Transition $t_2$ (Requires 1 from position $p_2$):**
    * $M(p_2) = 0$, which is not $\ge 1$.
    * **Result:** Disabled.

### Firing and State Evolution
If a transition is enabled, the tokens are removed from the input places and added to the output places, generating the next state.
* **Action:** Fire $t_1$.
* **Result:** `[0, 1, 0, 0, 1, 0]`

### Numerical Example (The $\omega$ Rule)
The $\omega$ appears when the system finds a recurring growth pattern relative to an older state (a "Grandparent" scenario).

1.  **Ancestor (Root):** `[1, 0, 0, 0]`
2.  **New Raw Calculation:** If firing transitions leads to a state that would be `[1, 0, 0, 1]`.
3.  **Comparison:** `[1, 0, 0, 1]` is greater than or equal to the Root `[1, 0, 0, 0]` in all positions, and strictly greater in the last position ($1 > 0$).
4.  **Final Node:** The algorithm replaces the growing value: `[1, 0, 0, ω]`.

---

## Qualitative Analysis

Properties are deduced directly from the structure of the generated graph using the following logic:

* **$\omega$ appears $\rightarrow$ The net is unbounded**
* **The graph is strongly connected $\rightarrow$ The net is reversible**
* **Every strongly connected component contains all the transitions $\rightarrow$ The net is live**
* **Exists a terminal node $\rightarrow$ Deadlock**

### Boundedness
* **Definition:** A net is bounded if the number of tokens in every place never exceeds a fixed, finite number.
* **Analysis:** If the symbol $\omega$ appears in any state, the net is declared **Unbounded**.

### Strongly Connected Components (SCC) and Reversibility

* **What is an SCC?**
    An **SCC** is a subset of states in a directed graph where **every state is reachable from every other state** within that subset. This represents a perfect, closed loop.

* **Analysis Logic for Reversibility:**
    A net is **Reversible** if the initial state ($M_0$) is reachable from every other state. This requires the **entire reachability graph** to form a **single Strongly Connected Component**.

### Liveness
Liveness is a stronger property related to the permanent capability of transitions.

* **1. Deadlock (No Liveness)**
    * **Definition:** A state where **no** transitions are enabled (a "terminal node" in the graph).

* **2. Strict Liveness**
    * **Definition:** The net is strictly Live if it is always possible to eventually fire **any** transition in the net, no matter what state the system is in.
    * **The Livelock Problem:** A net can be deadlock-free but still Not Live if it enters a **Livelock** scenario (a cycle where some transitions are permanently excluded).
    * **Detection:** The code verifies that **"Every strongly connected component contains all the transitions."** This ensures that any closed loop the system settles into still permits all defined functions.

---