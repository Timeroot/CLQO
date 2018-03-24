# CLQO
A Solver for Unconstrained Binary Quadratic Optimization (UBQO, BQO, QUBO) and Max 2-SAT, based on semidefinite relaxation with constraint learning.

# The problem
Binary Quadratic Optimization is maximizing a quadratic polynomial over variables each constrained to be +1 or -1. (Equivalently, they can be constrained to be 0 or 1, or any pair of values, given by a simple linear substitution.) In general, this problem is NP hard. However, it's also a problem arises in practice not infrequently, so research is largely focused in finding good approximate solutions. For some immediate reductions, it is very efficient at describing Max Independent Set, Max 2-SAT, Max Cut, and Max Clique.

# The solutions
Existing approaches mostly fall into three categories:
* Branch-and-bound algorithms. These try to explore one branch of the tree and get an upper- and lower-bound for that branch; then they explore the other branch. Recurse, exploring more branches, allowing tighter and tighter bounds. When the upper bound for one branch is below the lower bound of another, then the former can be discarded. Finding lower bounds is akin to finding good assignments; finding upper bounds is harder.
* Stochastic methods. These include global and local search methods. Global search methods include things like genetic algorithms, or just large sets of trial vectors with different distributions. Local search methods include simulated annealing, or the 'tabu' search, to tweak small clusters of variables.
* Semidefinite relaxation. While the objective function is quadratic in the underlying variables, one can add in quadratically many dummy variables for the objective which cause it to be linear, thus, easily optimized. The problem then is enforcing that the quadratic variables match up with the underlying linear variables.

This program focuses on a variant of the third.

# How the relaxation works
To enforce the correct relation between quadratic variables, semidefinite relaxation checks that the matrix of assigned products if ''positive semi-definite'', a necessary constraint. This is still easy to optimize with, because the constraint defines a convex set. Thus traditional linear optimization algorithms like the ellipsoid algorithm can quickly converge to the optimum.

However, this constraint is not quite sufficient. For example -- for 3 variables, the allowed space looks something like Reuleaux tetrahedron (although with straight "edges"), although the actually valid assignments correpsond to the corners of the inscribing tetrahedron. Thus, this is not a perfect description of the problem: it is a relaxation, with a "fake" larger maximum value. Although on many problems the linear objective will "slide" on the "puffed out" space to the corners, in some instances it will not. The approach is still attractive because it can typically be 'rounded' at the end to a valid solution, one which can be proven to a give a good approximation ratio to the optimum; additionally, the value of the un-rounded assignment is the global optimum of the relaxation gives an upper bound to the problem, so that you get good lower and upper bounds. This can naturally be coupled with Branch/Bound algorithms.

There are many good software packages that implement Semidefinite Programming (SDP) to solve this problem, such as CSDP.

# What CLQO does

In the Reuleaux triangle story, it should be apparent that adding in linear constraints corresponding to the sides of the tetrahedron will perfectly define the actual convex hull of the allowed assignments. Optimizing a linear functional will then necessarily give the global optimum. By standard linear programming algorithms, adding in linear constraints is fine. The problem is that with large numbers of variables, there become exponentially many such linear facets to add in. There are O(n^3) such constraints on 'local' sets of 3 variables, then O(n^5) on sets of 5, then O(n^6), and so on. While these linear constraints enforce local behavior, the semidefiniteness property is supposed to enforce global behavior. Althoug the first set of n^3 constraints help fairly well, this is still not very practical, and higher-orders do not even improve the appoximation ratio, unfortunately. (That means: they don't help much.) And there are O(n^n) constraints total, horrible!

So: we solve the problem _without_ the semidefiniteness constraint enforced. The resulting assignment violates this somehow; the idea is that now we can check _how_ we violate the semidefiniteness to find a relevant constraint. In general this is awfully hard -- you will not always find a constraint easily -- but it does allow you to quickly identify, "ah, there is a particular problem with this set of 3 variables, I will add in that constraint". And so on. In general, the space being searched is only n^2 dimensional, so the true optimum is determined by some set of only n^2 constraints. Again -- finding them all is not always possible. But finding many quickly is.

Thus, CLQO: the Constraint Learning Quadratic Optimizer. As far as this author is aware, this is unqiue in its approach.

The eventual goal is that in addition to learning these linear constraints, it will apply them with the semidefiniteness constraint for the final round-off to get a good solution. Currently the priorities are:
 * Implementing detection for a broad variety of constraint types (such as detailed at http://comopt.ifi.uni-heidelberg.de/software/SMAPO/cut/cut.html )
 * Intelligently recognizing when a previously-added constraint that is now slack is unlikely to be used any more, and removing it. The objective function still monotonically improves, so this does not harm the completeness of the solver, but can greatly improve speed. CLQO currently tries to do this, but further tuning is likely merited.
 * Calling a SDP solver periodically to check if the current constraints suffice, together with the SDP constraint, for a global optimum.
 * Detecting large batches of constraints in one call. Currently constraint detection takes far less time than the linear optimization, and could be batched nicely.
 * Intelligently switching between simplex and interior-point methods for optimization, depending on which is faster at that point. (For reasonable size problems, simplex seems to consistently outperform for the whole process.)
 * Trying to solve with an incremental linear solver, as opposed to a "fresh" solution each time. Might be irrelevant if clause turnover rates get high enough.
 * Connecting to a branch-and-bound solver for CLQO to be used as a branch-evaluation subroutine. Alternately, implementing brancn-and-bound within CLQO as appropriate: when a problem is not converging well, trying to branch once or twice and run programs on each separately.

CLQO exists and is usable, but is not easily so. Hopefully development will continue.
