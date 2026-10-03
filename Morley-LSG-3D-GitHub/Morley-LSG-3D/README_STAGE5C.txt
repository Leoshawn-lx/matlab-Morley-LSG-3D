Stage 5C -- lambda=1e5 solver-accuracy diagnostic

Purpose
-------
Stage 5B showed a clean lambda comparison in L2 and strain, but at N=8 the
high-derivative errors (grad-strain and three-bar norm) were about 1.6 times
the lambda=1 values, while ||div u_h|| also jumped relative to the N=4 direct
solve.  Because the lambda=1e5 right-hand side has norm O(1e10), a relative
residual O(1e-10) can still hide a non-negligible absolute residual.

This stage does NOT change the Morley JS finite element method.

Run
---
  RUN_STAGE5C_SOLVER_DIAGNOSTIC

Part A solves N=4, lambda=1e5 twice:
  (1) sparse direct LU -- reference linear algebra solution;
  (2) component-block GMRES plus residual-correction iterative refinement.
It compares the coefficient vectors and all error measures.

Part B solves N=8 for lambda=1 and 1e5 using the same block preconditioner,
but requires a final true residual <= 1e-12 and allows two residual-correction
steps.  It then reports the new lambda error ratios and N=4 -> N=8 rates.

Decision
--------
Do not start the iota sweep until this diagnostic has separated FE behavior
from Krylov-solver error at lambda=1e5.
