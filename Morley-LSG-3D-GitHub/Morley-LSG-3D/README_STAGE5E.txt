Stage 5E: exact-full-H energy-coordinate GMRES
===============================================
This stage does NOT change the Morley JS finite element formulation.

Why this stage exists
---------------------
Stage 5D showed that the exact full coupled symmetric-part preconditioner
reproduced the N=4 direct-LU FE error functionals to about 0.1%, even though
the coefficient-vector difference was larger than an unnecessarily strict
old threshold.  It also showed that a tiny ||Ax-b||/||b|| alone is not a
reliable accuracy certificate for the highly ill-conditioned lambda=1e5
system.

This stage uses a stronger coordinate transform.  After diagonal scaling,
H=(A+A')/2 is factored exactly as H=L*L'.  The full system is transformed to
B = L^{-1} A L^{-T}; its symmetric part is the identity in exact arithmetic.
GMRES is applied to this transformed system.

Run
---
RUN_STAGE5E_ENERGY_TRANSFORM

Part A validates the transformed solver on N=4 against direct LU.
Only if that passes does Part B attempt N=8 lambda=1e5.

Memory guard
------------
Before the N=8 exact Cholesky factorization, symbfact estimates fill.  If the
rough factor storage exceeds the configured guard, the script stops safely
before attempting Cholesky.  No FE result is then claimed.
