Stage 5B: N=8 lambda-robustness continuation
============================================

Run:
    RUN_STAGE5_LAMBDA_N8

Problem:
    mu = 1
    iota = 0.1
    lambda = 1 and 1e5
    N = 8 pure-MATLAB annular mesh

The corrected nonhomogeneous boundary RHS J(v,g) is retained.
The JS scheme is unchanged: P0 stress in J and P0 divergence in s.

Why lambda=1 is run again:
    Stage 5B compares both lambda values in the same code path and stores
    the two N=8 radial profiles directly. The converged lambda=1 free
    solution is then used only as the GMRES initial guess for lambda=1e5.

The file results/stage5a_verified_N4_baseline.csv contains the already
validated N=4 errors from Stage 5A. It is used only to print the N=4 -> N=8
observed rates for each lambda.

Acceptance:
    The finite element case is accepted only if the true residual of the
    original unscaled system satisfies ||A*x-b||/||b|| <= 1e-8.
