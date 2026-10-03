Stage 5F: fresh-session N=8 exact-H energy-transform solve

Why this stage exists
---------------------
Stage 5E Part A validated the exact-H energy-coordinate solver on N=4:
its FE error functionals matched the direct-LU reference very closely.
Stage 5E Part B did not fail mathematically; it was stopped by the memory
safety guard because only about 3.84 GB of physical memory was available,
while the symbolic N=8 Cholesky factor estimate was about 3.04 GB.

Run order
---------
1) Run RUN_STAGE5F_PREPARE_X0 once.
   It computes and saves the trusted N=8 lambda=1 free solution.

2) Close MATLAB completely.
   Close other memory-heavy applications if possible.

3) Start a fresh MATLAB session in this same folder.
   Run RUN_STAGE5F_N8_HIGHLAMBDA only.

The second script performs an early memory check.  If less than 7.5 GB is
available, it stops before assembly.  The exact-H solver also retains the
Stage-5E symbolic Cholesky guard, so it will not intentionally attempt an
unsafe factorization.

No finite-element bilinear form, boundary condition, projection, or
stabilization term is changed in Stage 5F.
