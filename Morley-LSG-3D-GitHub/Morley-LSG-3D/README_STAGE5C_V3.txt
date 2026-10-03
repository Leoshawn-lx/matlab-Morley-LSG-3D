Stage 5C-v3 fair solver calibration
===================================
Run:
    RUN_STAGE5C_N4_CALIBRATION

Why this version exists
-----------------------
The previous Stage 5C attempted an unnecessarily strict 1e-12 residual from a zero
initial guess at lambda=1e5.  It stagnated and therefore aborted before the useful
quantity -- the actual GMRES-vs-direct finite-element comparison -- was printed.

This version performs a fair calibration of the actual Stage-5B strategy:
1. Solve N=4, lambda=1 by direct LU and keep that free vector as x0.
2. Solve the SAME N=4, lambda=1e5 system by direct LU (reference).
3. Solve N=4, lambda=1e5 with the Stage-5B component-block GMRES policy,
   using the lambda=1 solution as x0.
4. Continue the diagnostic comparison even if GMRES misses the residual target.
5. Print vector differences, error ratios, div ratio and radial-profile difference.

No finite-element bilinear form, boundary condition, projection, stabilization, mesh,
or error definition has been changed.

Do NOT run the N=8 solver again until this N=4 calibration has been interpreted.
