Stage 5D: N=4 full-coupled preconditioner calibration

Run:
  RUN_STAGE5D_N4_FULLH_CALIBRATION

Purpose:
The previous component-block preconditioner reproduced L2/profile quantities but polluted
high-order quantities for lambda=1e5.  This diagnostic keeps the FE formulation unchanged
and tests whether the missing cross-component lambda coupling in the preconditioner is the
cause.

Two iterative candidates are compared with the same N=4 direct-LU reference:
  A. GMRES + exact sparse Cholesky of the full symmetric part H=(As+As')/2.
  B. GMRES + incomplete Cholesky of the SAME full coupled H.

If candidate B matches direct LU, it is suitable to try on N=8.  If only A matches, the
finite element formulation is fine but the deployable incomplete preconditioner still needs
tuning.
