# Project instructions for Codex

## Project identity

- Communicate with the user in Chinese unless they ask otherwise.
- This is the existing Stage5F 3D annular pull-out MATLAB project, using the Morley JS projection scheme and a radial Aifantis reference.
- Read `README.md`, `README_STAGE5F.txt`, and `config/stage5_pullout_config.m` before numerical work.
- The mesh is generated in MATLAB. Do not assume Gmsh, SP scripts, `run_main.m`, or `run_all_tests.m` exist.
- Keep the current directory structure; `setup_stage5.m` configures the required paths.

## Numerical meaning

- Preserve the finite element formulation, P0 stress projection in J, P0 divergence projection in s, boundary treatment, error definitions, and physical parameters unless the task requires a specific change.
- Do not substitute homogeneous boundary data for the current reference-derived boundary traces.
- Distinguish solver accuracy from finite element error. A small residual alone is not a proof of accuracy or parameter robustness.
- The bundled `results/stage5a_verified_N4_baseline.csv` is required input inherited from the source package. Preserve it. Do not fabricate or silently update baseline values.
- Changes to geometry, material parameters, iota, or mesh invalidate reuse of old baselines and initial guesses unless compatibility is established.

## Execution

- Check that MATLAB is actually installed and licensed before claiming that code can execute. If it is unavailable, state that only static checks were performed.
- Do not change MATLAB code to Octave or Python merely to obtain a passing run.
- For an explicitly requested first environment check, start with `cfg=setup_stage5();`. For a requested small numerical diagnostic use `RUN_STAGE5C_N4_CALIBRATION`.
- `RUN_STAGE5E_ENERGY_TRANSFORM` also attempts N=8; do not treat it as an N=4-only check.
- The Stage5F sequence is `RUN_STAGE5F_PREPARE_X0`, termination of that MATLAB process, and then `RUN_STAGE5F_N8_HIGHLAMBDA` in a fresh process using the same project directory.
- Retain the generated `results/stage5f_N8_lambda1_x0.mat` between those processes. It is intentionally excluded from Git.
- Do not launch N=8 or increase to N=16 as an incidental check for documentation or packaging changes. Larger numerical runs must be within the user's requested scope and resource budget.
- Keep memory guards and acceptance criteria. Report a stopped or unaccepted solve honestly; do not relax limits merely to complete a run.
- Current available-memory queries run only on Windows. On Linux, verify the actual memory/container limit before costly solves; the symbolic factor estimate is not a total peak-memory guarantee.

## Files and results

- Keep generated matrices, solutions, caches, and bulk output out of version control. Preserve the baseline CSV exception in `.gitignore`.
- Do not commit credentials, MATLAB license files, or machine-specific authentication/configuration.
- Treat historical stage notes as provenance; current Stage5F instructions take precedence for the current entry points.
- Record the actual command, environment, parameters, and validation status when reporting a numerical run. Never present inherited CSV values as newly computed results.
- The original-source hash list records the packaging snapshot; if source is intentionally changed later, describe the change rather than concealing it by claiming the snapshot is still identical.
