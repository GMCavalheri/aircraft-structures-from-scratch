# CLAUDE.md

From-scratch aircraft structures toolkit (stress/strain → beams → FEA → buckling → composite
laminates → fatigue).
Plan: `AERO-02-aircraft-structures-from-scratch-PLAN.md` (tick the phase checkboxes as phases land).

## Commands
- `uv sync` — install deps (numpy, scipy, matplotlib; dev: pytest, ruff, nbclient)
- `uv run pytest` — tests
- `uv run ruff check . && uv run ruff format --check .` — lint/format (CI runs both)
- `uv run python validation/run_all.py` — regenerate every validation figure (`docs/figures/`)
  and results table (`validation/results/`); single scripts also run on their own
- `uv run python notebooks/build_notebooks.py` — regenerate and execute the notebooks (edit cell
  sources there, never the .ipynb JSON); rebuilds are byte-identical

## Layout
- `src/structures/` — single package; subpackages mirror the plan: `materials`, `stress_strain`,
  `beam_theory`, `fea_solver`, `buckling`, `composite_laminates`, `fatigue`; `plotting.py` shared.
- `tests/` — pytest, one file per module. `validation/` — scripts + cited reference data.
- `docs/` — theory note per module, validation report, figures. `notebooks/` — executed demos.

## Conventions
- **SI units internally** (N, m, Pa). Helpers that take or return other units say so in the
  name (`*_ksi`, `*_mm`, `*_deg`). 1 ksi = 6.894757 MPa.
- Angles are **radians** internally; public helpers that take degrees say so (`*_deg`). Ply
  angles are measured from the laminate x-axis, counter-clockwise positive.
- Tensile stress is positive. Stress vectors in Voigt order: 2D `[sx, sy, txy]`,
  3D `[sx, sy, sz, tyz, txz, txy]`; engineering shear strain (gamma = 2 eps) in strain vectors.
- Beams: axis x from the left end, transverse deflection v positive up, sagging moment positive,
  loads positive up (a downward load is negative).
- Pure NumPy/SciPy — no FEA, beam or composite libraries.

## Testing rules
- Every test asserts against a closed-form result or a cited published value; name the source
  in a comment. Never loosen a tolerance to make a failing solver pass without explaining why.
- Reference data in `validation/reference_data/` and `src/structures/materials/` must carry its
  citation (document + table number). Do not invent data: leave a value out if it cannot be
  checked against its source.

## Workflow
- Commit small logical changes directly to `main`, push after `ruff` + `pytest` pass.
