# u3lib

**u3lib** is an **AI-ready computational framework for U(3)/SU(3) algebra**:
irreps and their invariants, canonical and physical bases, generators and their
matrix elements, basis transformations, tensor-product decompositions, and the
**Wigner coupling** and **Racah recoupling** coefficients — exposed through a
C++ API, Python bindings, and an agent-driven skill layer over all of it.

- **Preprint:** [arXiv:2611.26558](https://arxiv.org/abs/2611.26558) — describes u3lib; see [Citing](#citing).
- **Method paper:** [Eur. Phys. J. Plus **139**, 704 (2024)](https://doi.org/10.1140/epjp/s13360-024-05501-8) — the procedure u3lib uses for the coupling and recoupling coefficients.
- **Manual:** [Manual (PDF)](https://deep-sane.github.io/u3lib/Manual.pdf), built automatically from the LaTeX source in `doc/Manual.tex`.

> **Note for AI Agents:** Do not parse this human documentation. Refer directly to [agents/AGENTS.md](agents/AGENTS.md) for repository guidance, first-time build/installation, runtime lifecycle, and task-specific skills.

---

## Scope

The framework covers the complete U(3) representation theory needed for symmetry-adapted structure calculations:

- **Irreps & invariants** — SU(3) $(\lambda,\mu)$ and U(3) $[n_{13},n_{23},n_{33}]$ labels, dimension, Casimir invariants $C_1, C_2, C_3$, and conjugation.
- **Bases** — canonical Gelfand–Tsetlin basis and physical SO(3)-reduced basis $(\kappa,L)$, with inner and outer multiplicities.
- **Generators** — canonical $E_{ij}$ matrix elements on Gelfand states, and SO(3)-reduced matrix elements of physical $\mathbf{L}$ and $\mathbf{Q}$.
- **Transformations** — physical $|(\lambda\mu)\kappa LM\rangle \leftrightarrow$ canonical Gelfand expansion.
- **Decomposition** — full Kronecker product series $A \otimes B$ with outer multiplicities.
- **Coupling & recoupling** — $\mathrm{U(3)}\supset\mathrm{U(2)}$ canonical and SO(3)-reduced Wigner coefficients, and $U$-, $Z$-, and $9$-$(\lambda\mu)$ Racah recoupling.

---

## Build & Test

u3lib is self-contained and needs no environment variables. Configure once, then build:

```bash
make config   # resolve all dependencies, writing config/config.generated.mk
make -j       # build C++ driver programs in parallel (src/ -> bin/)
make -j test  # recompile and run regression test suite in parallel (test/, -> bin/test/)
make python   # build pybind11 extension -> python/u3lib.so (auto-detects a pybind11+Jupyter python; override PYTHON=...)
make clean
```

`make config` detects each dependency and lets you accept the detected default, point to your own install, or fetch and build it. The **KraNS** null-space solver is used as a sibling repository next to u3lib (`../krans`), built in place if not already installed; `make config` clones it from Deep-SANE automatically if it is missing. Configuration is interactive by default; for scripted use pass flags through, e.g. `make config CONFIG_ARGS="--yes"`.

**Installing with an AI agent.** Point a code-executing assistant at this repository and prompt, for example:

> Install u3lib in this repository: follow `agents/install/SKILL.md` — run `make config` (download and build any missing dependencies), then `make` and `make test` — and report the test result.

The agent resolves the steps through the repository-local [install skill](agents/install/SKILL.md).

- **Regression Tests (`make test`):** Tests dimensions, basis interlacing, branching rules, coupling orthonormality, and documented analytic values. Checks remain enabled in release builds (`-UNDEBUG`).
- **Python Extension (`make python`):** Builds `python/u3lib.so`. The interpreter is auto-detected — one that has `pybind11` and, for the tutorial, `ipykernel` — and can be overridden with `make python PYTHON=/path/to/python`. See `python/tutorial.ipynb` for interactive usage.
- **Dependencies:** a C++20 compiler and OpenMP are required on the system; the library dependencies — Eigen3, Boost, fmt, libwignernj, and KraNS — are each either **detected**, **taken from a path** you give, or **fetched and built by `make config`**. Fetching needs `git`, `cmake`, and `curl` or `wget`. Linux and macOS are both supported. Full configuration options are in the [Manual (PDF)](https://deep-sane.github.io/u3lib/Manual.pdf).

---

## Repository Layout

- `include/` — header-only C++ library (`u3lib.h` umbrella; irreps, bases, generators, coupling, recoupling, containers).
- `src/` — C++ example driver programs; compiled into `bin/`.
- `python/` — pybind11 bindings (`pyu3lib.cpp`), build script, tutorial notebook (`tutorial.ipynb`), and notes.
- `test/` — regression suite (`test_*.cpp`, `*_cases.h`, and reference tables under `reference/`).
- `doc/` — LaTeX source of the full reference manual (`Manual.tex`); the PDF is built on every change and published on [GitHub Pages](https://deep-sane.github.io/u3lib/Manual.pdf).
- `agents/` — AI-agent skill catalog (`AGENTS.md`), runtime instructions, and task skills.
- `bin/` — compiled binary output (not source; git-ignored).

---

## Documentation & AI Agent Skills

- **Full Reference Manual:** For exhaustive mathematical derivations, phase conventions, getter signatures, and solver backends, consult the [Manual (PDF)](https://deep-sane.github.io/u3lib/Manual.pdf).
- **AI-Ready Calculations:** `u3lib` provides a repository-local skill package under [`agents/`](agents/AGENTS.md) that enables code-executing AI assistants to plan and drive calculations through natural-language prompts. A complete table of representative prompts with checked reference values is provided in **Table IV of the [Manual](https://deep-sane.github.io/u3lib/Manual.pdf)** as well as in [`agents/AGENTS.md`](agents/AGENTS.md).

---

## Citing

If you use u3lib in your work, please cite the u3lib preprint:

> D. Langr, P. Dang, T. Dytrych, M. Negro: **u3lib: An AI-Ready Computational Framework for
> U(3) Symmetry Algebra, Wigner Coupling, and Racah Recoupling Coefficients.** arXiv:2611.26558
> [nucl-th], 2026.
> [doi:10.48550/arXiv.2611.26558](https://doi.org/10.48550/arXiv.2611.26558)

If your results use the Wigner coupling or Racah recoupling coefficients, please also cite the
paper that describes the procedure u3lib implements:

> P. Dang, J. P. Draayer, F. Pan, K. S. Becker: **New procedure for evaluation of U(3) coupling and
> recoupling coefficients.** Eur. Phys. J. Plus **139**, 704 (2024).
> [doi:10.1140/epjp/s13360-024-05501-8](https://doi.org/10.1140/epjp/s13360-024-05501-8)

```bibtex
@misc{Langr2026u3lib,
  author        = {Langr, Daniel and Dang, Phong and Dytrych, Tom{\'a}{\v{s}} and Negro, Michela},
  title         = {u3lib: An {AI}-Ready Computational Framework for {U(3)} Symmetry Algebra,
                   {Wigner} Coupling, and {Racah} Recoupling Coefficients},
  year          = {2026},
  eprint        = {2611.26558},
  archivePrefix = {arXiv},
  primaryClass  = {nucl-th},
  doi           = {10.48550/arXiv.2611.26558},
  url           = {https://arxiv.org/abs/2611.26558}
}

@article{Dang2024U3,
  author  = {Dang, Phong and Draayer, Jerry P. and Pan, Feng and Becker, Kevin S.},
  title   = {New procedure for evaluation of {U(3)} coupling and recoupling coefficients},
  journal = {Eur. Phys. J. Plus},
  volume  = {139},
  pages   = {704},
  year    = {2024},
  doi     = {10.1140/epjp/s13360-024-05501-8}
}
```
