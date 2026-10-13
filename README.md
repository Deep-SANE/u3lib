# u3lib

**u3lib** is an **AI-assisted computational framework for U(3)/SU(3) algebra**:
irreps and their invariants, canonical and physical bases, generators and their
matrix elements, basis transformations, tensor-product decompositions, and the
**Wigner coupling** and **Racah recoupling** coefficients — exposed through a
C++ API, Python bindings, and an agent-driven skill layer over all of it.

- **Paper:** [Eur. Phys. J. Plus **139**, 704 (2024)](https://doi.org/10.1140/epjp/s13360-024-05501-8)
- **Manual:** [Manual (PDF)](doc/Manual.pdf); LaTeX source in `doc/Manual.tex`.

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
make          # build C++ driver programs in src/ (-> bin/)
make test     # recompile and run regression test suite (test/, -> bin/test/)
make python   # build the pybind11 extension -> python/u3lib.so (PYTHON=/path/to/python)
make clean
```

`make config` detects each dependency and lets you accept the detected default, point to your own install, or — for the two numerical libraries — download the pinned upstream source and build it locally (into `external-libs/`, which is not part of the repository). It is interactive by default; for scripted use pass flags through, e.g. `make config CONFIG_ARGS="--yes"`.

**Installing with an AI agent.** Point a code-executing assistant at this repository and prompt, for example:

> Install u3lib in this repository: follow `agents/install/SKILL.md` — run `make config` (download and build any missing dependencies), then `make` and `make test` — and report the test result.

The agent resolves the steps through the repository-local [install skill](agents/install/SKILL.md).

- **Regression Tests (`make test`):** Tests dimensions, basis interlacing, branching rules, coupling orthonormality, and documented analytic values. Checks remain enabled in release builds (`-UNDEBUG`).
- **Python Extension (`make python`):** Builds `python/u3lib.so`. See `python/tutorial.ipynb` for interpreter requirements and interactive usage.
- **Dependencies:** a C++20 compiler and OpenMP are required on the system; the library dependencies — Eigen3, Boost (headers), fmt, the **libwignernj** SU(2) Wigner backend, and the **sns** null-space solver (dense-QR fallback) — are each either detected, taken from a path you give, or **downloaded and built by `make config`** (into `external-libs/`, which is not committed). Fetching needs `git`, `cmake`, and `curl` or `wget`. Linux and macOS are both supported. Full configuration options are in [doc/Manual.pdf](doc/Manual.pdf).

---

## Repository Layout

- `include/` — header-only C++ library (`u3lib.h` umbrella; irreps, bases, generators, coupling, recoupling, containers).
- `src/` — C++ example driver programs; compiled into `bin/`.
- `python/` — pybind11 bindings (`pyu3lib.cpp`), build script, tutorial notebook (`tutorial.ipynb`), and notes.
- `test/` — regression suite (`test_*.cpp`, `*_cases.h`, and reference tables under `reference/`).
- `doc/` — full reference manual (`Manual.pdf` and LaTeX source).
- `agents/` — AI-agent skill catalog (`AGENTS.md`), runtime instructions, and task skills.
- `bin/` — compiled binary output (not source; git-ignored).

---

## Documentation & AI Agent Skills

- **Full Reference Manual:** For exhaustive mathematical derivations, phase conventions, getter signatures, and solver backends, consult [doc/Manual.pdf](doc/Manual.pdf).
- **AI-Assisted Calculations:** `u3lib` provides a repository-local skill package under [`agents/`](agents/AGENTS.md) that enables code-executing AI assistants to plan and drive calculations through natural-language prompts. A complete table of representative prompts with checked reference values is provided in **Table IV of the [Manual](doc/Manual.pdf)** as well as in [`agents/AGENTS.md`](agents/AGENTS.md).
