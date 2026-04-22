# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**divref-wf** creates a DivRef-style resource bundle (FASTA sequences and DuckDB index) for common human variation, based on the [DivRef project](https://github.com/e9genomics/human-diversity-reference). It wraps original Python scripts with improved typing, parameterization, and unit tests, and adds Snakemake workflow orchestration.

## Environment & Commands

This project uses **two package managers**:

- **`uv`** — manages the `divref` Python package under `./divref/`
- **`pixi`** — manages the full workspace (Python toolkit + Snakemake + Hail)

### Python Toolkit (`./divref/`)

```bash
uv run --directory divref poe fix-and-check-all   # Fix then check everything (required before commit)
uv run --directory divref poe check-all            # Check format, lint, tests, and types
uv run --directory divref poe fix-all              # Auto-fix formatting and linting
uv run --directory divref pytest                   # Run tests
uv run --directory divref pytest tests/test_main.py::test_name  # Run a single test
uv run --directory divref mypy divref/             # Type-check only
```

### Workspace (Pixi)

```bash
pixi run fix-and-check-all   # Fix and check toolkit + Snakemake linting
pixi run lint --check        # Validate Snakemake files with snakefmt
pixi run download            # Run the Snakemake download workflow
pixi run setup-gcs           # Download GCS connector JAR (required once for Hail on GCS)
```

## Architecture

### Repository Layout

```
divref/                    # Python package (uv-managed)
  divref/
    main.py                # CLI entry point; registers tools in _tools list
    alias.py               # HailPath type alias (str; accepts local, gs://, hdfs://)
    defaults.py            # Package-wide constants: POPULATIONS, REFERENCE_GENOME, freq thresholds
    hail.py                # Hail initialization with GCS connector setup
    haplotype.py           # Shared Hail utilities for haplotype sequence/windowing
    tools/                 # One module per CLI subcommand
  tests/                   # pytest tests
  pyproject.toml           # Package deps, ruff/mypy/pytest config
workflows/                 # Snakemake workflows
  generate_divref.smk      # Main workflow (extract → haplotypes → reference download)
  create_test_data.smk     # Generates gnomAD subset for unit tests
  config/config.yml        # Workflow configuration (chromosomes, populations, paths)
pixi.toml                  # Workspace config (snakemake + hail environments)
```

### CLI Pattern

Tools are plain functions registered in `main.py`; **defopt** auto-generates the CLI from their docstrings. To add a new tool:

1. Create `divref/tools/<name>.py` with a keyword-only function and Google docstring
2. Import and add it to `_tools` in `main.py`

```bash
divref <tool-name> --arg value   # Invokes the registered tool
```

### Tool Pipeline (execution order)

The tools implement a data pipeline:
1. `extract_gnomad_afs` / `extract_gnomad_single_afs` → per-population allele frequency Hail table
2. `extract_sample_metadata` → simplified sample→population mapping table
3. `create_gnomad_sites_vcf` → VCF of variants above AF threshold (uses output of step 1)
4. `compute_haplotypes` → groups phased variants into haplotype windows using Hail
5. `compute_haplotype_statistics` → haplotype count distributions
6. `compute_variation_ratios` → per-sample variant counts at multiple freq thresholds
7. `create_fasta_and_index` → FASTA sequences + DuckDB index (final deliverable)
8. `remap_divref` → maps haplotype coordinates back to reference genome (post-CALITAS step)

`extract_gnomad_single_afs` is an alternative to `extract_gnomad_afs` supporting both gnomAD v4.1 (JOINT) and v3.1.2 (HGDP+1KG) table schemas.

### Key Shared Modules

**`haplotype.py`**
- `get_haplo_sequence(context_size, variants)` — builds haplotype sequence strings with flanking reference context; handles SNPs, insertions, deletions
- `split_haplotypes(ht, window_size)` — splits multi-variant haplotypes at gaps ≥ `window_size` bases; discards sub-haplotypes with <2 variants
- `variant_distance(v1, v2)` — reference bases between two variants (accounts for indel length)

**`compute_haplotypes.py` two-window strategy**: To avoid systematic edge artefacts, the tool runs two overlapping window passes (offset by `window_size / 2`) and unions the results. Intermediate `.1.ht` / `.2.ht` files are cleaned up after the merge.

**`hail.py`**: `hail_init(gcs_credentials_path)` — sets `GOOGLE_APPLICATION_CREDENTIALS`, verifies GCS connector JAR (installed via `pixi run setup-gcs`), then calls `hl.init()` with Spark GCS config.

**`defaults.py`**: `POPULATIONS`, `REFERENCE_GENOME`, `VARIATION_RATIO_FREQUENCY_THRESHOLDS` — defaults shared across tools.

### Data Models (`remap_divref.py`)

Pydantic `frozen=True` models: `Variant`, `ReferenceMapping`, `Haplotype` — used for type-safe coordinate remapping. `Haplotype` uses field aliases to match mixedCase column names in the DuckDB index created by `create_fasta_and_index`.

## Git Workflow

### Commit Granularity

Commit after completing one of:
- A single function/method implementation
- One refactoring step (rename, extract, move)
- A bug fix with its regression test
- A documentation update

**Size guidelines:**
- Per commit: 100–300 lines preferred, 400 max
- Per PR: No hard limit, but consider splitting if >800 lines or >5 unrelated files

**Good commit scope examples:**
- `Add FastaIndex.validate() method`
- `Rename species_map → species_to_ref_fasta_map`
- `Fix off-by-one in BED coordinate parsing`

### Commit Messages

Use [Conventional Commits](https://www.conventionalcommits.org/) for commit messages and PR
titles+bodies. Common types: `feat`, `fix`, `chore`, `docs`, `refactor`, `test`.

```
<type>: <imperative description> (<72 chars total)

Detailed body explaining:
- What changed
- Why (link issues with "Closes #123" or "Related to #456")
- Any non-obvious implementation choices
```

### Commit Rules
- Run `uv run poe fix-and-check-all` before each commit; all checks must pass
- No merge commits
- Do not rebase without explicit user approval
- **Never mix formatting and functional changes.** If unavoidable, isolate formatting into separate commits at start or end of branch.

### Branch Hygiene
- Use `.gitignore` liberally
- Never commit: IDE files, personal test files, local debug data, commented-out code

## Coding Conventions

### Organization
- Extract logic into small–medium functions with clear inputs/outputs
- Scope variables tightly; limit visibility to where needed
- Use block comments for visual separation when function extraction isn't practical

### Naming
- Meaningful names, even if long: `species_to_ref_fasta_map` not `species_map`
- Short names only for tight scope (loop indices, single-line lambdas)
- Signal behavior in function names: `to_y()`, `is_valid()` → returns value; `update_x()` → side effect

## Testing

### Principles
- Generate test data programmatically; avoid committing test data files
- Test behavior, not implementation—tests should survive refactoring
- Cover: expected behavior, error conditions, boundary cases
- Scale rigor to code longevity: thorough for shared code, lighter for one-off scripts

### Coverage Expectations
- New public functions: at least one happy-path test + one error case
- Bug fixes: add a regression test that would have caught the bug
- Performance-critical code: include benchmark or explain in PR why not needed

## Documentation Maintenance

When modifying code, update as needed:
- [ ] Docstrings (if signature or behavior changed)
- [ ] README.md (if usage patterns changed)
- [ ] Migration notes (if breaking change)

## Python-Specific

### Pragmatism
- Balance functional, OOP, and imperative—use what's clearest
- When in doubt, prefer pure functions and immutable data
- Know your utility libraries; contribute upstream rather than writing one-offs

### Style
- Heavier use of classes and type annotations than typical Python
- Prefer `@dataclass(frozen=True)` and Pydantic models with `frozen=True`

### Functions
- Functions should have **either** returns **or** side effects, not both
- Exceptions: logging, caching (where side effect is performance-only)

### Documentation
- Google-style docstrings with `Args:`, `Returns:`, `Yields:`, and `Raises:` blocks
- Docstrings are required on all public functions/classes
- Code comments should explain non-obvious choices and complex logic

### Typing
- **Required:** Type annotations on all function parameters and returns
- **Parameters:** Accept the most general type practical (e.g., `Iterable` over `List`)
- **Returns:** Return the most specific type without exposing implementation details
- Annotate locals when: they become return values, or called function lacks hints
- Use type aliases or `NewType` for complex structures
- Avoid `Any`—prefer type alias or `TypeVar`
- Avoid `cast()` and `type: ignore`—prefer alternatives, but when unavoidable (e.g., incorrect upstream stubs), document the reason inline.

