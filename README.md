# genEXE

Miscellaneous standalone scripts for genomics bioinformatics tasks.

## Structure

| Directory | Contents |
|-----------|----------|
| `bash/`   | Shell scripts for data processing, pipeline automation, and system tasks |
| `R/`      | R scripts for statistical analysis, visualization, and genomics workflows |
| `python/` | Python scripts for data manipulation, parsing, and analysis |
| `data/`   | Small reference files, test datasets, and configuration files |
| `docs/`   | Additional documentation and notes |

## Scripts

<!-- Scripts are listed here as they are added -->

### Bash

| Script | Description |
|--------|-------------|
| [pvar_merge_check.sh](bash/pvar_merge_check.sh) | Compare two PLINK2 pvar files for merge compatibility. Outputs a TSV cataloguing all variants with match codes: `2` perfect, `1` semi-match (diff ALT), `-1` REF mismatch, `0` missing in one dataset. |

### R

*(none yet)*

### Python

*(none yet)*

## Usage

Each script is standalone. Refer to the header comments within each file for dependencies, usage, and examples.
