# Investigation Report: preprocess_genome_tags Parameter

## Summary
After a comprehensive investigation of the `long-read-qc` repository, **no evidence of a parameter named `preprocess_genome_tags` was found** in any part of the repository's history or current codebase.

## Investigation Methodology

### 1. Current Codebase Search
- **Action**: Searched all files for `preprocess_genome_tags`
- **Command**: `grep -r "preprocess_genome_tags" .`
- **Result**: No matches found

### 2. Git History Analysis
- **Repository Status**: Successfully unshallowed repository to access complete history
- **Total Commits Analyzed**: 235 commits
- **Searches Performed**:
  - `git log --all --full-history -S "preprocess_genome_tags"` - No results
  - `git log --all -p | grep "preprocess_genome_tags"` - No results
  - Searched commit messages for related terms like "preprocess", "genome", and "tag"

### 3. Related Parameters Found
While `preprocess_genome_tags` does not exist, the following related parameters were found in the codebase:

#### Genome-related Parameters (from `nextflow.config` and `nextflow_schema.json`):
- **`map_genome`** (boolean, default: false) - Map all reads to genome
- **`genomes`** (string, default: "s3://pioneer-data/genomes/") - Genomes directory path

#### Available Parameters in the Pipeline:
```
--samplesheet <file>      Path to the sample sheet
--outdir <dir>            Output directory
--tech <str>              Sequencing technology (map-ont/map-pb/map-hifi)
--map_genome <bool>       Map all reads to genome
--error_rate <float>      Error rate for barcode searching
--min_overlap <int>       Minimum overlap for barcode searching
--min_bc_len <int>        Minimum barcode length
--max_bc_len <int>        Maximum barcode length
--meta_ovlp <int>         Overlap bp for sourmash
--sourmash_db <file>      Path to the sourmash database
--taxonomy <file>         Path to the taxonomy database
--cores <int>             Number of cores to use
```

### 4. GitHub Issues and Pull Requests
- **Action**: Searched all issues and PRs for `preprocess_genome_tags`
- **Result**: No historical issues or PRs mentioning this parameter

### 5. Relevant Commits Found
While searching for related terms, the following potentially relevant commits were identified:
- **c91d0ab**: "Merge pull request #28 from Pioneer-Research-Labs/2025_05_07_tagmentation_issues"
- **3d45386**: "Adding flag for mapping all reads to genome and corresponding summary tables"
- **525dcf0**: "Added tests for mapping all reads to genome"

None of these commits reference `preprocess_genome_tags`.

## Conclusion

The parameter `preprocess_genome_tags` **does not exist** and **has never existed** in the `long-read-qc` repository based on:
1. Complete absence from current codebase
2. No appearances in git history spanning 235 commits
3. No mentions in GitHub issues or pull requests
4. No references in configuration files or documentation

## Possible Scenarios

1. **Parameter may have been confused with another parameter**: Perhaps `map_genome` or another genome-related parameter was intended?
2. **Parameter from a different project**: The parameter might exist in a different repository or pipeline
3. **Planned but never implemented**: The parameter may have been discussed but never actually added to the codebase
4. **Typo in the original query**: The parameter name might have been misspelled

## Recommendations

If you're looking for functionality to preprocess or tag genomes, consider:
- Using the existing `map_genome` parameter for genome mapping functionality
- Checking if this parameter exists in a related pipeline or different repository
- Clarifying the intended functionality to suggest alternative approaches
