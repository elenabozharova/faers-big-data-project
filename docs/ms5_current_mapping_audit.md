# MS5 current drug-vocabulary and ChEMBL mapping audit

## Scope and conclusion

This is an inspection of the repository at the MS4B-to-MS5 boundary. No
mapping was generated or changed. MS5 must consume only globally resolved
DRUG rows. The current generator is not yet a multi-year MS5 generator: its
default input is one raw, unfiltered 2025 Q4 archive, and no production mapping
artifact is currently present under `data/processed/faers/`. Consequently,
coverage of the corrected 2023–2025 vocabulary cannot currently be measured
or claimed.

## Current selected-name logic

The active downstream selector is `_choose_raw_drug_name()` in
`src/process_faers_case_level.py` (around line 831). Equivalent logic exists in
`map_faers_to_chembl.build_unique_name_table()` and
`drug_normalization.run.select_drug_field()`.

1. Missing `prod_ai` or `drugname` columns become an all-empty string Series.
2. Values are null-filled with `""`, cast to string, and trimmed.
3. `prod_ai` is meaningful when its length is greater than two and the value,
   after removing literal periods, is not numeric.
4. Meaningful `prod_ai` wins; otherwise `drugname` is selected.
5. When both are populated, a meaningful `prod_ai` wins without comparing the
   two values. When both are blank, the chosen name is blank.

`build_case_drugs()` preserves `drugname`, `prod_ai`, and `chosen_raw_name` in
`case_drugs_long.parquet`. The resolved MS4B rows retain the two original raw
fields but do not yet contain `chosen_raw_name`. The rule is independent of
year and resolution mode.

## Current normalization contract

The ChEMBL path calls `normalize_drug_names.normalize_drug_string()`
(`src/normalize_drug_names.py`, around line 120). This is distinct from the
RxNorm-oriented code under `src/drug_normalization/`; that alternative removes
salt variants and performs NFKD processing, but is not called by the active
ChEMBL pipeline.

The active function executes these transformations in order:

1. Falsey and non-string inputs return `""`.
2. The string is lowercased and trimmed. There is no explicit Unicode
   normalization or diacritic removal.
3. Numeric dosage patterns are replaced with spaces before combination
   splitting. The pattern accepts bare integers and decimal/range/ratio forms,
   with optional `mg`, `mcg`, `ug`, `ml`, `g`, `l`, `iu`, `meq`, or `%`.
4. `/`, `+`, and the whole words `and`, `with`, `plus`, and `without` become
   the combination separator `|`.
5. All punctuation other than `|` is replaced with spaces. Despite a source
   comment saying word-internal hyphens are preserved, the implemented regex
   removes hyphens too.
6. Each combination segment is processed separately. Standalone hyphens are
   removed, whitespace is split, bare digit tokens are removed, and configured
   noise tokens are removed. Noise includes measurement units, dosage forms,
   routes, release modifiers, pack/presentation descriptors, and
   `nos`/`nec`/`hcpcs`.
7. Empty segments are dropped and remaining segments are joined with `|`.

Salt names such as hydrochloride, succinate, and sodium are intentionally
retained. Parentheses have no special rule: punctuation is removed but their
contents remain. Combination products remain component strings rather than a
single product identifier. There is no explicit brand-name transformation or
special-token rejection beyond the noise-token list.

The implementation is deterministic, row-order/year independent, blank-safe,
and appears idempotent for outputs produced by its own transformations. The
repository has examples and downstream regression tests, but no focused unit
suite exhaustively proving Unicode behavior, idempotence, punctuation edge
cases, collision safety, or every noise token. The separate RxNorm smoke test
does not establish the active ChEMBL normalizer's full contract.

## Current ChEMBL mapping contract

### Paths and construction

The configured output directory is `data/processed/faers/`:

- `chembl_name_index.parquet`: cached ChEMBL preferred names and synonyms;
- `faers_drug_to_chembl.parquet`: one row per normalized FAERS name;
- `unmatched_drugs.csv`, `ambiguous_drugs.csv`, and `mapping_stats.json`.

None of these artifacts is present in the inspected worktree. The runtime
loader is `load_drug_mapping_file()` in `src/process_faers_case_level.py`. It
joins on mapping column `normalized` versus row column
`normalized_drug_name`; duplicate normalized keys are warned about and the
first is retained.

The generated mapping schema consists of `raw_name`, `raw_ai`, `field_used`,
`normalized`, `faers_count`, `chembl_id`, `match_method`, `match_score`,
`candidate_count`, `is_ambiguous`, `candidates`, and `match_detail`.
Combination rows may also populate `component_methods`.

### Matching order and values

The methods, in execution order, are:

- `exact_pref`: exact normalized preferred-name match, score 1.0;
- `exact_synonym`: exact normalized synonym match, score 1.0;
- `combo_all`: every `|` component matches exactly and unambiguously, score
  1.0, component IDs stored in `candidates`, but no single `chembl_id`;
- `fuzzy`: unmatched names only, at least one shared token followed by
  `difflib.SequenceMatcher.ratio()`, default threshold 0.85;
- `unmatched`;
- `empty_normalized`.

Exact and fuzzy candidate IDs and scores are preserved. Multiple unique
candidates set `is_ambiguous=True` and `chembl_id=None`. There is no manual
mapping tier. One normalized key is intended to occupy one artifact row; it
may retain multiple candidate IDs but cannot yield multiple accepted IDs. One
ChEMBL ID can be associated with many normalized FAERS names.

`build_case_drugs()` emits these row statuses:

- `mapped`;
- `unmapped_blank_name`;
- `unmapped_no_match`;
- `unmapped_ambiguous`;
- `excluded_fuzzy` (the default unless `include_fuzzy=True`).

There is an implementation caveat: `is_ambiguous` is evaluated only for rows
that already have a non-null `chembl_id`, while ambiguous mapping rows set
`chembl_id=None`. Thus ambiguous rows currently fall into
`unmapped_no_match`, making `unmapped_ambiguous` effectively unreachable for
artifacts produced by the current generator. Likewise, `combo_all` has no
single `chembl_id` and is treated as `unmapped_no_match` downstream. This
should be resolved deliberately during MS5 design rather than silently
altered.

Configuration defaults to database `chembl_36`, so ChEMBL 36 is the intended
source. Neither the cache nor mapping output embeds or proves a ChEMBL
version. No mapping-artifact or normalization version is recorded.

## Mapping-population provenance

`map_faers_to_chembl.py` reads exactly one ZIP/member selected by
`FAERS_ZIP_PATH` and `FAERS_DRUG_FILE`. Repository defaults are
`data/raw/faers/faers_ascii_2025Q4.zip` and `ASCII/DRUG25Q4.txt`. It reads raw
DRUG directly, without MS4A/MS4B filtering. Therefore a mapping generated
with repository defaults would represent one raw quarter, including
nonselected versions—not the globally corrected multi-year population.

The build procedure is reproducible in code if the configured ChEMBL database
or a valid name-index cache is available. However, no artifact, build
manifest, input hash, timestamped metadata, mapping version, or ChEMBL version
is present, so the provenance of any mapping used outside this worktree cannot
be proven. It is also impossible to determine whether such an external
artifact contains names first observed outside its original source population.

## Corrected DRUG inputs available

All three files were written by the MS4B case-resolution-only production run;
only retained/resolved boundary artifacts are present in each annual output
directory.

| Analytical year | Resolved DRUG path | Rows | Distinct CASEIDs | Assigned year column | Raw name fields | Mapping fields |
|---:|---|---:|---:|---|---|---|
| 2023 | `data/processed/faers_case_level_annual/2023/resolved_drug_rows.parquet` | 4,763,980 | 1,361,395 | No | `prod_ai`, `drugname` | None |
| 2024 | `data/processed/faers_case_level_annual/2024/resolved_drug_rows.parquet` | 5,008,214 | 1,308,773 | No | `prod_ai`, `drugname` | None |
| 2025 | `data/processed/faers_case_level_annual/2025/resolved_drug_rows.parquet` | 4,942,013 | 1,295,857 | No | `prod_ai`, `drugname` | None |

Every file also contains `primaryid`, `caseid`, `role_cod`, `drug_seq`,
physical `source_quarter`, and the remaining raw DRUG attributes. The full
common schema is `primaryid`, `drug_seq`, `role_cod`, `drugname`, `prod_ai`,
`val_vbm`, `route`, `dose_vbm`, `cum_dose_chr`, `cum_dose_unit`, `dechal`,
`rechal`, `lot_num`, `exp_dt`, `nda_num`, `dose_amt`, `dose_unit`,
`dose_form`, `dose_freq`, `source_quarter`, and authoritative `caseid`. No
chosen, normalized, status, method, or ChEMBL columns exist.

## Proposed MS5 input contract

Minimum columns should be `caseid`, `primaryid`, `prod_ai`, and `drugname`,
with `role_cod` retained when available. Add `assigned_report_year` from the
validated annual directory/year being read, then validate every CASEID against
the included MS4A assignment for that year. Do not infer analytical year from
physical `source_quarter`, and do not read raw or legacy annual DRUG rows.

## Proposed multi-year vocabulary schema

Use one row per nonblank `normalized_drug_name`:

- `normalized_drug_name`: active normalizer output and unique key;
- `first_seen_year` / `last_seen_year`: minimum/maximum analytical year;
- `years_seen`: sorted analytical-year list;
- `raw_row_count`: corrected resolved DRUG rows producing the key;
- `distinct_case_count`: distinct authoritative CASEIDs producing it;
- `distinct_raw_name_count`: distinct selected raw names producing it;
- `sample_raw_names`: deterministic, bounded, lexically sorted examples.

Blank selected or normalized names should be retained in quality counts but
not represented as a normal vocabulary key unless a documented sentinel is
adopted.

## Proposed versioned mapping schema

Retain supported current fields and add explicit provenance:

- `normalized_drug_name`, `chembl_id`, `mapping_status`, `mapping_method`;
- `mapping_score` (current `match_score`);
- `mapping_candidates` (current `candidates` plus candidate count/detail);
- `first_seen_year`, `last_seen_year`, `years_seen`, `raw_row_count`,
  `distinct_case_count` from the vocabulary;
- proposed additions: `mapping_artifact_version`, `chembl_version`, and
  `normalization_version`.

The three version fields require new explicit versioning; the current artifact
cannot supply them. Combination mappings need a deliberate representation
because a component-ID list is not a single molecular `chembl_id`.

## Proposed coverage and quality reports

1. Name coverage: mapped normalized keys divided by eligible normalized keys,
   split by status/method and year presence.
2. Row coverage: corrected resolved rows with an accepted mapping divided by
   eligible corrected resolved rows.
3. Exact distribution of mapping statuses and methods, including ambiguous,
   fuzzy-excluded, combinations, and blanks.
4. Top unmapped names by row count and distinct CASEID count.
5. Normalization collisions: normalized keys with multiple materially
   different selected raw names, frequencies, and samples.
6. ChEMBL fragmentation: ChEMBL IDs reached through multiple normalized names,
   stratified by exact/synonym/fuzzy method.
7. Case usability per analytical year: cases with zero, exactly one, and at
   least two distinct accepted mapped drugs. This is an audit only; do not
   generate pairs in MS5.

## Confirmed risks

- The existing generator's default one-quarter raw population is too narrow
  and semantically wrong for MS5.
- No local artifact exists, so current full-population coverage is zero in
  operational terms and quantitatively unknown in mapping terms.
- Aggressive removal of bare numbers and punctuation can create collisions;
  `without` becoming a combination separator is particularly questionable.
- Salt forms remain fragmented unless ChEMBL synonyms reconcile them.
- Brand coverage depends entirely on ChEMBL synonym contents.
- Combination matches are recorded but unusable by current single-ID logic.
- Fuzzy matching selects the highest scoring normalized name without a
  documented margin over the runner-up; ties depend on candidate iteration
  order and need deterministic review.
- Ambiguous status handling is inconsistent between generator and consumer.
- No ChEMBL, mapping-artifact, or normalization version is persisted.
- First-occurrence provenance in `build_unique_name_table()` depends on source
  row order; future multi-year samples must be deterministic.

## Recommended MS5 implementation sequence

1. Freeze and test the selected-name and normalization contracts, including
   collision and idempotence audits, without changing behavior silently.
2. Read only the three globally resolved DRUG files; attach and validate
   analytical year against MS4A assignments.
3. Build a deterministic multi-year vocabulary and its collision report.
4. Version the normalization, ChEMBL snapshot, input hashes, and artifact.
5. Apply exact preferred/synonym mapping first and report coverage.
6. Decide explicit policies for ambiguous, combination, salt, and fuzzy cases;
   keep fuzzy results reviewable and excluded by default.
7. Produce name-, row-, status-, fragmentation-, and case-usability reports.
8. Validate reproducibility and only then wire the versioned artifact into
   downstream case-drug processing. Pair construction belongs to MS6.

## Unresolved questions

- Which normalization behavior should be frozen where comments and regex
  behavior disagree, especially hyphens and `without`?
- Should multi-component products remain component lists, map to a ChEMBL
  product entity, or be excluded from single-molecule analyses?
- What deterministic fuzzy tie/margin and manual-review policy is acceptable?
- Is ChEMBL 36 the required frozen release, and how will its database/cache
  provenance be hashed and recorded?
- Should salts remain distinct vocabulary keys while sharing an ingredient
  mapping, or should MS5 add an explicitly versioned ingredient-normalization
  layer?
- What artifact-version convention and bounded sample size should be used?
