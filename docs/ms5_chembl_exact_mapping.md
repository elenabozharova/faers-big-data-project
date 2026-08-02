# MS5 versioned ChEMBL exact mapping

## Scope

MS5C–F build a verified ChEMBL 36 preferred/synonym name index, map the
globally corrected multi-year FAERS vocabulary with exact-only rules, report
coverage, and expose explicit statuses to the existing consumer. This stage
does not call RxNorm, use fuzzy/manual mappings, create ChEMBL combination
entities, split combinations into case drugs, generate pairs, or run the full
annual downstream pipeline.

## ChEMBL 36 index

`src/build_chembl_name_index.py` connects using the existing `config.py`
settings. Before labelling output as ChEMBL 36 it requires:

- the active database name to match `chembl_36`;
- `molecule_dictionary`, `molecule_synonyms`, and `version` in the active
  schema;
- the installed `version` table contents to identify release 36.

Preferred names come from `molecule_dictionary.chembl_id/pref_name`.
Synonyms come from `molecule_synonyms.synonyms`, joined to
`molecule_dictionary` by `molregno`. Both use the frozen
`normalize_drug_string()` function. Exact duplicate relationships are removed,
but every distinct normalized-name-to-ChEMBL-ID relationship survives so
ambiguity can be detected later.

The index schema is `normalized_chembl_name`, `chembl_id`, `source_name`,
`name_type` (`preferred` or `synonym`), `source_table`, `chembl_version`, and
`normalization_version`. It is deterministically ordered and versioned as:

```text
data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet
```

Rows are streamed from PostgreSQL, normalized in batches, deduplicated in
temporary SQLite, and streamed into a staged Parquet file. Passwords,
usernames, hosts, and connection strings are not written to metadata.

## Exact mapping policy

`src/build_faers_chembl_mapping.py` reads the existing 94,725-row vocabulary
and versioned index. For ordinary names, preferred candidates take precedence
over synonyms:

- one ID: `mapped` with `exact_pref` or `exact_synonym`;
- multiple IDs: `unmapped_ambiguous`, null `chembl_id`, and all sorted
  candidates;
- no candidates: `unmapped_no_match` / `unmatched`.

No fuzzy mapping is executed or accepted.

Combination names are split in their existing normalized component order.
Every component follows the same preferred-then-synonym exact policy:

- all uniquely mapped: `mapped_combination_components` / `combo_all_exact`;
- any ambiguous component: `unmapped_combination_ambiguous`;
- otherwise any unmatched component: `unmapped_combination_partial`.

Component IDs and canonical JSON details are retained, but `chembl_id` remains
null. Component matches are not exposed as independent case drugs.

The final mapping keeps exactly one row per vocabulary key and carries all
vocabulary frequency/year fields plus ChEMBL, normalization, vocabulary,
index, and mapping artifact versions.

## Consumer contract

`build_case_drugs()` now recognizes the versioned `normalized_drug_name`,
`mapping_status`, and `mapping_method` columns. Only `mapping_status=mapped`
with a nonblank `chembl_id` produces `mapped_drug_id`. Explicit ambiguous,
unmatched, and combination statuses remain visible in `case_drugs_long`; a
component-only combination never becomes an ordinary mapped case drug.

Legacy artifacts remain explicit: the old `normalized`/`match_method` schema
is accepted, and an old ambiguous row is now recognized from `is_ambiguous`
even though its `chembl_id` is correctly null. Old `combo_all` rows are not
silently reinterpreted because that schema lacks the explicit component
contract.

## Coverage reports

After a verified index exists, the mapper scans only projected annual columns
in PyArrow batches, reapplies the frozen selector/normalizer, and uses temporary
SQLite for exact distinct-case status counts. It writes:

- `ms5_chembl_mapping_coverage_by_year.csv`, including a true vocabulary-union
  `ALL` row;
- `ms5_chembl_mapping_status_by_year.csv`;
- `ms5_chembl_unmatched_names_top.csv`;
- `ms5_chembl_ambiguous_mappings.csv`;
- `ms5_chembl_combination_mapping.csv`;
- `ms5_chembl_mapping_summary.json`.

Single-ID name/row coverage excludes component-only combinations. Fully
recognized combinations and all unmapped categories are reported separately.
The optional case-usability report is deferred because producing it would add
another large exact case/drug aggregation; no drug pairs are needed or created.

## Commands

Set the local PostgreSQL password in the process environment, then build the
verified index:

```powershell
$env:PGPASSWORD = "<local ChEMBL database password>"
python src/build_chembl_name_index.py
```

Only after successful verification and index generation:

```powershell
python src/build_faers_chembl_mapping.py `
  --years 2023-2025 `
  --vocabulary data/interim/drug_mapping/multi_year_drug_vocabulary_2023_2025_norm_v1.parquet `
  --name-index data/interim/drug_mapping/chembl36_name_index_norm_v1_index_v1.parquet `
  --output data/processed/mappings/faers_to_chembl_2023_2025_chembl36_norm_v1_exact_v1.parquet
```

Existing outputs require `--overwrite`. All artifacts are staged before
publication and path metadata is repository-relative POSIX.

## Production status

The local `chembl_36` database was verified from both the active database name
and its `version` table, which identifies ChEMBL Release 36. The installed
schema contains the required `molecule_dictionary`, `molecule_synonyms`, and
`version` tables.

The production index contains 150,008 distinct name relationships, 102,527
normalized names, 50,614 ChEMBL IDs, and 3,160 normalized names associated
with more than one ID. The exact mapper emitted one row for each of the 94,725
eligible FAERS vocabulary names. It produced 8,624 single-ID mappings covering
13,138,103 of 14,714,102 eligible resolved DRUG rows (89.2892%). It also
reported 1,345 fully recognized combinations, 325 ordinary ambiguous names,
826 ambiguous combinations, 6,248 partial combinations, and 77,357 unmatched
names. Ordinary ambiguity and combination ambiguity remain separate in both
the coverage report and summary.
