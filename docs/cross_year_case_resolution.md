# Cross-year FAERS case resolution policy

Independent yearly resolution can count a CASEID again when a later archive
contains only a follow-up. Milestone 4A builds a DEMO-only global audit index;
Milestone 4B integrates its selections into annual DRUG/REAC processing.

The project freezes assignment to the initial FDA receipt year instead of
assigning every observed year, the latest-version year, or shifting a case to
its first loaded partition. `assigned_report_year` uses the earliest valid
`INIT_FDA_DT`, then earliest valid `FDA_DT`, then first source partition.
Assignments outside the requested window remain auditable but excluded.

A row is eligible when `source_year <= assigned_report_year`; FDA_DT is not an
eligibility condition. Winners rank valid/highest CASEVERSION, latest source
year and quarter, latest valid FDA_DT, numeric then lexical PRIMARYID, and
stable provenance. Missing versions are not zero. A case with no eligible row
is excluded rather than shifted.

Cross-CASEID PRIMARYID reuse is critical and quarantined. Other date, version,
chronology, and tie conditions are reported. A false
Follow-up metrics have distinct scopes:

- `has_additional_version_in_loaded_data` means the CASEID has more than one
  deduplicated DEMO observation, including same-year versions.
- `has_later_year_followup_in_loaded_data` means it appears in multiple loaded
  source years and at least one row has `source_year > assigned_report_year`.
- `has_later_followup_in_loaded_data` is the legacy all-status audit predicate
  `source_year > assigned_report_year`. Outside-window and no-eligible cases
  may satisfy it even when they appear in only one loaded source year.
- `cases_removed_as_cross_year_followups` counts annual appearances, not unique
  CASEIDs globally, so one CASEID can contribute in multiple later years.

## Milestone 4B integration boundary

For analytical year Y, the annual pipeline uses only included MS4A assignments
whose `assigned_report_year` is Y. It filters raw DRUG and REAC rows by the
assignment's `selected_primaryid` and attaches the assignment CASEID as the
authoritative identifier. Selected source year and quarter route reads to the
physical archive only; they do not determine analytical inclusion. Later
follow-up PRIMARYIDs are excluded.

Included cases remain in `retained_cases.parquet` even when DRUG, REAC, or
both tables have no matching row; those conditions are counted in the MS4B
data-quality report. Follow-up audit fields are never inclusion criteria.
Reaction normalization, drug normalization and ChEMBL mapping, case-level
deduplication, pair construction, and support aggregation remain unchanged.

## Missing selected REAC audit

The MS4B production integration found three selected 2024 cases without a
resolved REAC row. A read-only, chunked audit scanned all 12 loaded 2023–2025
REAC members using the production delimiter, string typing, column
normalization, and identifier normalization. All three selected PRIMARYIDs
were absent from both their expected 2024 Q3 member and every other loaded
REAC member; this is not a routing, parsing, or identifier-normalization error.

Two logical cases have REAC data under later, nonselected 2025 PRIMARYIDs.
Those rows remain excluded because source-partition eligibility selected the
2024 version and MS4B never substitutes or borrows another version's data.
The third logical case has no REAC rows under any DEMO version loaded for the
study. The cases remain part of the globally retained MS4A population but do
not contribute to case_reactions, drug-pair–reaction records, or reaction
support counts.

All assignment statuses are included in these audit counts unless a report
field explicitly says otherwise. A false loaded-data flag proves nothing about
unloaded years.
