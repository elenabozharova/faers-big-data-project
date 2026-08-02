# Cross-year FAERS case resolution policy

Independent yearly resolution can count a CASEID again when a later archive
contains only a follow-up. Milestone 4A builds a DEMO-only global audit index;
Milestone 4B will integrate its selections into annual DRUG/REAC processing.

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

All assignment statuses are included in these audit counts unless a report
field explicitly says otherwise. A false loaded-data flag proves nothing about
unloaded years.
