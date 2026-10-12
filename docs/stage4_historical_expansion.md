# Stage 4 - Historical expansion (2013-2022)

## Purpose

The 2023-2025 pilot cannot provide complete three-year lag features for the
temporal recurrence experiment: no target year in that window has three full
years of history behind it (2024 has one prior year, 2025 has two). Stage 4
extends the acquisition and case-resolution boundary backward so a target
year can eventually have `count_t1`, `count_t2`, and `count_t3` all computed
from real observed support rather than nulled out as unavailable.

## What Stage 4 actually required (vs. what needed to be built)

Every stage from MS4A through MS7B already accepts an arbitrary `--years`
list or range with no hardcoded window in code -- the vocabulary builder's
own `README.md` reproduction command already shows `--years 2013-2025`. The
**only** structural blocker was that `config/faers_sources.csv` had rows for
2023-2025 only. No pipeline source code required changes.

Two policy questions did need real (not assumed) answers, both confirmed
against actual downloaded 2016 data and the real MS4A resolver:

1. **Is 2013-2022 the same FAERS ASCII format as 2023-2025?** Yes for 2016
   (spot-checked; see below). FDA switched from "Legacy AERS" to "FAERS"
   format after 2012 Q3 (the last AERS file, partial to 27 Aug 2012; the first
   FAERS file is 2012 Q4), so 2013-2022 should be uniformly compatible, but only
   2016 had been verified at this point. The later full raw profile
   (`reports/data_quality/raw_profile_2013_2025/`) found one schema change
   inside the window, at 2014 Q3 (new fields, e.g. DRUG `prod_ai`).
2. **Does global CASEID resolution still behave correctly when the loaded
   year range widens?** Yes, with one non-obvious caveat documented below.

## Supported range after this change

`config/faers_sources.csv` now has verified rows for **all of 2013-2025**
(52 quarters). Every row was confirmed by an actual successful download and
ZIP-member validation, not inferred from a pattern.

## Manifest URL provenance

FDA's export page does not follow one fixed URL/casing pattern across years
(2023-2025 already mixes lowercase `q1/q2` with uppercase `Q3/Q4` within the
same year). A first WebFetch of FDA's export page produced a suspiciously
uniform 40-row table for 2013-2022 that was **not trusted and not used
directly** -- it reproduced the exact 2023-2025 casing quirk in a way
consistent with model pattern-completion rather than genuine extraction. It
was used only as a starting *candidate* per year, and every single year was
then verified by an actual download attempt (`download_faers.py` fails loud
on a wrong URL or a missing/mismatched archive member, so a bad guess is
never silently accepted). Confirmed casing per year:

| Years | Q1 | Q2 | Q3 | Q4 |
|---|---|---|---|---|
| 2013-2018 | q | q | q | q |
| 2019-2021 | Q | Q | Q | Q |
| 2022 | q | q | Q | Q |
| 2023-2024 | q | q | q | Q |
| 2025 | q | q | q | Q |

The WebFetch candidate table's own guesses turned out right for every year
except that it is not something to trust blind -- one bug in the process
that regenerated the manifest programmatically briefly mis-cased 2024 Q3 by
copying 2023's pattern instead of 2024's actual pattern; caught via `git
diff` against the original file before any download was attempted, not
after. That's the whole reason each year gets a live download, not just a
pattern match.

One genuine one-off historical naming exception, found only because the
download failed loudly: **2018 Q1**'s DEMO file is named
`DEMO18Q1_new.txt`, not `DEMO18Q1.txt` -- FDA republished a corrected file
for that quarter under a different name. Fixed in the manifest row's
`expected_inner_demo` field with a note explaining why.

## Schema compatibility (2016, verified)

`ascii/DEMO16Q1.txt` header:
`primaryid$caseid$caseversion$i_f_code$event_dt$mfr_dt$init_fda_dt$fda_dt$rept_cod$...`

`ascii/DRUG16Q1.txt` header:
`primaryid$caseid$drug_seq$role_cod$drugname$prod_ai$val_vbm$route$...`

`ascii/REAC16Q1.txt` header:
`primaryid$caseid$pt$drug_rec_act`

All columns required by `process_faers_case_level.py`'s loaders
(`caseid`, `primaryid`, `caseversion`, `fda_dt`, `init_fda_dt` for DEMO;
`primaryid`/`drugname`/`prod_ai`/`role_cod` for DRUG; `primaryid`/`pt` for
REAC) are present. The loaders already fail loud (`_require_columns`) on a
genuinely incompatible schema, so no new coercion logic was added or needed.

**2013 (oldest year) spot-check**, to bookend 2016: `DEMO13Q1.txt` has all
required columns. `DRUG13Q1.txt` is missing `prod_ai` entirely (that column
was added to the FAERS DRUG table at some point after 2013), and
`REAC13Q1.txt` is missing `drug_rec_act`. Neither breaks anything: `prod_ai`
and `drug_rec_act` are both in `process_faers_case_level.py`'s
*warn-and-default* list, not its *required* list -- `build_case_drugs`
already fills a missing `prod_ai` column with `""` before use
(`process_faers_case_level.py` around line 896-898). This is exactly the
graceful-degradation behavior the loaders were already designed for; no
change was needed to handle it.

## Global CASEID resolution behavior after widening the window

Two real MS4A runs (both isolated to
`data/interim/case_resolution_stage4_smoke/`, never touching the production
2023-2025-only index) provide before/after evidence:

**Run 1** -- `--years 2016 2023 2024 2025` (partial, 4 years):

| Year | Included cases | vs. documented 2023-2025-only production |
|---:|---:|---:|
| 2016 | 1,164,345 | n/a |
| 2023 | 1,361,395 | identical |
| 2024 | 1,308,772 | -1 |
| 2025 | 1,295,857 | identical |

**Run 2** -- `--years 2013 2014 ... 2025` (full, all 13 years):

| Year | Included cases | vs. documented 2023-2025-only production | vs. Run 1 |
|---:|---:|---:|---:|
| 2013 | 738,292 | n/a | n/a |
| 2014 | 832,546 | n/a | n/a |
| 2015 | 1,162,457 | n/a | n/a |
| 2016 | 1,164,345 | n/a | identical to Run 1 |
| 2017 | 1,211,308 | n/a | n/a |
| 2018 | 1,405,225 | n/a | n/a |
| 2019 | 1,413,737 | n/a | n/a |
| 2020 | 1,443,823 | n/a | n/a |
| 2021 | 1,535,120 | n/a | n/a |
| 2022 | 1,511,560 | n/a | n/a |
| 2023 | 1,361,395 | identical | identical |
| 2024 | 1,308,772 | -1 | identical |
| 2025 | 1,295,856 | -1 | **-1 vs. Run 1's 1,295,857** |

All 12 published invariants passed in both runs. 2023 stays byte-identical
to production in both. 2024 stays off by exactly -1 in both (a stable,
reproducible effect, not noise). **2025 was unaffected in Run 1 but shifted
by -1 once the full 2013-2025 range was loaded in Run 2** -- direct evidence
that even a year far from the newly added history can shift slightly once
the complete range is present, not just years adjacent to the new data. This
is consistent with the documented resolution policy, not a bug: a case's
`assigned_report_year` is the earliest valid `INIT_FDA_DT` among *all loaded*
observations, so adding earlier years can occasionally supply a valid
`INIT_FDA_DT` for a case that previously had to fall back to a less precise
assignment method, in any year, not only the newly loaded ones.

Full-range run cost: 19,478,153 DEMO rows across 52 quarters, 16,779,363
unique CASEIDs, 16,384,436 included, in 5,583 seconds (~93 minutes) total
runtime -- useful as a sizing reference for the eventual production
promotion.

`excluded_no_eligible_version` grew substantially for the newly-added years
once the full range replaced the partial one (e.g. 2016: 133 in Run 1 vs.
21,719 in Run 2). This is the same INIT_FDA_DT-visibility mechanism working
in the other direction: a case that Run 1 classified as
`assigned_outside_study_window` (because its true origin year wasn't in the
requested set) can become `excluded_no_eligible_version` in Run 2 once its
origin year *is* requested but no eligible row happens to be loaded for it.
`assigned_outside_study_window` fell from 90,172 (Run 1) to 113,905 across a
much larger population (Run 2) -- proportionally much smaller, as expected
once most of the historical window is actually loaded.

**Non-obvious policy consequence confirmed by this run and now pinned down
by `test_stage4_multi_year_gap_caseversion_resolves_to_single_case`** in
`tests/test_global_case_index.py`: eligibility requires
`source_year <= assigned_year`. Because `INIT_FDA_DT` is invariant across a
real case's versions, a version reported in a *later* year than the case's
true original year is never eligible to win, no matter how high its
`CASEVERSION` is -- it is recorded as an audited followup only. Resolution
freezes to the best eligible version reported at or before the case's true
origin year. This was already the documented policy
(`docs/cross_year_case_resolution.md`); Stage 4 doesn't change it, but it
means a case first reported in, say, 2016 and only amended in 2019+ will
keep resolving to its 2016 version even after 2019 is loaded.

2016 alone (without 2013-2015 loaded) loses 87,770 of its 1,252,115 raw
CASEIDs to `assigned_outside_study_window`: their `INIT_FDA_DT` predates
2016, which the resolver can see even without the earlier archive loaded
(the field is carried on the row itself). This is direct evidence that
**Stage 4 should be done as one contiguous 2013-2022 expansion**, not
year-by-year, or a meaningful slice of genuinely-2016-and-later cases will
be artificially excluded until 2013-2015 are also loaded.

## Reproduction

```powershell
# Verify one year's manifest rows before trusting them (never skip this)
.venv\Scripts\python.exe src\download_faers.py --years 2013 --dry-run
.venv\Scripts\python.exe src\download_faers.py --years 2013

# Full-range cross-year smoke test -- ALWAYS use isolated --out-dir/--report-dir,
# never the defaults, or you will overwrite the production 2023-2025 index
.venv\Scripts\python.exe src\build_global_case_index.py `
  --years 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 `
  --out-dir data/interim/case_resolution_stage4_smoke `
  --report-dir reports/data_quality/stage4_smoke

# MS4B case-resolution boundary for the newly-added years, same isolation rule
.venv\Scripts\python.exe src\process_faers_case_level_annual.py `
  --years 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 --case-resolution-only `
  --global-case-assignments data/interim/case_resolution_stage4_smoke/global_case_assignments.parquet `
  --report-dir reports/data_quality/stage4_smoke `
  --out-dir data/processed/faers_case_level_annual_stage4_smoke
```

Note: on this machine, plain `python` on PATH resolves to the system
interpreter, not the project's `.venv` (which has pandas/pyarrow/duckdb/
psycopg2 installed) -- use `.venv\Scripts\python.exe` explicitly, or the
scripts fail with `ModuleNotFoundError` before doing anything.

## Not yet done

- All 2013-2025 archives are downloaded and verified, but MS4A/MS4B have
  only been run against the full 2013-2025 window in an **isolated sandbox
  directory** (`data/interim/case_resolution_stage4_smoke/`), never against
  the production output paths. Promoting this to the real production
  `data/interim/case_resolution/` index is a separate, explicit step -- see
  `reports/data_quality/stage4_smoke/` for the sandbox run's results once it
  finishes.
- MS4B has not yet been run for the 10 newly-added years (2013-2015,
  2017-2022, plus 2016) to produce `resolved_drug_rows.parquet` /
  `resolved_reac_rows.parquet` -- this is required before MS5 can be
  rebuilt, since MS5 reads MS4B's per-year output.
- MS5 (vocabulary + ChEMBL mapping) has not been rebuilt for any historical
  year. A live ChEMBL PostgreSQL connection was confirmed reachable and
  working in this environment, but the rebuild itself was not run.
- MS6, MS7A, and MS7B have not been rerun for any historical year.
- The first target year with three complete lags, once 2013-2022 is fully
  loaded, is **2016** (`2016-3=2013 >= first_year 2013`) -- computed directly
  from `build_temporal_features.py`'s `lag()` logic, not assumed, and
  confirmed by the new `test_stage4_three_complete_lags` test.

## Promotion and downstream completion

All of the above is now done. The sandbox's validated full 2013-2025 MS4A
result was promoted to the production `data/interim/case_resolution/` index
(backup preserved at `data/interim/case_resolution_pre_stage4_backup/`), and
validated in place -- SHA-256 integrity, schema, CASEID uniqueness, and all
13 years' included-case counts matched the sandbox's Run 2 table above
exactly. No MS4A recomputation was needed.

MS4B, MS5, MS6, and MS7A were then built for 2013-2022, additively --
2023-2025's artifacts were left completely untouched (verified by checksum).
MS7B was rebuilt across the full 2013-2025 range, and MS7C's baselines and
logistic regression were rerun across every target year with three real lags
(2016-2025). Full results, methodology, and interpretation are in
[Historical Baseline v1](historical_baseline_v1.md).

Three genuine pre-existing pipeline bugs were found and fixed while building
this dataset, all newly exposed by the wider year range (each has a
regression test):

1. `integrate_global_case_assignments.py` raised a false-positive "later
   follow-up PRIMARYIDs reached downstream" assertion for ~340/16.4M included
   cases (2013-2021 only) where FAERS redundantly re-exported the identical
   `(CASEID, CASEVERSION)` row in a later year -- the case's own winning
   PRIMARYID was mistaken for a leaked followup of itself.
2. `faers_utils.py`'s FAERS ASCII table loader assumed UTF-8; 2019 Q3's DRUG
   file contains a Windows-1252 byte (`0xAD`, soft hyphen) that isn't valid
   UTF-8. Now falls back to `cp1252` only on a confirmed UTF-8 decode failure.
3. `build_multi_year_drug_vocabulary.py` hard-required a `prod_ai` column
   that doesn't exist at all in 2013's FAERS DRUG table (added in a later
   FAERS schema revision) -- fixed to tolerate its absence exactly like
   `process_faers_case_level.py`'s `_choose_raw_drug_name()` already does.

Two additional issues were fixed in the *shared, multi-year* output paths
(`rebuild_global_downstream.py`'s `annual_summary_all_years.json` and
`build_gold_annual_pair_reactions.py`'s whole-`output_root` staging) so that
building new years additively can never destroy already-published years --
the Gold-dataset one in particular would have deleted the entire validated
2023-2025 Gold dataset the first time an additive year range was published
with `--overwrite`. See [Historical Baseline v1](historical_baseline_v1.md)
for full detail.
