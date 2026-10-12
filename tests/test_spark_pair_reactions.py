import sys
import unittest
from pathlib import Path

import pandas as pd
import pyarrow as pa

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import process_faers_case_level as pipeline  # noqa: E402
import spark_pair_reactions as spark_job  # noqa: E402
from normalize_drug_names import normalize_drug_string  # noqa: E402


class SparkBatchFunctionTests(unittest.TestCase):
    """The mapInArrow functions must agree with the local pipeline row for row (no Spark needed)."""

    def test_drug_name_batches_match_local_choice_and_normalization(self):
        prod_ai = ["ASPIRIN", "  ", "", "12.5", "AB", "½½½", "METFORMIN HCL", "x.y"]
        drugname = ["BAYER 81MG", "Tylenol 500 mg", "LIPITOR", "WARFARIN", "AMOXICILLIN/CLAV", "", "", "  "]
        frame = pd.DataFrame({"prod_ai": prod_ai, "drugname": drugname})
        chosen = pipeline._choose_raw_drug_name(frame)
        expected = [normalize_drug_string(name) if name.strip() else None for name in chosen]
        batch = pa.RecordBatch.from_pydict({"prod_ai": prod_ai, "drugname": drugname})
        (out,) = spark_job.drug_name_batches([batch])
        self.assertEqual(out.column(2).to_pylist(), expected)

    def test_reaction_pt_batches_match_build_case_reactions(self):
        pts = ["  Nausea ", "HEAD	ACHE", "", "drug  ineffective", "Rash Pruritic"]
        resolved = pd.DataFrame({"caseid": [str(i) for i in range(len(pts))], "primaryid": "1", "pt": pts})
        local, _ = pipeline.build_case_reactions(resolved)
        (out,) = spark_job.reaction_pt_batches([pa.RecordBatch.from_pydict({"pt": pts})])
        for caseid, value in zip(resolved["caseid"], out.column(1).to_pylist()):
            survives = value is not None and value not in spark_job.REACTION_PT_EXCLUDE
            self.assertEqual(local.set_index("caseid")["reaction_pt"].get(caseid), value if survives else None)


if __name__ == "__main__":
    unittest.main()
