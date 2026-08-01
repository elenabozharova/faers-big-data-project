import csv, hashlib, json, subprocess, sys, tempfile, unittest, zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
import profile_faers_raw as p

DEMO=" PRIMARYID $CASEID$CASEVERSION$FDA_DT$INIT_FDA_DT\n1$10$1$20240101$20231231\n1$10$1$bad$\n2$$x$$20240230\n"
DRUG="PRIMARYID$CASEID$DRUGNAME$PROD_AI$ROLE_COD$DRUG_SEQ\n1$10$A$$PS$1\n1$10$$$C$2\n1$10$B$AI$C$2\n"
REAC="PRIMARYID$CASEID$PT\n1$10$Headache\n1$10$ Headache \n1$10$\n"

class ProfileTests(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory(); self.root=Path(self.t.name); self.raw=self.root/"raw"; self.raw.mkdir(); self.out=self.root/"out"; self.man=self.root/"m.csv"
  self.write_zip(); self.write_manifest()
 def tearDown(self): self.t.cleanup()
 def write_zip(self, demo=DEMO, drug=DRUG, reac=REAC):
  with zipfile.ZipFile(self.raw/"q1.zip","w") as z:
   if demo is not None:z.writestr("ASCII/DEMO24Q1.txt",demo)
   if drug is not None:z.writestr("ASCII/DRUG24Q1.txt",drug)
   if reac is not None:z.writestr("ASCII/REAC24Q1.txt",reac)
 def write_manifest(self):
  fields=list(__import__('download_faers').MANIFEST_COLUMNS)
  row={"year":2024,"quarter":1,"source_system":"FAERS","format":"ASCII","url":"https://x/q1.zip","source_filename":"q1.zip","local_filename":"q1.zip","expected_inner_demo":"ASCII/DEMO24Q1.txt","expected_inner_drug":"ASCII/DRUG24Q1.txt","expected_inner_reac":"ASCII/REAC24Q1.txt","notes":""}
  with self.man.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerow(row)
 def args(self,*x): return ["--years","2024","--quarters","1","--manifest",str(self.man),"--raw-dir",str(self.raw),"--out-dir",str(self.out),*x]
 def rows(self):
  with (self.out/"raw_schema_profile.csv").open() as f:return list(csv.DictReader(f))
 def test_valid_archive_order_schema_and_quality(self):
  before=hashlib.sha256((self.raw/"q1.zip").read_bytes()).hexdigest(); self.assertEqual(p.main(self.args("--chunk-size","1")),0); rows=self.rows()
  self.assertEqual([r["table"] for r in rows],["DEMO","DRUG","REAC"]); self.assertEqual([int(r["row_count"]) for r in rows],[3,3,3])
  demo,drug,reac=rows; self.assertEqual(json.loads(demo["raw_column_names"])[0]," PRIMARYID "); self.assertEqual(json.loads(demo["normalized_column_names"])[0],"primaryid")
  self.assertEqual((demo["missing_caseversion"],demo["invalid_caseversion"]),("0","1")); self.assertEqual((demo["missing_fda_dt"],demo["invalid_fda_dt"]),("1","1")); self.assertEqual((demo["missing_init_fda_dt"],demo["invalid_init_fda_dt"]),("1","1"))
  self.assertEqual((demo["duplicate_primaryid_rows"],demo["duplicate_caseid_caseversion_rows"]),("1","1")); self.assertEqual(drug["duplicate_primaryid_drug_seq_rows"],"1"); self.assertEqual(drug["missing_both_drugname_and_prod_ai"],"1")
  self.assertEqual(reac["duplicate_primaryid_pt_rows"],"1"); self.assertEqual(reac["missing_pt"],"1"); self.assertEqual(before,hashlib.sha256((self.raw/"q1.zip").read_bytes()).hexdigest()); self.assertEqual(list(self.raw.iterdir()),[self.raw/"q1.zip"])
 def test_chunk_sizes_identical_and_signature_stable(self):
  self.assertEqual(p.main(self.args("--chunk-size","1")),0); a=self.rows(); self.assertEqual(p.main(self.args("--chunk-size","99")),0); b=self.rows()
  keys=["row_count","schema_signature",*p.METRICS]; self.assertEqual([[r[k] for k in keys] for r in a],[[r[k] for k in keys] for r in b]); self.assertNotEqual(a[0]["schema_signature"],a[1]["schema_signature"])
 def test_missing_column_null_metric_and_warning(self):
  self.write_zip(demo="PRIMARYID$CASEID$FDA_DT$INIT_FDA_DT\n1$1$20240101$20240101\n"); self.assertEqual(p.main(self.args()),0); r=self.rows()[0]
  self.assertEqual(r["missing_caseversion"],""); self.assertIn("caseversion",json.loads(r["missing_expected_columns"])); self.assertTrue(json.loads(r["warnings"]))
 def test_failure_continues_and_nonzero(self):
  self.write_zip(reac=None); self.assertEqual(p.main(self.args()),1); rows=self.rows(); self.assertEqual(len(rows),3); self.assertEqual([r["profile_status"] for r in rows],["profiled","profiled","failed"])
 def test_sample_marked_and_separate_default(self):
  old=p.DEFAULT_OUT_DIR; p.DEFAULT_OUT_DIR=self.root/"default"
  try:
   args=["--years","2024","--quarters","1","--manifest",str(self.man),"--raw-dir",str(self.raw),"--max-rows-per-table","1"]
   self.assertEqual(p.main(args),0); self.assertTrue((p.DEFAULT_OUT_DIR/"sample/raw_schema_profile.csv").exists()); self.assertFalse((p.DEFAULT_OUT_DIR/"raw_schema_profile.csv").exists())
  finally:p.DEFAULT_OUT_DIR=old
 def test_drift_summary(self):
  base=p.blank(__import__('download_faers').Source(2024,1,"FAERS","ASCII","https://x","q1.zip","q1.zip","d","r","x"),"DEMO",Path("q"),None)
  a=base|{"profile_status":"profiled","normalized_column_names":["a"],"raw_column_names":["A"],"schema_signature":p.signature(["a"]),"missing_expected_columns":[]}
  b=dict(a,quarter=2,normalized_column_names=["a","b"],raw_column_names=["a","B"],schema_signature=p.signature(["a","b"]))
  _,drift,first,forms,_=p.schema_analysis([a,b]); self.assertEqual(drift[0]["columns_added"],["b"]); self.assertEqual(first["DEMO.b"],"2024Q2"); self.assertEqual(forms["DEMO"]["a"],["A","a"])
 def test_persisted_paths_are_portable_and_loading_still_works(self):
  self.assertEqual(p.repository_path(ROOT/"config/faers_sources.csv"),"config/faers_sources.csv")
  self.assertEqual(p.repository_path(self.root/"private/q1.zip"),"q1.zip")
  self.assertEqual(p.main(self.args()),0)
  rows=self.rows()
  with (self.out/"raw_schema_summary.json").open(encoding="utf-8") as f: summary=json.load(f)
  paths=[r["zip_path"] for r in rows]+[summary[k] for k in ("source_manifest","raw_directory","output_profile_csv")]
  for value in paths:
   self.assertNotRegex(value,r"^[A-Za-z]:"); self.assertNotIn("\\",value); self.assertNotIn(str(self.root),value)
  self.assertEqual([r["row_count"] for r in rows],["3","3","3"])
 def test_invocation_from_different_working_directory_is_stable(self):
  command=[sys.executable,str(ROOT/"src/profile_faers_raw.py"),*self.args()]
  first=subprocess.run(command,cwd=ROOT,check=False,capture_output=True,text=True)
  self.assertEqual(first.returncode,0,first.stderr); csv_first=(self.out/"raw_schema_profile.csv").read_text(encoding="utf-8")
  second=subprocess.run(command,cwd=self.root,check=False,capture_output=True,text=True)
  self.assertEqual(second.returncode,0,second.stderr); csv_second=(self.out/"raw_schema_profile.csv").read_text(encoding="utf-8")
  self.assertEqual(csv_first,csv_second); self.assertNotIn("\\",csv_second); self.assertNotRegex(csv_second,r"[A-Za-z]:\\")

if __name__=="__main__":unittest.main()
