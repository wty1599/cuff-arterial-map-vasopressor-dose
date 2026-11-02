"""Synthetic regression checks for the release boundary. No clinical inputs."""
from pathlib import Path
import csv,io,json,subprocess,sys,tempfile,unittest,zipfile
import small_cell_gate as gate
def csvbytes(header,rows,delimiter=','):
 s=io.StringIO(newline='');w=csv.writer(s,delimiter=delimiter);w.writerow(header);w.writerows(rows);return s.getvalue().encode()
class ReleaseGateTests(unittest.TestCase):
 def blocked(self,b,name):self.assertTrue(gate.scan_bytes(b,name)[0])
 def allowed(self,b,name):self.assertFalse(gate.scan_bytes(b,name)[0])
 def test_patient_threshold(self):self.blocked(csvbytes(['n_pairs','n_patients','mean_map'],[[40,3,70]]),'patients.csv')
 def test_pair_threshold(self):self.blocked(csvbytes(['n_pairs','difference'],[[4,8]]),'pairs.csv')
 def test_boundary_allowed(self):self.allowed(csvbytes(['n_pairs','n_patients','mean_map'],[[10,10,70]]),'safe.csv')
 def test_missing_measurement_is_not_released_as_a_value(self):self.allowed(csvbytes(['n_pairs','estimate'],[['<10','suppressed']]),'masked.csv')
 def test_masked_count_retaining_measurement_is_blocked(self):self.blocked(csvbytes(['n_pairs','estimate'],[['<10',.2]]),'masked_bad.csv')
 def test_counts_without_measurements(self):self.allowed(csvbytes(['file','count'],[['a',2]]),'inventory.csv')
 def test_count_aliases(self):
  for field in ['frequency','frequency_n','sample_size','n','pairs','patients','denominator_pairs','cell_pairs']:
   with self.subTest(field=field):self.blocked(csvbytes([field,'reading'],[[3,70]]),'alias.csv')
 def test_patient_identifiers(self):
  for field in ['patient_id','subject_id','patientunitstayid','pid','sid','cid','row_id']:
   with self.subTest(field=field):self.blocked(csvbytes([field,'value'],[['synthetic',70]]),'id.csv')
 def test_time_aliases(self):
  for field in ['timestamp','time','charttime','cuff_time','offset']:
   with self.subTest(field=field):self.blocked(csvbytes([field,'mean_map'],[['synthetic',70]]),'time.csv')
 def test_clock_schema_parameter(self):self.allowed(json.dumps({'mimic_observation_time':'charttime without storetime','adult_age':18}).encode(),'params.json')
 def test_json_records(self):self.blocked(json.dumps([{'n_pairs':3,'difference':5}]).encode(),'records.json')
 def test_json_vectors(self):self.blocked(json.dumps({'n_pairs':[3,20],'mean_map':[70,71]}).encode(),'vectors.json')
 def test_json_columns_matrix(self):self.blocked(json.dumps({'columns':['n_pairs','mean_map'],'data':[[3,70]]}).encode(),'matrix.json')
 def test_numeric_dynamic_keys(self):
  b=json.dumps({'314159265':{'n_pairs':3,'mean_map':70}}).encode();f,e=gate.scan_bytes(b,'keys.json');self.assertTrue(f);self.assertNotIn('314159265',json.dumps([f,e]))
 def test_numeric_root_vector(self):self.blocked(b'[1,2,3]','vector.json')
 def test_numeric_root_matrix(self):self.blocked(b'[[1,2],[3,4]]','matrix.json')
 def test_numeric_header(self):self.blocked(b'1,2\n3,4\n','no_header.csv')
 def test_duplicate_csv_header(self):self.blocked(b'n_pairs,n_pairs,mean_map\n3,20,70\n','duplicate.csv')
 def test_ragged_extra_field(self):self.blocked(b'n_pairs,mean_map\n20,70,3\n','ragged.csv')
 def test_ragged_missing_field(self):self.blocked(b'n_pairs,mean_map\n20\n','ragged.csv')
 def test_duplicate_json_key(self):self.blocked(b'{"n_pairs":3,"n_pairs":20,"mean_map":70}','duplicate.json')
 def test_delimiter_renaming(self):self.blocked(b'n_pairs,mean_map\n3,70\n','renamed.tsv')
 def test_utf16(self):self.blocked('n_pairs,mean_map\n3,70\n'.encode('utf-16'),'utf16.csv')
 def test_json_masquerading_as_csv(self):self.blocked(b'{"n_pairs":3,"mean_map":70}','renamed.csv')
 def test_text_masquerading(self):self.blocked(b'n_pairs,mean_map\n3,70\n','renamed.txt')
 def test_text_json_masquerading(self):self.blocked(b'{"n_pairs":3,"mean_map":70}','renamed.txt')
 def test_zip_recursion(self):
  b=io.BytesIO()
  with zipfile.ZipFile(b,'w') as z:z.writestr('data/unsafe.csv',b'n_pairs,mean_map\n3,70\n')
  self.blocked(b.getvalue(),'payload.zip')
 def test_zip_restricted_member(self):
  b=io.BytesIO()
  with zipfile.ZipFile(b,'w') as z:z.writestr('restricted/readme.md','synthetic')
  self.blocked(b.getvalue(),'payload.zip')
 def test_xlsx(self):
  import openpyxl
  b=io.BytesIO();w=openpyxl.Workbook();w.active.append(['n_pairs','mean_map']);w.active.append([3,70]);w.save(b);self.blocked(b.getvalue(),'unsafe.xlsx')
 def test_parquet(self):
  import pyarrow as pa,pyarrow.parquet as pq
  b=io.BytesIO();pq.write_table(pa.table({'n_pairs':[3],'mean_map':[70.]}),b);self.blocked(b.getvalue(),'unsafe.parquet')
 def test_parse_failure_closed(self):self.blocked(b'not a valid parquet','unreadable.parquet')
 def test_staged_bytes_are_checked(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);subprocess.run(['git','init','-q',str(p)],check=True,capture_output=True)
   f=p/'data.csv';f.write_bytes(b'n_pairs,mean_map\n3,70\n');subprocess.run(['git','-C',str(p),'add','data.csv'],check=True,capture_output=True)
   f.write_bytes(b'n_pairs,mean_map\n20,70\n');self.assertTrue(gate.staged(p)[0])
 def test_output_does_not_echo_values(self):
  raw=b'n_pairs,mean_map\n3,123.456789\n';f,e=gate.scan_bytes(raw,'echo.csv');self.assertNotIn('123.456789',json.dumps([f,e]))
 def test_pre_push_checks_commit_even_if_worktree_is_safe(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);cmd=['git','-C',str(p)];subprocess.run(['git','init','-q',str(p)],check=True,capture_output=True)
   subprocess.run(cmd+['config','user.name','Synthetic Gate Test'],check=True,capture_output=True);subprocess.run(cmd+['config','user.email','synthetic@example.invalid'],check=True,capture_output=True)
   f=p/'data.csv';f.write_bytes(b'n_pairs,mean_map\n3,70\n');subprocess.run(cmd+['add','data.csv'],check=True,capture_output=True);subprocess.run(cmd+['commit','-qm','synthetic'],check=True,capture_output=True)
   oid=subprocess.run(cmd+['rev-parse','HEAD'],check=True,capture_output=True).stdout.decode().strip();f.write_bytes(b'n_pairs,mean_map\n20,70\n')
   self.assertTrue(gate.pre_push(p,f'refs/heads/main {oid} refs/heads/main '+('0'*40))[0])
 def test_cache_keeps_current_locator(self):
  b=b'n_pairs,mean_map\n3,70\n';gate.scan_bytes(b,'first.csv');f,e=gate.scan_bytes(b,'second.csv');self.assertEqual(f[0]['path'],'second.csv')
if __name__=='__main__':unittest.main()
