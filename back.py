import subprocess
import sys
protein_file=sys.argv[1]
cutoff_num_contact=sys.argv[2]
cutoff_perc_contact=sys.argv[3]
id_job=sys.argv[4]
email=sys.argv[5]
subprocess.Popen( ["python lig_infer.py "+protein_file+" "+cutoff_num_contact+" "+cutoff_perc_contact+" "+id_job+" "+email], shell=True)
