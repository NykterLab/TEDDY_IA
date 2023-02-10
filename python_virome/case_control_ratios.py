import sys
import os
import psycopg2
import shannoni

conn = psycopg2.connect("dbname=TEDDY")
c = conn.cursor()


case_f = {}
control_f = {}
c.execute("select case_ind, outcome, num_files from microbiome_files where tissue = 'stool' and otype = 'cv' order by case_ind")

rows = c.fetchall()
for r in rows:
	ind = r[0]
	outcome = r[1]
	numfiles = int(r[2])
	if (outcome == 1):
		case_f[ind] = numfiles
	if (outcome == 0):
		control_f[ind] = numfiles

case_control_ratio = {}
balanced = 0
for ind in case_f:
	if (control_f.get(ind) == None):
		print "case_ind 0 files for control, maybe exlude", ind
		continue
	case_microbiome = case_f[ind] 
	control_microbiome = control_f[ind]
	if (control_microbiome >= case_microbiome):
		balanced += 1
	case_control_ratio[ind] = case_microbiome*1.0/control_microbiome

print len(case_control_ratio), balanced 
print ("10 95 99", case_control_ratio["10"], case_control_ratio["95"], case_control_ratio["99"])

conn.close()
