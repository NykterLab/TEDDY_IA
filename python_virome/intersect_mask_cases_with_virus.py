import csv
"""
case_ind,ia,ia_age,cd,t1d,hla,hla_discrepant,mask,iaa,gada,ia2a,country,sample_mask,due_num,draw_mon,sample_age_days,diversity,cvb_exp,cva_exp,echo_exp,rhino_exp,adenoa_exp,adenoc_exp,adenod_exp,adenob_exp,adeno_exp,ev_exp,crassphage_exp,noro_exp,cytomegalo_exp,parvo_exp,cvb1_exp,cvb2_exp,cvb3_exp,cvb4_exp,cvb5_exp,bifido_phage,papiloma,n_samples
3,1,1645,NA,FALSE,DR4/DR3,N,314979,TRUE,FALSE,FALSE,FIN,9467264,3,1,114,0.191444,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15
"""

#store all pairs
"""
LM8-945-015:mp183 linjake$ head ia_rnaseq_3mon_paired.csv
,case_ind,age_IA,case_mask_id,case_draw_dte_agedys,case_sample_id,control_mask_id,control_draw_dte_agedys,control_sample_id,due_num
,217,1215,935244,1124,9407979,500882,1082,6517442,36
,212,946,617795,825,1152770,734670,825,7587501,27
,210,272,441651,191,4654972,659851,190,7082152,6
,218,748,296606,671,3217433,275258,657,1537843,21
,133,460,268602,391,4166861,788365,402,3887174,12
,137,1395,869968,1325,7845085,328344,1306,7943944,42
,135,1028,603749,948,2619331,691030,924,2613218,30
,134,279,945130,202,4438099,233859,178,9449560,6
,225,379,270387,288,2962359,204528,302,3549419,9
"""

pairedf = ["ia_rnaseq_3mon_paired.csv", "ia_rnaseq_6mon_paired.csv", "ia_rnaseq_9mon_paired.csv", "ia_rnaseq_12mon_paired.csv", "ia_rnaseq_15mon_paired.csv", "ia_rnaseq_18mon_paired.csv", "ia_rnaseq_21mon_paired.csv", "ia_rnaseq_24mon_paired.csv", "ia_rnaseq_27mon_paired.csv"]

case_ctl_rnaseq = {}
for p in pairedf:
	pf = open(p, "r")
	reader = csv.DictReader(pf)
	for r in reader:
		mask_sample = r["case_mask_id"] + ":" + r["case_sample_id"]
		case_ctl_rnaseq[mask_sample] = r["control_mask_id"] + ":" + r["control_sample_id"]
	pf.close()

f = open("./report_plasma_virus_ia_samples_acc.csv", "r")
print(len(case_ctl_rnaseq), "total paired RNASeq files")

reader = csv.DictReader(f)
env_cases = {} 
for r in reader:
	mask = r["mask"]
	ev = int(r["ev_exp"])
	ia = int(r["ia"])
	if (ia == 1 and ev > 0):
		env_cases[mask] = r
f.close()

print("subjects with enterovirus found in plasma", len(env_cases))

f = open("./report_virus_ia_samples_acc.csv", "r")
stool_env_cases = {}
reader = csv.DictReader(f)
for r in reader:
        mask = r["mask"]
        ev = int(r["ev_exp"])
        ia = int(r["ia"])
        if (ia == 1 and ev > 0):
                stool_env_cases[mask] = r
f.close()
print("subjects with enterovirus found in stool", len(stool_env_cases))

plasma_stool_env = set(set(env_cases.keys()) & set(stool_env_cases.keys()))
print(plasma_stool_env)

import sys

fin = sys.argv[1]
"""
"","case__mask_id","case_sample_id_0","case_sample_id_1"
"1","296606","3217433","9368401"
"""
fo = open(fin, "r")
foo = open(fin.split(".")[0] + "_virus_sample_masks.csv", "w")
foo.write("case_mask,case_sample_id_0,case_sample_id_1,control_mask,control_sample_id_0,control_sample_id_1\n")

reader = csv.DictReader(fo, quotechar='"')
pi = 0
for r in reader:
	smask= (r["case__mask_id"])
	if (stool_env_cases.get(smask) != None or env_cases.get(mask) != None):
		pi+=1
		control_sample0 = case_ctl_rnaseq[smask + ":" + r["case_sample_id_0"]].split(":")
		control_sample1 = case_ctl_rnaseq[smask + ":" + r["case_sample_id_1"]].split(":")
		print (" threshold found with enterovirus", smask, r["case_sample_id_0"],r["case_sample_id_1"], "matched controls", control_sample0[0], control_sample0[1],control_sample1[1])
		foo.write(",".join([ smask, r["case_sample_id_0"],r["case_sample_id_1"],control_sample0[0], control_sample0[1],control_sample1[1]]) + "\n")
fo.close()
outn = foo.name
foo.close()
print(outn, "stored paired case and control sample masks", pi)

