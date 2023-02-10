import csv
import numpy
"""
head annotation_pairedMax_all_data.txt 
"Status"	"Time"	"Number_IA"	"T1D"	"Endotype"
"7014485"	"Case"	3	1	1	"IAA 1st"
"1782944"	"Case"	3	1	0	"IAA 1st"
"3815259"	"Case"	3	3	0	"IAA 1st"
"1362053"	"Case"	3	2	0	"Other"
"9936260"	"Case"	3	1	0	"IAA 1st"
"5343323"	"Case"	3	1	0	"IAA 1st"
"2215043"	"Case"	3	3	0	"IAA 1st"
"9611347"	"Case"	3	2	0	"IAA 1st"
"5695346"	"Case"	3	1	1	"IAA 1st"
(base) lljali@wks-100210-mac TEDDY_HII % grep 7014485 Mp183_sample_mapping_v2.csv 
(base) lljali@wks-100210-mac TEDDY_HII % grep 7014485 MP183_RNASeq_Sample_Mapping_File_Masked_v1.csv 
590368;7014485;12;11;380;RNA;RNA-Seq
"""

ncc = {}

"""
(base) lljali@wks-100210-mac paired_max % wc -l *
     159 ia_rnaseq_12mon_paired.csv
     140 ia_rnaseq_15mon_paired.csv
     112 ia_rnaseq_18mon_paired.csv
      99 ia_rnaseq_21mon_paired.csv
      80 ia_rnaseq_24mon_paired.csv
      64 ia_rnaseq_27mon_paired.csv
      48 ia_rnaseq_30mon_paired.csv
      47 ia_rnaseq_33mon_paired.csv
      30 ia_rnaseq_36mon_paired.csv
     209 ia_rnaseq_3mon_paired.csv
     209 ia_rnaseq_6mon_paired.csv
     175 ia_rnaseq_9mon_paired.csv
    1372 total
(base) lljali@wks-100210-mac paired_max % vi ia_rnaseq_3mon_paired.csv
(base) lljali@wks-100210-mac paired_max % pwd
/Users/lljali/data/TEDDY/TEDDY_HII/paired_max
(base) lljali@wks-100210-mac paired_max % head ia_rnaseq_3mon_paired.csv
,case_ind,age_IA,case_mask_id,case_draw_dte_agedys,case_sample_id,control_mask_id,control_draw_dte_agedys,control_sample_id,due_num
,344,15,590368,380,7014485,642022,381,8592804,12
,345,35,583295,1003,1782944,848030,985,6076356,33
,346,27,682042,730,3815259,711669,755,7546875,24
,347,28,629627,743,1362053,585002,728,2253606,24
,340,6,231441,126,9936260,651328,133,7767514,3
,341,17,991505,451,5343323,504610,457,5114656,15
,343,15,268779,369,2215043,632318,359,9531389,12
,348,15,844268,351,9611347,491018,368,7867781,12
,349,12,954008,301,5695346,349394,314,1603739,9
(base) lljali@wks-100210-mac paired_max % 
"""

ncc_tracked = {}
for i in [3,6,9,12]:
	f = "./paired_max/ia_rnaseq_" + str(i) +"mon_paired.csv"
	fo = open(f, "r")
	#fo.close()
	reader = csv.DictReader(fo, delimiter=",")
	for r in reader:
		case_ind = r["case_ind"] 
		if (ncc_tracked.get(case_ind) == None):
			ncc_tracked[case_ind] = r["age_IA"]
			ncc[r["case_mask_id"]] = int(r["age_IA"])

rna_sample_subject = {}

import numpy
print("ncc age_IA", len(ncc_tracked))
#, numpy.median(ncc.values())
ages = []
for s in ncc:
	ages.append(ncc[s])

print(numpy.median(ages), "median")

f = open("MP183_RNASeq_Sample_Mapping_File_Masked_v1.csv", "r")
for l in f.readlines():
	tk = l.strip().split(";")
	dtype = tk[-1]
	if (dtype == "RNA-Seq"):
		rna_sample_subject[tk[1]] = tk[0]

f.close()
print("N Samples", len(rna_sample_subject))

"""
(base) lljali@wks-100210-mac TEDDY_HII % head ../gad_only.csv
"","case_ind","mask","ia","ia_age","iaa","gada","ia2a","t1d"
"1",52,892347,1,1403,FALSE,TRUE,FALSE,FALSE
"98",197,403715,1,1517,FALSE,TRUE,FALSE,FALSE
(base) lljali@wks-100210-mac TEDDY_HII % head ../iaa_only.csv
"","case_ind","mask","ia","ia_age","iaa","gada","ia2a","t1d"
"58",325,587714,1,363,TRUE,FALSE,FALSE,FALSE
"202",402,478935,1,460,TRUE,FALSE,FALSE,FALSE
(base) lljali@wks-100210-mac TEDDY_HII % head ../2plus_ia.csv 
"","case_ind","mask","ia","ia_age","iaa","gada","ia2a","t1d"
"16",27,451065,1,1211,TRUE,TRUE,TRUE,FALSE
"""

mask_seroage = {}
seroages = []
for fn in ["../gad_only.csv", "../iaa_only.csv", "../2plus_ia.csv"]:
	fo = open(fn, "r")
	reader = csv.DictReader(fo, delimiter=",")
	for r in reader:
		mask_seroage[r["mask"]] = int(r["ia_age"])
		seroages.append(int(r["ia_age"]))
	fo.close()
print("N mask_Seroage", len(mask_seroage), "median age", numpy.median(seroages))

proc = {}
sero = {}
sero["IAA 1st"] = []
sero["Other"] = []
sero["GADA 1st"] = []
err = {}
f = open("annotation_pairedMax_all_data.txt", "r")
for l in f.readlines():
	l = l.replace('"', '')
	tk = l.strip().split("\t")
	sample = tk[0]
	if (rna_sample_subject.get(sample) != None):
		subject = rna_sample_subject[sample]
		if (proc.get(subject) == None):
			proc[subject] = tk[-1]
			if (mask_seroage.get(subject) == None):
				err[subject] = 1
				continue
			onsetage = mask_seroage[subject]
			st = tk[-1]
			sero[st].append(onsetage)
f.close()
			
print("N Proc subjects",len(proc), "err", len(err))
#print(len(sero["IAA 1st"] ), "IAA 1st", numpy.median(sero["IAA 1st"]))
q75, q25 = numpy.percentile(sero["IAA 1st"], [75 ,25])
print(len(sero["IAA 1st"] ), "IAA 1st", numpy.median(sero["IAA 1st"]), q25, q75)

#print(q25,)
q75, q25 = numpy.percentile(sero["GADA 1st"], [75 ,25])
print(len(sero["GADA 1st"] ), "GADA 1st", numpy.median(sero["GADA 1st"]), q25, q75)
q75, q25 = numpy.percentile(sero["Other"], [75 ,25])
print(len(sero["Other"] ), "Other", numpy.median(sero["Other"]), q25, q75)

