import sys
import os
import psycopg2
import shannoni

conn = psycopg2.connect("dbname=TEDDY")
c = conn.cursor()


hlas = eval("{-1:'HLA*Results*Pending', 0:'Not*Eligible', 1:'DR4*030X/0302*DR3*0501/0201', 2:'DR4*030X/0302*DR4*030X/0302', 3:'DR4*030X/0302*DR4*030X/020X', 4:'DR4*030X/0302*DR8*0401/0402', 5:'DR4*030X/0302*DR1*0101/0501', 6:'DR4*030X/0302*DR13*0102/0604', 7:'DR4*030X/0302*DR4*030X/0304', 8:'DR4*030X/0302*DR9*030X/0303', 9:'DR3*0501/0201*DR3*0501/0201', 10:'DR3*0501/0201*DR9*030X/0303', 99:'Results*Under*Review'}")

"""
~/Desktop/TEDDY/t1d_list1/plasma_wgs/case:head sample_mask_microbiome_map.tsv 
7784061	bin0001	Plasma_WGS	343828	107

create table vipie_findings (case_ind text, outcome integer, ia_t1d_list text, tissue text, microbiome text, mask_id text, sample_mask_id text, agedays integer, acc text, genus_desc text, tax_desc text, hits float, primary key(case_ind, mask_id, sample_mask_id, tissue, microbiome, acc));
"""

if (len(sys.argv) < 3):	
	print "need design, tissue, [optional cv, pv]"
	sys.exit()
qexcluded = {}
lowss = []
design = sys.argv[1]
tissue = sys.argv[2]
hit_cutoff = 2
if (tissue == "stool"):
	hit_cutoff = 5

typespecific = "stool_cv"
otypes = "'wgs', 'cv'"
if (tissue == "stool"):
	otypes = "'cv', 'pv'"
alltypes = 0 
if (len(sys.argv) == 4):
	otypes = "'" + sys.argv[3] + "'"
	typespecific = sys.argv[3] + "_"
        tissue_otype = tissue + "_" + sys.argv[3]
else:
	typespecific = ""
	alltypes = 1	
	#print "Need to fix harmonizing for cv and pv, and right now stool tissue only"
	#sys.exit()

tsne = open(typespecific + design + "_" + tissue + "_numeric_viral_exposureds_harm.csv", "w")
outf = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_harm.csv", "w")
adj_paired = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_adj_paired_harm.csv", "w")
sample_writer = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_samples_harm.csv", "w")
itable = "mp183_case_cntrl_t1d_list1"
if (design.find("ia") != -1):
	itable = "mp183_case_cntrl_ia_list1"
_ind_map = {}

#stool cv
cv_case_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/stool_cv/case/sample_mask_microbiome_map_harmony.tsv", "r")
cv_control_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/stool_cv/control/sample_mask_microbiome_map_harmony.tsv", "r")
harmonized_control = {}
harmonized_case = {}
harmonized_control["stool_cv"] = {}
harmonized_case["stool_cv"] = {}

for l in cv_control_sample_f.readlines():
	tk = l.strip().split("\t")
	if (tk[0].find("sample") == -1):
		ind = tk[7]
		sample_mask = tk[0]
		if (harmonized_control["stool_cv"].get(ind) == None):
			harmonized_control["stool_cv"][ind] = "'%s'" %sample_mask
		else:
			harmonized_control["stool_cv"][ind] = harmonized_control["stool_cv"][ind] + ", '%s'" %sample_mask	
for l in cv_case_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_case["stool_cv"].get(ind) == None):
                        harmonized_case["stool_cv"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_case["stool_cv"][ind] = harmonized_case["stool_cv"][ind] + ", '%s'" %sample_mask
cv_case_sample_f.close()
cv_control_sample_f.close()

#repeat stool pv
pv_case_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/stool_pv/case/sample_mask_microbiome_map_harmony.tsv", "r")
pv_control_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/stool_pv/control/sample_mask_microbiome_map_harmony.tsv", "r")
harmonized_control["stool_pv"] = {}
harmonized_case["stool_pv"] = {}
for l in pv_control_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_control["stool_pv"].get(ind) == None):
                        harmonized_control["stool_pv"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_control["stool_pv"][ind] = harmonized_control["stool_pv"][ind] + ", '%s'" %sample_mask   

for l in pv_case_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_case["stool_pv"].get(ind) == None):
                        harmonized_case["stool_pv"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_case["stool_pv"][ind] = harmonized_case["stool_pv"][ind] + ", '%s'" %sample_mask

pv_case_sample_f.close()
pv_control_sample_f.close()

#plasma cv
plasma_cv_case_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/plasma_cv/case/sample_mask_microbiome_map_harmony.tsv", "r")
plasma_cv_control_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/plasma_cv/control/sample_mask_microbiome_map_harmony.tsv", "r")
harmonized_control["plasma_cv"] = {}
harmonized_case["plasma_cv"] = {}
for l in plasma_cv_control_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_control["plasma_cv"].get(ind) == None):
                        harmonized_control["plasma_cv"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_control["plasma_cv"][ind] = harmonized_control["plasma_cv"][ind] + ", '%s'" %sample_mask

for l in plasma_cv_case_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_case["plasma_cv"].get(ind) == None):
                        harmonized_case["plasma_cv"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_case["plasma_cv"][ind] = harmonized_case["plasma_cv"][ind] + ", '%s'" %sample_mask
plasma_cv_case_sample_f.close()
plasma_cv_control_sample_f.close()

plasma_wgs_case_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/plasma_wgs/case/sample_mask_microbiome_map_harmony.tsv", "r")
plasma_wgs_control_sample_f = open("/Users/jakelin/Desktop/TEDDY/ia_list1/plasma_wgs/control/sample_mask_microbiome_map_harmony.tsv", "r")
harmonized_control["plasma_wgs"] = {}
harmonized_case["plasma_wgs"] = {}
for l in plasma_wgs_control_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_control["plasma_wgs"].get(ind) == None):
                        harmonized_control["plasma_wgs"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_control["plasma_wgs"][ind] = harmonized_control["plasma_wgs"][ind] + ", '%s'" %sample_mask

for l in plasma_wgs_case_sample_f.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                sample_mask = tk[0]
                if (harmonized_case["plasma_wgs"].get(ind) == None):
                        harmonized_case["plasma_wgs"][ind] = "'%s'" %sample_mask
                else:
                        harmonized_case["plasma_wgs"][ind] = harmonized_case["plasma_wgs"][ind] + ", '%s'" %sample_mask
plasma_wgs_case_sample_f.close()
plasma_wgs_control_sample_f.close()


expairs = 0
c.execute("select mask_id, case_ind, outcome from " + itable + " order by case_ind, outcome" )
srows = c.fetchall()
header = "mask_id,case_ind,t1d_outcome,ia_outcome,hla_category,gender,fdr,ia_first_diag_agedys,t1d_diag_agedys,delta_ia_t1d,persist_conf_gad_agedys,persist_conf_ia2a_agedys,persist_conf_miaa_agedys,enterovirus_priorx,enterovirus_year1,coxsackiea_priorx,coxsackiea_year1,coxsackieb_priorx,coxsackieb_year1,coxsackieb3_priorx,coxsackieb3_year1,coxsackie_notb3_priorx,coxsackie_notb3_year1,cytomegalovirus_priorx,cytomegalovirus_year1,norovirus_priorx,norovirus_year1,crassphage_priorx,crassphage_year1,parvovirus_priorx,parvovirus_year1,rotavirus_priorx,rotavirus_year1,rhinovirus_priorx,rhinovirus_year1,shannon_year1,shannon_overall,ia_pattern,country,celiac_outcome,samples\n"
outf.write(header)
adj_paired.write(header)
#adj_sample.write(header)

multiple_ev = []
multiple_cmv = []
multiple_noro = []
multiple_crass = []
multiple_parvo = []
multiple_rota = []



"""
_cmv_case = 0
_cmv_control = 0
_ev_case = 0
_ev_control = 0
_noro_case = 0
_noro_control = 0
_crass_case = 0
_crass_control = 0
_parvo_case = 0
_parvo_control = 0
_rota_case = 0
_rota_control = 0
"""
_virus_case = {}
_virus_control = {}

_virus_case["ev"] = 0
_virus_control["ev"] = 0
_virus_case["cmv"] = 0
_virus_control["cmv"] = 0
_virus_case["noro"] = 0
_virus_control["noro"] = 0
_virus_case["crass"] = 0
_virus_control["crass"] = 0
_virus_case["rota"] = 0
_virus_control["rota"] = 0
_virus_case["parvo"] = 0
_virus_control["parvo"] = 0
_virus_case["coxa"] = 0
_virus_control["coxa"] = 0
_virus_case["coxb"] = 0
_virus_control["coxb"] = 0
_virus_case["coxb3"] = 0
_virus_control["coxb3"] = 0
_virus_case["coxb3n"] = 0
_virus_control["coxb3n"] = 0
_virus_case["rhino"] = 0
_virus_control["rhino"] = 0


case_f = {}
control_f = {}
c.execute("select case_ind, outcome, num_files from microbiome_files_qc where tissue = '" + tissue + "' and otype in (" + otypes + ") order by case_ind")
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
        case_control_ratio[ind] = case_microbiome*1.0/(control_microbiome)


#given case_ind, return sample_mask ids as an in clause str
def get_harmonized_sample_mask(tissue_otype, case_ind, outcome):
	if (harmonized_control[tissue_otype].get(case_ind) == None or harmonized_case[tissue_otype].get(case_ind) == None):
		return ""
	hsm = harmonized_control[tissue_otype][case_ind]
	if (outcome == "1"):
		hsm = harmonized_case[tissue_otype][case_ind]
	#samples = str(len(hsm.split(",")))
	return hsm

def parse_rows(_rows, _virus, outcome, mask_id, case_ind, out):
	_prior = 0;_1year = 0;_case=0;_control=0;
	if (len(_rows) > 0):
                _fev = 1
                _age0_exposed = int(_rows[0][0])
                _age0_hits = int(_rows[0][1])
                _evprior = 0
                if (_age0_exposed < 365 and _age0_hits > 1):
                        _1year = 1#_age0_exposed
                _se = {}
                for _e in _rows:
                        _hits = int(_e[1])
                        if (_hits < hit_cutoff):
                                continue
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
			out.write(",".join([case_ind, mask_id, str(outcome), smid, str(_e[0]), str(_hits), _virus]) + "\n")
                _prior = len(_se)
                if (outcome == "1"):
                	_virus_case[_virus] = _virus_case[_virus] + 1#         	_ev_case += 1
                else:
                        _virus_control[_virus] = _virus_control[_virus] + 1
	return (_1year, _prior)

for r in srows:
	_ind_map[r[0]] = r[1]
	outcome = r[2]
	case_ind = r[1]
	mask_id = r[0]	
	
	c.execute("select hla_category, fdr, t1d_diag_agedys, persist_conf_gad_agedys, persist_conf_ia2a_agedys, persist_conf_miaa_agedys, t1d, sex, country, celiac_disease from mp183_onetime where mask_id = '%s'" %(mask_id) )
	srow = c.fetchone()
	ss = [str(s) for s in srow]
	hla = hlas[int(srow[0])]
	fdr = "0"
        ev_cv = "0"
	t1d_outcome = "0"
	celiac_outcome = "0"
	if (srow[-1] == "TRUE"):
		celiac_outcome = "1"
	delta_ia_t1d = 0
	ia_pattern = "control"
	gad_age = srow[3]
	if (gad_age == 0):
		gad_age = 1496 + 500	
        ia2a_age = srow[4]
	if (ia2a_age == 0):
                ia2a_age = 1496 + 500
        miaa_age = srow[5]
	if (miaa_age == 0):
                miaa_age = 1496 + 500
	if (outcome == "1"):
		delta_ia_t1d = 1496*2
		min_ia_age = min(gad_age, ia2a_age, miaa_age)
		if (min_ia_age == gad_age):
			ia_pattern = "GAD"
		if (min_ia_age == ia2a_age):
                        ia_pattern = "IA2A"
		if (min_ia_age == miaa_age):
                        ia_pattern = "mIAA"

        ia_first_age = min(gad_age, ia2a_age, miaa_age)#srow[3],srow[4],srow[5])
	if (outcome == "0"):
		ia_first_age = 1496 + 500
		delta_ia_t1d = 1496*2
	if (srow[6] == "TRUE"):
                t1d_outcome = "1"
		delta_ia_t1d = srow[2] - ia_first_age
	if (srow[1] == "TRUE"):
		fdr = "1"
	t1d_age = srow[2]
	if (srow[2] == 0):
		t1d_age = 1496*2
        gender = srow[7]
        cc = srow[8]
        #1=US, 2=FIN, 3=GER, 4=SWE
        country = "US"
        if (cc == "2"):
            country = "FIN"
        if (cc == "3"):
            country = "GER"
        if (cc == "4"):
            country = "SWE" 
	#delta_ia_t1d = 9999
	#delta_ia_t1d = srow[2] - ia_first_age
	key_days = [str(d) for d in [ia_first_age, t1d_age, gad_age, ia2a_age, miaa_age, delta_ia_t1d]]
	_fev = 0
	_hsm = ""
	if (alltypes == 1):
		if (tissue == "stool"):
			_hsm_cv = get_harmonized_sample_mask("stool_cv", case_ind, outcome)
			if (len(_hsm_cv) < 3):
				_hsm_cv = ""			
			_hsm_pv = get_harmonized_sample_mask("stool_pv", case_ind, outcome) 
			if (len(_hsm_pv) < 3):
                        	_hsm_pv = ""
		if (tissue == "plasma"):
			_hsm_cv = get_harmonized_sample_mask("plasma_cv", case_ind, outcome)
                	_hsm_pv = get_harmonized_sample_mask("plasma_wgs", case_ind, outcome)
		_hsm = _hsm_cv + "," + _hsm_pv
	else:
		_hsm = get_harmonized_sample_mask(tissue_otype, case_ind, outcome)
	if (_hsm[0:1] == ","):
		_hsm = _hsm[1:]
	if (_hsm[-1:] == ","):
                _hsm = _hsm[:-1]
	_in_hsm = "and sample_mask_id in (" + _hsm + ")"
	if (_hsm == "" or _hsm == ","):#.find("NONE") != -1):
		expairs += 1
		qexcluded[case_ind] = 1
		print "exclude no samples in ", case_ind
		continue
	len_hsm = str(len(_hsm.split("',")))
	if (int(len_hsm) < 4):
		lowss.append(case_ind)
	print tissue, otypes, " case_ind, outcome ", case_ind, outcome, _in_hsm, len_hsm
	enterovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%enterovirus%' and lower(genus_desc) not like '%rhinovirus%') " + _in_hsm + " order by agedays asc, hits desc"
	c.execute(enterovirus_q)
	enterovirus_rows = c.fetchall()
	_found = len(enterovirus_rows)
	res = parse_rows(enterovirus_rows,"ev",outcome, mask_id, case_ind, sample_writer)	
	_evprior_1year = res[0]
	_evprior = res[1]
	
        #_fcoxb = 0	
        coxb_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(coxb_q)
        coxb_rows = c.fetchall()
        res = parse_rows(coxb_rows,"coxb",outcome, mask_id, case_ind, sample_writer)
        _coxbprior_1year = res[0]
        _coxbprior = res[1]

	coxa_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie a%' or lower(genus_desc) like '%coxsackievirus a%') " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(coxa_q)
        coxa_rows = c.fetchall()
        res = parse_rows(coxa_rows,"coxa",outcome, mask_id, case_ind, sample_writer)
        _coxaprior_1year = res[0]
        _coxaprior = res[1]

        #_fcoxb3 = 0
        coxb3_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b3%' or lower(genus_desc) like '%coxsackievirus b3%') " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(coxb3_q)
        coxb3_rows = c.fetchall()
        res = parse_rows(coxb3_rows,"coxb3",outcome, mask_id, case_ind, sample_writer)
        _coxb3prior_1year = res[0]
        _coxb3prior = res[1]

        coxb3n_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') and lower(genus_desc) not like '%coxsackie b3%' and lower(genus_desc) not like '%coxsackievirus b3%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(coxb3n_q)
        coxb3n_rows = c.fetchall()
        res = parse_rows(coxb3n_rows,"coxb3n",outcome, mask_id, case_ind, sample_writer)
        _coxb3nprior_1year = res[0]
        _coxb3nprior = res[1]

        _fcmv = 0
	cytomegalovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%cytomegalo%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(cytomegalovirus_q)
        cmvvirus_rows = c.fetchall()
        res = parse_rows(cmvvirus_rows,"cmv",outcome, mask_id, case_ind, sample_writer)
        _cmvprior_1year = res[0];_cmvprior = res[1]
	#;_cmv_case = res[2];_cmv_control = res[3]

	#rhinovirus - most common, common cold
	rhinovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rhinovirus%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(rhinovirus_q)
        rhinovirus_rows = c.fetchall()
        res = parse_rows(rhinovirus_rows,"rhino",outcome, mask_id, case_ind, sample_writer)
        _rhinoprior_1year = res[0];_rhinoprior = res[1]

	norovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in ( " + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%norovirus%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(norovirus_q)
        norovirus_rows = c.fetchall()
	res = parse_rows(norovirus_rows,"noro",outcome, mask_id, case_ind, sample_writer)
        _noroprior_1year = res[0];_noroprior = res[1]
	#;_noro_case = res[2];_noro_control = res[3]
	crassphage_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%crassphage%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(crassphage_q)
        crass_rows = c.fetchall()
	res = parse_rows(crass_rows,"crass",outcome, mask_id, case_ind, sample_writer)
        _crassprior_1year = res[0];_crassprior = res[1]
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%parvo%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(_q)
        parvo_rows = c.fetchall()
        res = parse_rows(parvo_rows,"parvo",outcome, mask_id, case_ind, sample_writer)
        _parvoprior_1year = res[0];_parvoprior = res[1]
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rotavirus%' " + _in_hsm + " order by agedays asc, hits desc"
        c.execute(_q)
	#print(_q)
        rota_rows = c.fetchall()
        res = parse_rows(rota_rows,"rota",outcome, mask_id, case_ind, sample_writer)
        _rotaprior_1year = res[0];_rotaprior = res[1]

	_q = "select avg(shannon) from vipie_findings_sample_mask_summary where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and agedays < 365 "  + _in_hsm 
        c.execute(_q)
        _row = c.fetchone()
	year1shannon = _row[0]
	if (year1shannon == None):
                year1shannon = 0
	_q = "select avg(shannon) from vipie_findings_sample_mask_summary where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "'" + _in_hsm 
        c.execute(_q)
        _row = c.fetchone()
        allshannon = _row[0]
	
	if (allshannon == None):
		allshannon = 0
	#if (case_control_ratio.get(case_ind) ==  None or (case_control_ratio.get(case_ind) < .8) or (case_control_ratio.get(case_ind) > 1.2)):#) and (control_f[case_ind] < 5)):
	#	qexcluded[case_ind] = 1
	#	expairs += 1
	#	continue
	case2control = 1.0#case_control_ratio[case_ind]
	#if (outcome == "0"):
	#	case2control = 1 - case2control
	virus_corrected = [str(x/case2control) for x in [_evprior, _coxaprior, _coxbprior, _coxb3prior, _coxb3nprior, _cmvprior, _noroprior, _crassprior, _parvoprior, _rotaprior, _rhinoprior]]	
	
	line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], str(_evprior), str(_evprior_1year), str(_coxaprior), str(_coxaprior_1year), str(_coxbprior), str(_coxbprior_1year), str(_coxb3prior), str(_coxb3prior_1year), str(_coxb3nprior), str(_coxb3nprior_1year),str(_cmvprior), str(_cmvprior_1year), str(_noroprior), str(_noroprior_1year), str(_crassprior), str(_crassprior_1year), str(_parvoprior), str(_parvoprior_1year),str(_rotaprior), str(_rotaprior_1year), str(_rhinoprior), str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome, len_hsm])

		
	outf.write(line + "\n")

        line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], virus_corrected[0], str(_evprior_1year), virus_corrected[1], str(_coxaprior_1year), virus_corrected[2], str(_coxbprior_1year), virus_corrected[3], str(_coxb3prior_1year), virus_corrected[4], str(_coxb3nprior_1year), virus_corrected[5], str(_cmvprior_1year), virus_corrected[6], str(_noroprior_1year), virus_corrected[7], str(_crassprior_1year), virus_corrected[8], str(_parvoprior_1year), virus_corrected[9], str(_rotaprior_1year), virus_corrected[10], str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome, len_hsm])
	adj_paired.write(line + "\n")
	#numfiles = 0
	#if (outcome == "1"):
	#	numfiles = case_f[case_ind]
	#if (outcome == "0"):
        #        numfiles = control_f[case_ind]
        #virus_corrected = [str(x/(numfiles*1.0)) for x in [_evprior, _coxaprior, _coxbprior, _coxb3prior, _coxb3nprior, _cmvprior, _noroprior, _crassprior, _parvoprior, _rotaprior, _rhinoprior]]
	#line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], virus_corrected[0], str(_evprior_1year), virus_corrected[1], str(_coxaprior_1year), virus_corrected[2], str(_coxbprior_1year), virus_corrected[3], str(_coxb3prior_1year), virus_corrected[4], str(_coxb3nprior_1year), virus_corrected[5], str(_cmvprior_1year), virus_corrected[6], str(_noroprior_1year), virus_corrected[7], str(_crassprior_1year), virus_corrected[8], str(_parvoprior_1year), virus_corrected[9], str(_rotaprior_1year), virus_corrected[10], str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome])
        #adj_sample.write(line + "\n")

	tsne.write(" ".join([str(_evprior), str(_coxaprior), str(_coxbprior), str(_coxb3prior), str(_coxb3nprior), str(_cmvprior), str(_noroprior), str(_crassprior), str(_parvoprior),str(_rotaprior), str(_rhinoprior), str(year1shannon), str(allshannon)]) + "\n")

conn.close()
print "case", _virus_case
print "control", _virus_control
sample_writer.close()
outf.close()
tsne.close()
adj_paired.close() 
#adj_sample.close() 
print "excluded pairs", expairs/2
print qexcluded
print "low_sample_subjects", len(lowss)/2
print lowss
