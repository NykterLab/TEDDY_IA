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

if (len(sys.argv) != 3):	
	print "need design, tissue"
	sys.exit()

design = sys.argv[1]
tissue = sys.argv[2]
hit_cutoff = 2
if (tissue == "stool"):
	hit_cutoff = 5

#otype = sys.argv[3]

outf = open(design + "_" + tissue + "_viral_exposured_agedays_coxb3_ignore.csv", "w")
itable = "mp183_case_cntrl_t1d_list1"
if (design.find("ia") != -1):
	itable = "mp183_case_cntrl_ia_list1"
_ind_map = {}


c.execute("select mask_id, case_ind, outcome from " + itable + " order by case_ind, outcome" )
srows = c.fetchall()
header = "mask_id,case_ind,t1d_outcome,ia_outcome,hla_category,gender,fdr,ia_first_diag_agedys,t1d_diag_agedys,delta_ia_t1d,persist_conf_gad_agedys,persist_conf_ia2a_agedys,persist_conf_miaa_agedys,enterovirus_priorx,enterovirus_year1,coxsackieb_priorx,coxsackieb_year1,coxsackieb3_priorx,coxsackieb3_year1,coxsackie_notb3_priorx,coxsackie_notb3_year1,cytomegalovirus_priorx,cytomegalovirus_year1,norovirus_priorx,norovirus_year1,crassphage_priorx,crassphage_year1,parvovirus_priorx,parvovirus_year1,rotavirus_priorx,rotavirus_year1,rhinovirus_priorx,rhinovirus_year1,ev_cv,shannon_year1,shannon_overall,ia_pattern,country\n"
outf.write(header)

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
_virus_case["coxb"] = 0
_virus_control["coxb"] = 0
_virus_case["coxb3"] = 0
_virus_control["coxb3"] = 0
_virus_case["coxb3n"] = 0
_virus_control["coxb3n"] = 0
_virus_case["rhino"] = 0
_virus_control["rhino"] = 0

_virus_case_overlapped = {}
_virus_control_overlapped = {}

def parse_rows(_rows, _virus, outcome, interval=90):
	_prior = 0;_1year = 0;_case=0;_control=0;
        ages = []
        overlapped = 0
	if (len(_rows) > 0):
                _fev = 1
                _age0_exposed = int(_rows[0][0])
                _age0_hits = int(_rows[0][1])
                _evprior = 0
                if (_age0_exposed < 365 and _age0_hits > 1):
                        _1year = 1#_age0_exposed
                _se = {}
                prev_age = 0
                mask_id = _rows[0][-1]
		overlapped = 0
		current_age = 0
		prev_sample = ""
                for _e in _rows:
                        _hits = int(_e[1])
                        _age = int(_e[0])
                        current_age = _age
			current_sample = _e[3]
                        if (prev_age != 0 and prev_sample != current_sample and prev_sample != ""):
                        	delta_age = current_age - prev_age
                                if (delta_age <= 90):
                                    overlapped += 1
                        if (_hits < hit_cutoff):
                                continue
                        
                        smid = _e[3]
                        ages.append(_age)
                        if (_se.get(smid) == None):
                                _se[smid] = float(_e[1])
                        else:
                                _se[smid] = _se[smid] + float(_e[1])
                        prev_age = current_age
			prev_sample = current_sample
                _prior = len(_se)
                if (outcome == "1"):
                	_virus_case[_virus] = _virus_case[_virus] + 1
			if (overlapped > 0):
				_virus_case_overlapped[mask_id] = overlapped
                else:
                        _virus_control[_virus] = _virus_control[_virus] + 1
			if (overlapped > 0):
				_virus_control_overlapped[mask_id] = overlapped
                
	return (_1year, _prior)

for r in srows:
	_ind_map[r[0]] = r[1]
	outcome = r[2]
	case_ind = r[1]
	mask_id = r[0]	
	
	c.execute("select hla_category, fdr, t1d_diag_agedys, persist_conf_gad_agedys, persist_conf_ia2a_agedys, persist_conf_miaa_agedys, t1d, sex, country from mp183_onetime where mask_id = '%s'" %(mask_id) )
	srow = c.fetchone()
	ss = [str(s) for s in srow]
	hla = hlas[int(srow[0])]
	fdr = "0"
        ev_cv = "0"
	t1d_outcome = "0"
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
        otypes = "'wgs', 'cv'"
	if (tissue == "stool"):
		otypes = "'wgs', 'cv', 'pv'"	
	enterovirus_q = "select agedays, hits, genus_desc, sample_mask_id, mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%enterovirus%' order by agedays asc"
	c.execute(enterovirus_q)
	enterovirus_rows = c.fetchall()
	_found = len(enterovirus_rows)
	#res = parse_rows(enterovirus_rows,"ev",outcome)	
	_evprior_1year = 0#res[0]
	_evprior = 0#res[1]
	
	res = [0,0]
        #_fcoxb = 0	
        coxb_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') order by agedays asc"
        c.execute(coxb_q)
        coxb_rows = c.fetchall()
        #res = parse_rows(coxb_rows,"coxb",outcome)
        _coxbprior_1year = res[0]
        _coxbprior = res[1]

        #_fcoxb3 = 0
        coxb3_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b3%' or lower(genus_desc) like '%coxsackievirus b3%') order by agedays asc"
        c.execute(coxb3_q)
        coxb3_rows = c.fetchall()
        #res = parse_rows(coxb3_rows,"coxb3",outcome)
        _coxb3prior_1year = res[0]
        _coxb3prior = res[1]

        coxb3n_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') and lower(genus_desc) not like '%coxsackie b3%' and lower(genus_desc) not like '%coxsackievirus b3%' order by agedays asc"
        c.execute(coxb3n_q)
        coxb3n_rows = c.fetchall()
        #res = parse_rows(coxb3n_rows,"coxb3n",outcome)
        _coxb3nprior_1year = res[0]
        _coxb3nprior = res[1]

        _fcmv = 0
	cytomegalovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%cytomegalo%' order by agedays asc"
        c.execute(cytomegalovirus_q)
        cmvvirus_rows = c.fetchall()
        res = parse_rows(cmvvirus_rows,"cmv",outcome)
        _cmvprior_1year = res[0];_cmvprior = res[1]

	#rhinovirus - most common, common cold
	rhinovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rhinovirus%' order by agedays asc"
        c.execute(rhinovirus_q)
        rhinovirus_rows = c.fetchall()
        #res = parse_rows(rhinovirus_rows,"rhino",outcome)
        _rhinoprior_1year = res[0];_rhinoprior = res[1]

	norovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in ( " + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%norovirus%' order by agedays asc"
        c.execute(norovirus_q)
        norovirus_rows = c.fetchall()
	#res = parse_rows(norovirus_rows,"noro",outcome)
        _noroprior_1year = res[0];_noroprior = res[1]
	crassphage_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%crassphage%' order by agedays asc"
        c.execute(crassphage_q)
        crass_rows = c.fetchall()
	#res = parse_rows(crass_rows,"crass",outcome)
        _crassprior_1year = res[0];_crassprior = res[1]
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%parvo%' order by agedays asc"
        c.execute(_q)
        parvo_rows = c.fetchall()
        #res = parse_rows(parvo_rows,"parvo",outcome)
        _parvoprior_1year = res[0];_parvoprior = res[1]
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rotavirus%' order by agedays asc"
        c.execute(_q)
        rota_rows = c.fetchall()
        #res = parse_rows(rota_rows,"rota",outcome)
        _rotaprior_1year = res[0];_rotaprior = res[1]

	_q = "select avg(shannon) from vipie_findings_sample_mask_summary where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and agedays < 365"
        c.execute(_q)
        _row = c.fetchone()
	year1shannon = _row[0]
	if (year1shannon == None):
                year1shannon = 0
	_q = "select avg(shannon) from vipie_findings_sample_mask_summary where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "'"
        c.execute(_q)
        _row = c.fetchone()
        allshannon = _row[0]
	if (allshannon == None):
		allshannon = 0
	line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], str(_evprior), str(_evprior_1year), str(_coxbprior), str(_coxbprior_1year), str(_coxb3prior), str(_coxb3prior_1year), str(_coxb3nprior), str(_coxb3nprior_1year),str(_cmvprior), str(_cmvprior_1year), str(_noroprior), str(_noroprior_1year), str(_crassprior), str(_crassprior_1year), str(_parvoprior), str(_parvoprior_1year),str(_rotaprior), str(_rotaprior_1year), str(_rhinoprior), str(_rhinoprior_1year), ev_cv, str(year1shannon), str(allshannon), ia_pattern, country])
	outf.write(line + "\n")

"""	
print "xev"
print("\n".join(multiple_ev))
print "xcmv"
print("\n".join(multiple_cmv))
print "xnoro"
print("\n".join(multiple_noro))
print "xcrass"
print("\n".join(multiple_crass))
print "xparvo"
print("\n".join(multiple_parvo))
print "xrota"
print("\n".join(multiple_rota))
"""

conn.close()
print "case", _virus_case
print "control", _virus_control
print len(_virus_case_overlapped), _virus_case_overlapped
print len(_virus_control_overlapped), _virus_control_overlapped
 
