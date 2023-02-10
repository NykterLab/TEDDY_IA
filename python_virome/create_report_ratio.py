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
design = sys.argv[1]
tissue = sys.argv[2]
hit_cutoff = 2
if (tissue == "stool"):
	hit_cutoff = 5

typespecific = ""
otypes = "'wgs', 'cv'"
if (tissue == "stool"):
	otypes = "'cv', 'pv'"
if (len(sys.argv) >= 4):
	otypes = "'" + sys.argv[3] + "'"
	typespecific = sys.argv[3] + "_"

label = ""
if (len(sys.argv) == 5):
        label = sys.argv[4]

tsne = open(typespecific + design + "_" + tissue + "_numeric_viral_exposureds_ratio"+label+".csv", "w")
outf = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_ratio"+label+".csv", "w")
adj_paired = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_adj_paired_ratio"+label+".csv", "w")
adj_sample = open(typespecific + design + "_" + tissue + "_viral_exposured_agedays_adj_sample_ratio"+label+".csv", "w")
itable = "mp183_case_cntrl_t1d_list1"
if (design.find("ia") != -1):
	itable = "mp183_case_cntrl_ia_list1"
_ind_map = {}

"""
#stool cv
case_sample_map = open("/Users/jakelin/Desktop/TEDDY/%s/%s_%s/case/sample_mask_microbiome_map.tsv" %(design, tissue, otype), "r")
control_sample_map = open("/Users/jakelin/Desktop/TEDDY/%s/%s_%s/control/sample_mask_microbiome_map.tsv" %(design, tissue, otype), "r")
#stool pv
control_samples = {}
for l in control_sample_map.readlines():
	tk = l.strip().split("\t")
	if (tk[0].find("sample") == -1):
		ind = tk[7]
		if (control_samples.get(ind) == None):
			control_samples[ind] = "'" + "'"

#for every case_ind in case and control, there will be a set of sample_mask_ids correspond to the list input files
"""

expairs = 0
c.execute("select mask_id, case_ind, outcome from " + itable + " order by case_ind, outcome" )
srows = c.fetchall()
header = "mask_id,case_ind,t1d_outcome,ia_outcome,hla_category,gender,fdr,ia_first_diag_agedys,t1d_diag_agedys,delta_ia_t1d,persist_conf_gad_agedys,persist_conf_ia2a_agedys,persist_conf_miaa_agedys,enterovirus_priorx,enterovirus_year1,coxsackiea_priorx,coxsackiea_year1,coxsackieb_priorx,coxsackieb_year1,coxsackieb3_priorx,coxsackieb3_year1,coxsackie_notb3_priorx,coxsackie_notb3_year1,cytomegalovirus_priorx,cytomegalovirus_year1,norovirus_priorx,norovirus_year1,crassphage_priorx,crassphage_year1,parvovirus_priorx,parvovirus_year1,rotavirus_priorx,rotavirus_year1,rhinovirus_priorx,rhinovirus_year1,shannon_year1,shannon_overall,ia_pattern,country,celiac_outcome\n"
outf.write(header)
adj_paired.write(header)
adj_sample.write(header)

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


def parse_rows(_rows, _virus, outcome):
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
        #otypes = "'wgs', 'cv'"
	#if (tissue == "stool"):
	#	otypes = "'wgs', 'cv', 'pv'"	
	enterovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%enterovirus%' and lower(genus_desc) not like '%rhinovirus%') order by agedays asc"
	c.execute(enterovirus_q)
	enterovirus_rows = c.fetchall()
	_found = len(enterovirus_rows)
	res = parse_rows(enterovirus_rows,"ev",outcome)	
	_evprior_1year = res[0]
	_evprior = res[1]
	
        #_fcoxb = 0	
        coxb_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') order by agedays asc"
        c.execute(coxb_q)
        coxb_rows = c.fetchall()
        res = parse_rows(coxb_rows,"coxb",outcome)
        _coxbprior_1year = res[0]
        _coxbprior = res[1]

	coxa_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie a%' or lower(genus_desc) like '%coxsackievirus a%') order by agedays asc"
        c.execute(coxa_q)
        coxa_rows = c.fetchall()
        res = parse_rows(coxa_rows,"coxa",outcome)
        _coxaprior_1year = res[0]
        _coxaprior = res[1]

        #_fcoxb3 = 0
        coxb3_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b3%' or lower(genus_desc) like '%coxsackievirus b3%') order by agedays asc"
        c.execute(coxb3_q)
        coxb3_rows = c.fetchall()
        res = parse_rows(coxb3_rows,"coxb3",outcome)
        _coxb3prior_1year = res[0]
        _coxb3prior = res[1]

        coxb3n_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and (lower(genus_desc) like '%coxsackie b%' or lower(genus_desc) like '%coxsackievirus b%') and lower(genus_desc) not like '%coxsackie b3%' and lower(genus_desc) not like '%coxsackievirus b3%' order by agedays asc"
        c.execute(coxb3n_q)
        coxb3n_rows = c.fetchall()
        res = parse_rows(coxb3n_rows,"coxb3n",outcome)
        _coxb3nprior_1year = res[0]
        _coxb3nprior = res[1]

        _fcmv = 0
	cytomegalovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%cytomegalo%' order by agedays asc"
        c.execute(cytomegalovirus_q)
        cmvvirus_rows = c.fetchall()
        res = parse_rows(cmvvirus_rows,"cmv",outcome)
        _cmvprior_1year = res[0];_cmvprior = res[1]
	#;_cmv_case = res[2];_cmv_control = res[3]

	#rhinovirus - most common, common cold
	rhinovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rhinovirus%' order by agedays asc"
        c.execute(rhinovirus_q)
        rhinovirus_rows = c.fetchall()
        res = parse_rows(rhinovirus_rows,"rhino",outcome)
        _rhinoprior_1year = res[0];_rhinoprior = res[1]

	"""
	_cmvfound = len(cmvvirus_rows)
        _cmvprior = 0
        _cmvprior_1year = 0
        if (_cmvfound > 1):
		_fcmv = 1
                _age0_exposed = int(cmvvirus_rows[0][0])
                _cmvprior = 1
		#_cmvprior_1year = 1#_age0_exposed
                if (_age0_exposed < 365):
                        _cmvprior_1year = 1#_age0_exposed
                #here is the piece that can calculate acute versus mild infection
                _se = {}
                for _e in cmvvirus_rows:
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
                if (len(_se) > 1):
                        _cmvprior = len(_se)
                        multiple_cmv.append(outcome + " " + mask_id + " multiple_cmv_exposures " + str(_cmvprior))
		if (outcome == "1"):
			_cmv_case += 1
		else:
			_cmv_control += 1
		if (_fev == 1):
			ev_cv = "1"
	"""
	norovirus_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in ( " + otypes + ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%norovirus%' order by agedays asc"
        c.execute(norovirus_q)
        norovirus_rows = c.fetchall()
	res = parse_rows(norovirus_rows,"noro",outcome)
        _noroprior_1year = res[0];_noroprior = res[1]
	#;_noro_case = res[2];_noro_control = res[3]
        """_norofound = len(norovirus_rows)
        _noroprior = 0
        _noroprior_1year = 0
        if (_norofound > 1):
                _age0_exposed = int(norovirus_rows[0][0])
                _noroprior = 1
		#_noroprior_1year = _age0_exposed
                if (_age0_exposed < 365):
                        _noroprior_1year = 1#_age0_exposed
                #here is the piece that can calculate acute versus mild infection
                _se = {}
                for _e in norovirus_rows:
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
                if (len(_se) > 1):
                        _noroprior = len(_se)
                        multiple_noro.append(outcome + " " + mask_id + " multiple_noro_exposures " + str(_noroprior))
                if (outcome == "1"):
                        _noro_case += 1
                else:
                        _noro_control += 1
		if (_fev == 1 or _fcmv == 1):
                        ev_cv = "1"
	"""
	crassphage_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%crassphage%' order by agedays asc"
        c.execute(crassphage_q)
        crass_rows = c.fetchall()
	res = parse_rows(crass_rows,"crass",outcome)
        _crassprior_1year = res[0];_crassprior = res[1]
	"""
	crass_rows = c.fetchall()
        _crassfound = len(crass_rows)
        _crassprior = 0
        _crassprior_1year = 0
        if (_crassfound > 1):
                _age0_exposed = int(crass_rows[0][0])
                _crassprior = 1
                if (_age0_exposed < 365):
                        _crassprior_1year = 1#_age0_exposed
                #here is the piece that can calculate acute versus mild infection
                _se = {}
                for _e in crass_rows:
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
                if (len(_se) > 1):
                        _crassprior = len(_se)
                        multiple_crass.append(outcome + " " + mask_id + " multiple_crass_exposures " + str(_crassprior))
                if (outcome == "1"):
                        _crass_case += 1
                else:
                        _crass_control += 1
                if (_fev == 1 or _fcmv == 1):
                        ev_cv = "1"
	"""
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%parvo%' order by agedays asc"
        c.execute(_q)
        parvo_rows = c.fetchall()
        res = parse_rows(parvo_rows,"parvo",outcome)
        _parvoprior_1year = res[0];_parvoprior = res[1]
	"""
	_parvofound = len(parvo_rows)
        _parvoprior = 0
        _parvoprior_1year = 0
        if (_parvofound > 1):
                _age0_exposed = int(parvo_rows[0][0])
                _parvoprior = 1
                if (_age0_exposed < 365):
                        _parvoprior_1year = 1#_age0_exposed
                _se = {}
                for _e in parvo_rows:
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
                if (len(_se) > 1):
                        _parvoprior = len(_se)
                        multiple_parvo.append(outcome + " " + mask_id + " multiple_parvo_exposures " + str(_parvoprior))
                if (outcome == "1"):
                        _parvo_case += 1
                else:
                        _parvo_control += 1
                if (_fev == 1 or _fcmv == 1):
                        ev_cv = "1"
	"""
	_q = "select agedays, hits, genus_desc, sample_mask_id  from vipie_findings where lower(tissue) = '" + tissue + "' and lower(otype) in (" + otypes +  ") and mask_id = '" + mask_id + "' and lower(genus_desc) like '%rotavirus%' order by agedays asc"
        c.execute(_q)
        rota_rows = c.fetchall()
        res = parse_rows(rota_rows,"rota",outcome)
        _rotaprior_1year = res[0];_rotaprior = res[1]
	"""
	_rotafound = len(rota_rows)
        _rotaprior = 0
        _rotaprior_1year = 0
        if (_rotafound > 1):
                _age0_exposed = int(rota_rows[0][0])
                _rotaprior = 1
                if (_age0_exposed < 365):
                        _rotaprior_1year = 1#_age0_exposed
                _se = {}
                for _e in rota_rows:
                        smid = _e[-1]
                        if (_se.get(smid) == None):
                                _se[_e[-1]] = float(_e[1])
                        else:
                                _se[_e[-1]] = _se[_e[-1]] + float(_e[1])
                if (len(_se) > 1):
                        _rotaprior = len(_se)
                        multiple_rota.append(outcome + " " + mask_id + " multiple_rota_exposures " + str(_rotaprior))
                if (outcome == "1"):
                        _rota_case += 1
                else:
                        _rota_control += 1
                if (_fev == 1 or _fcmv == 1):
                        ev_cv = "1"
	"""

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
	if (case_control_ratio.get(case_ind) ==  None or (case_control_ratio.get(case_ind) < .7) or (case_control_ratio.get(case_ind) > 1.3)):#) and (control_f[case_ind] < 5)):
		qexcluded[case_ind] = 1
		expairs += 1
		continue
	case2control = case_control_ratio[case_ind]
	#if (outcome == "0" and case2control != 1):
	#	case2control = 1 - case2control
	virus_corrected = [str(x/case2control) for x in [_evprior, _coxaprior, _coxbprior, _coxb3prior, _coxb3nprior, _cmvprior, _noroprior, _crassprior, _parvoprior, _rotaprior, _rhinoprior]]	
	
	line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], str(_evprior), str(_evprior_1year), str(_coxaprior), str(_coxaprior_1year), str(_coxbprior), str(_coxbprior_1year), str(_coxb3prior), str(_coxb3prior_1year), str(_coxb3nprior), str(_coxb3nprior_1year),str(_cmvprior), str(_cmvprior_1year), str(_noroprior), str(_noroprior_1year), str(_crassprior), str(_crassprior_1year), str(_parvoprior), str(_parvoprior_1year),str(_rotaprior), str(_rotaprior_1year), str(_rhinoprior), str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome])

		
	outf.write(line + "\n")

        line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], virus_corrected[0], str(_evprior_1year), virus_corrected[1], str(_coxaprior_1year), virus_corrected[2], str(_coxbprior_1year), virus_corrected[3], str(_coxb3prior_1year), virus_corrected[4], str(_coxb3nprior_1year), virus_corrected[5], str(_cmvprior_1year), virus_corrected[6], str(_noroprior_1year), virus_corrected[7], str(_crassprior_1year), virus_corrected[8], str(_parvoprior_1year), virus_corrected[9], str(_rotaprior_1year), virus_corrected[10], str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome])
	adj_paired.write(line + "\n")
	numfiles = 0
	if (outcome == "1"):
		numfiles = case_f[case_ind]
	if (outcome == "0"):
                numfiles = control_f[case_ind]
        virus_corrected = [str(x/(numfiles*1.0)) for x in [_evprior, _coxaprior, _coxbprior, _coxb3prior, _coxb3nprior, _cmvprior, _noroprior, _crassprior, _parvoprior, _rotaprior, _rhinoprior]]
	line = ",".join([mask_id, case_ind, t1d_outcome, outcome, hla, gender, fdr, key_days[0], key_days[1], key_days[5], key_days[2], key_days[3], key_days[4], virus_corrected[0], str(_evprior_1year), virus_corrected[1], str(_coxaprior_1year), virus_corrected[2], str(_coxbprior_1year), virus_corrected[3], str(_coxb3prior_1year), virus_corrected[4], str(_coxb3nprior_1year), virus_corrected[5], str(_cmvprior_1year), virus_corrected[6], str(_noroprior_1year), virus_corrected[7], str(_crassprior_1year), virus_corrected[8], str(_parvoprior_1year), virus_corrected[9], str(_rotaprior_1year), virus_corrected[10], str(_rhinoprior_1year), str(year1shannon), str(allshannon), ia_pattern, country, celiac_outcome])
        adj_sample.write(line + "\n")

	tsne.write(" ".join([str(_evprior), str(_coxaprior), str(_coxbprior), str(_coxb3prior), str(_coxb3nprior), str(_cmvprior), str(_noroprior), str(_crassprior), str(_parvoprior),str(_rotaprior), str(_rhinoprior), str(year1shannon), str(allshannon)]) + "\n")
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
outf.close()
tsne.close()
adj_paired.close() 
adj_sample.close() 
print "excluded pairs", expairs/2
print qexcluded
