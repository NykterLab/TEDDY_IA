inf = open("oneyear_plasma_avg_shannon.csv","r")
lines = inf.readlines()
outcome0 = {}
outcome1 = {}

for l in lines:
	tk = l.strip().split(",")
	outcome = tk[2]
	ind = tk[1]
	if (outcome == "1"):
		outcome1[ind] = tk
	else:
		outcome0[ind] = tk
inf.close()

in1only = {}
shannon0_larger = 0
for _ind in outcome1:
	ind1 = outcome1[_ind]
	if (outcome0.get(_ind) != None):
		ind0 = outcome0[_ind]
		if (float(ind0[-1]) > float(ind1[-1])):
			shannon0_larger +=1
	else:
		in1only[_ind] = ind1	

print shannon0_larger, len(outcome1), len(in1only)

in0only = {}
for _ind in outcome0:
        ind1 = outcome0[_ind]
        if (outcome1.get(_ind) == None):
                in0only[_ind] = ind1

print "1 year shannon avg: control>case, case_incident_count, control_incident_count, shannon_diversity in case incidents only, in control incidents only,"
print shannon0_larger, len(outcome1), len(outcome0), len(in1only), len(in0only)
		
