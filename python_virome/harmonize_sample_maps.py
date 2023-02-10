#stool cv
import sys

tissue = sys.argv[1]
otype = sys.argv[2]
case_sample_map = open("/Users/jakelin/Desktop/TEDDY/ia_list1/%s_%s/case/sample_mask_microbiome_map.tsv" %(tissue, otype), "r")
control_sample_map = open("/Users/jakelin/Desktop/TEDDY/ia_list1/%s_%s/control/sample_mask_microbiome_map.tsv"  %(tissue, otype), "r")
case_sample_harmony = open("/Users/jakelin/Desktop/TEDDY/ia_list1/%s_%s/case/sample_mask_microbiome_map_harmony.tsv"  %(tissue, otype), "w")
control_sample_harmony = open("/Users/jakelin/Desktop/TEDDY/ia_list1/%s_%s/control/sample_mask_microbiome_map_harmony.tsv"  %(tissue, otype), "w")


#stool pv
days = 90
if (tissue == 'plasma'):
	days = 180
control_samples = {}
control_lines = {}
for l in control_sample_map.readlines():
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                if (control_samples.get(ind) == None):
                        control_samples[ind] = tk[0] + "_" + tk[4]
		else:
			control_samples[ind] = control_samples[ind] + ";" + tk[0] + "_" + tk[4]
		control_lines[tk[0] + "_" + tk[4]] = l

control_sample_map.close()

case_samples = {} 
case_lines = case_sample_map.readlines()
for l in case_lines:
        tk = l.strip().split("\t")
        if (tk[0].find("sample") == -1):
                ind = tk[7]
                if (case_samples.get(ind) == None):
                        case_samples[ind] = l
                else:
                        case_samples[ind] = case_samples[ind] + ";" + l

case_used = {}
for ind in control_samples:
	control_files = control_samples[ind].split(";")
	if (case_samples.get(ind) == None):
		print "no case_samples, skip", ind
		continue
	case_lines = case_samples[ind].split(";")
	cc = 0
	for s in control_files:		
		control_age = int(s.split("_")[1])
		if (ind == "99"):
			print ind, s, case_used
		case_used = {}	
		if (cc < len(case_lines)):
			case_line = case_lines[cc]
			case_age = case_line.split("\t")[4]
			if (abs(int(case_age) - control_age) <= days):
				case_sample_harmony.write(case_line)
				control_sample_harmony.write(control_lines[s])
				cc += 1
			else:
				previous_age = -1
				for ccc in range(cc+1, len(case_lines)):
					case_line = case_lines[ccc]					
                        		case_age = case_line.split("\t")[4]
                        		if (abs(int(case_age) - control_age) <= days):
                                		case_sample_harmony.write(case_line)
						control_sample_harmony.write(control_lines[s])
						cc = ccc+1
						break
case_sample_harmony.close()
 
control_sample_harmony.close()

