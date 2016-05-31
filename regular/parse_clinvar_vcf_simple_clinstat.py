#!/usr/bin/python
import re
import os
import pdb
import copy

#orig_f = '20151110_clinvar_original.vcf';
#clinvar_f  = open('multiple_variant_same_loc', 'r')
#clinvar_f  = open('test', 'r')
#clinvar_f = open('20151110_clinvar_original.vcf.test', 'r')
clinvar_f = open('clinvar_20160203.vcf', 'r')

clin = {}
#cln_allele = {}

for aline in clinvar_f:
	if(not re.match('^#', aline)):
		aline = re.sub('\s+$', '', aline)
		#aline = repr(aline)
		#aline = re.sub(r'\\x2c_', r'\\\\x2c_', aline)
		#aline = re.sub(r'\\x3d_', r'\\\\x3d_', aline)
		content = aline.split('\t')
		aVar = {}
		var_chr = content[0]
		var_p   = content[1]
		var_rsid  = content[2]
		var_ref = content[3]
		var_alt = content[4]
		if( var_chr== 'X'):
			var_chr = '23'
		elif (var_chr =='Y'):
			var_chr = '24'
		elif(var_chr =='MT'):
			var_chr = '25'
		
		var_id = ':'.join([var_chr, var_p ])
		
		alt_allele = re.split(',', var_alt)
		num_allele = len(re.split(',', var_alt))
		grp_size = []
		total_size = 0
		cln_allele = {}
		info = content[-1].split(';')
		for ele in info:  ## parse the info column of VCF file
			pair = ele.split('=')
			if(len(pair) ==2):
				## deal with multiple annotation for one allele
				##specific for clinvar
				if(pair[0] == 'CLNREVSTAT'):
					#tmp = {}
					grp_sig = re.split('[,|]', pair[1])
					for s in grp_sig:
						tmp = {}
						tmp[pair[0]] = s
							
						if var_id not in clin.keys():
							clin[var_id] = []
							clin[var_id].append(tmp)
						else:
							clin[var_id].append(tmp)

#json_f = open('multiple_json', 'r')
dig_log = open('comparison_log', 'w')
dig_log.write('json_line(from mongo)\torig_cont(from vcf)\n')

json_f = open('hg38_clinvar.regular.json', 'r')
header = json_f.readline()

for al in json_f:
	al = al.strip()
	al = re.sub('\s', '', al)
	
	al = re.sub('"', '', al)
	pos = re.search('p:(\d+),c:(\d+)', al).groups()
	#if pos[0] ==24 and pos[1] == 624389:
	#	print('here')
	vid = ':'.join([pos[1], pos[0]])
	json_val = ''
	al = re.sub('^{', '', al)
	al = re.sub('}$', '', al)
	fun = re.findall('{(.*?)}', al)## fetech one fun block
	outstr = ''
	orig_cont= ':'.join([pos[1],pos[0]])
	json_val = ''
	alt_alle_idx = 0
	for anno in clin[':'.join([pos[1],pos[0]])]:
		asJson_cont = ''
		orig_cont = orig_cont+ ' ' + anno['CLNREVSTAT']
	
	for fidx in range(1,len(fun)):
		field = fun[fidx].split(',')
		for fd in field:# for each sub-block in f
			kv = fd.split(':')
			if kv[0] == 'CLNREVSTAT':
				f_cont = kv[1]
				json_val = json_val + ' ' + f_cont
				break
	json_line = ''.join([ pos[1],":",pos[0] ,json_val ])
	dig_log.write(json_line+'\t'+orig_cont)	
	if(not json_line == orig_cont):
		dig_log.write('\t'+pos[1]+':'+pos[0] +' -> wrong')
	dig_log.write('\n')

dig_log.close()
