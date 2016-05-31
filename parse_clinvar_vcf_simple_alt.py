#!/usr/bin/python
import re
import os
import pdb
import copy

#orig_f = '20151110_clinvar_original.vcf';
#clinvar_f  = open('multiple_variant_same_loc', 'r')
#clinvar_f  = open('test', 'r')
#clinvar_f = open('20151110_clinvar_original.vcf.test', 'r')
clinvar_f = open('clinvar_20160203_papu_noY.vcf', 'r')
clv_o = open('clinvar_20160203_papu_noY_GRCh38.vcf', 'w')
clin = {}
#cln_allele = {}

conv ={}
convt = open('GCF_000001405.28.assembly.txt', 'r')
for cl in convt:
	if not cl.startswith('#'):
		aline = cl.strip().split('\t')
		conv[aline[6]] = aline[-1]
convt.close()

for aline in clinvar_f:
	if(not re.match('^#', aline)):
		aline = re.sub('\s+$', '', aline)
		#aline = repr(aline)
		#aline = re.sub(r'\\x2c_', r'\\\\x2c_', aline)
		#aline = re.sub(r'\\x3d_', r'\\\\x3d_', aline)
		content = aline.split('\t')
		aVar = {}
		clv_o.write(conv[content[0]])
		clv_o.write('\t')
		clv_o.write('\t'.join(content[1:]))
		clv_o.write('\n')
		alt_loc = conv[content[0]].split('_')
		var_chr = alt_loc[0]
		full_chr = conv[content[0]]
		#var_chr = re.sub('chr', '', var_chr)
		var_ncbi = ''
		if len(alt_loc) >1:
			var_ncbi = alt_loc[1]
		var_ct = ''
		if len(alt_loc) >2:
			var_ct = alt_loc[2]
		var_p   = content[1]
		var_rsid  = content[2]
		var_ref = content[3]
		var_alt = content[4]
		if( var_chr == 'chrX' or var_chr == 'X'):
			var_chr = 'chr23'
		elif (var_chr =='chrY' or var_chr == 'Y'):
			var_chr = 'chr24'
		elif(var_chr =='MT'):
			var_chr = 'chr25'
		
		#var_id = ':'.join([full_chr, var_p ])
		var_id = var_chr+'_'+var_ncbi+'_'+var_ct+':'+var_p
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
				if(pair[0] == 'CLNSIG'):
					#tmp = {}
					grp_sig = re.split('[,|]', pair[1])
					for s in grp_sig:
						tmp = {}
						if s == '0':
							s = 'Uncertain significance'
						elif s == '1':
							s = 'not provided'
						elif s == '2':
							s = 'Benign'
						elif s == '3':
							s = 'Likely benign'
						elif s == '4':
							s = 'Likely pathogenic'
						elif s == '5':
							s = 'Pathogenic'
						elif s == '6':
							s = 'drug response'
						elif s == '7':
							s = 'histocompatibility'
						elif s == '255':
							s = 'other'

						tmp[pair[0]] = s
							
						if var_id not in clin.keys():
							clin[var_id] = []
							clin[var_id].append(tmp)
						else:
							clin[var_id].append(tmp)

varid_o = open('ori_vcf_var_id', 'w')
for v in clin.keys():
	varid_o.write(v+'\n')
varid_o.close()
clv_o.close()

#json_f = open('multiple_json', 'r')
dig_log = open('comparison_log', 'w')
dig_log.write('json_line(from mongo)\torig_cont(from vcf)\n')

json_f = open('hg38_clinvar.alt_rnd_un.json', 'r')
#header = json_f.readline()

for al in json_f:
	al = re.sub('\s+$', '', al)
	
	al = re.sub('"', '', al)
	pos = re.search('c:(\d+),ct:([a-zA-Z]+),ncbi:([0-9a-zA-Z]+),p:(\d+)', al).groups()
	#if pos[0] ==24 and pos[1] == 624389:
	#	print('here')
	vid = 'chr'+pos[0]+'_'+ pos[2]+'_'+ pos[1]+':'+pos[3]
	json_val = ''
	al = re.sub('^{', '', al)
	al = re.sub('}$', '', al)
	fun = re.findall('{(.*?)}', al)## fetech one fun block
	outstr = ''
	#orig_cont= ':'.join([pos[0],pos[1], pos[2],pos[3]])
	orig_cont = vid
	json_val = ''
	alt_alle_idx = 0
	for anno in clin[vid]:
		asJson_cont = ''
		orig_cont = orig_cont+ ' ' + anno['CLNSIG']
	
	for fidx in range(1,len(fun)):
		field = fun[fidx].split(',')
		for fd in field:# for each sub-block in f
			kv = fd.split(':')
			if kv[0] == 'CLNSIG':
				f_cont = kv[1]
				json_val = json_val + ' ' + f_cont
				break
	json_line = ''.join([ vid ,json_val ])
	dig_log.write(json_line+'\t'+orig_cont)	
	if(not json_line == orig_cont):
		dig_log.write('\t'+pos[0]+':'+pos[1] +'    -> wrong')
	dig_log.write('\n')

dig_log.close()
