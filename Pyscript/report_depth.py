def report_depth(table,title,seq,step,ori):
	fh_depth = open('test.depth', 'w')
	standard_length = 22
	depth_sum = {}
	reports = []
	four_bases = ('A','T', 'C', 'G')
	reports = table[seq]
	if ori == 'f':
		for_rest_len = standard_length - step
		for x in range(0,for_rest_len):
			real_position = x + 1
			fh_depth.write(title + "\t" + str(real_position) + "\t") 
			for y in range(0,len(reports)):
				if reports[y][x] in depth_sum.keys():
					depth_sum[reports[y][x]] += 1
				else:
					depth_sum[reports[y][x]] = 0

			for base in four_bases:
				if base in depth_sum.keys():
					fh_depth.write(base + ":" + str(depth_sum[base]) + "\t")
				else:
					fh_depth.write(base + ":0" + "\t")
			fh_depth.write("\n")
			depth_sum = {}
		
	else:
		rev_rest_len = standard_length -1
		for x in range(step,rev_rest_len):
			real_position = standard_length - step + x + 1
			fh_depth.write(title + "\t" + real_position + "\t")
			for y in range(0,len(reports)):
				if reports[y][x] in depth_sum.keys():
					depth_sum[reports[y][x]] += 1
				else:
					depth_sum[reports[y][x]] = 0

			for base in four_bases:
				if base in depth_sum.keys():
					fh_depth.write(base + ":" + str(depth_sum[base]) + "\t")
				else:
					fh_depth.write(base + ":0" + "\t")
			fh_depth.write("\n")
			depth_sum = {}


consensus_depth = {'TAAAAATAATTATAAAATAACT': 
[['A','T','C','G','C','G','G','G','G','A','A','C','A','C','C','G','T','A','G','C','T','A'],
['A','T','C','G','C','G','G','C','G','A','A','C','A','C','G','G','T','A','G','C','A','A'],
['A','T','C','C','C','G','G','G','G','A','A','C','A','C','C','G','T','A','G','C','T','A'],
['A','T','C','G','C','G','G','G','G','A','A','C','A','C','C','G','T','A','G','C','T','N'],
['A','T','C','G','C','G','A','G','G','A','A','C','A','C','C','G','T','A','G','C','A','A']]}

report_depth(consensus_depth,"text","TAAAAATAATTATAAAATAACT",5,'f')
