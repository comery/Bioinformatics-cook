#!/usr/bin/env python

import os,sys,time
from optparse import OptionParser
import subprocess

#PATH = os.path.dirname(sys.argv[0])

#bwa = PATH + '/bwa'

def work(cmd):

	os.system(cmd)
	sys.stderr.write('[%s] Finish CMD: %s \n' % (time.ctime( ),cmd))
	
#def bowtie(PEargs):
	

# deal with bwa aln and sampe
# should give all the usable arguments or can read form configure files


def bwaAlignPE(PEargs):
	
	if len(PEargs) < 1:
		print "usage: PE -h, --help"

		sys.exit(0)

	usage = "usage: PE [options] ref.fa xx_1.fq xx_2.fq"
	parser = OptionParser(usage)
	parser.add_option("-d", "--outdir", dest="outputdir",
			help="out put directory")
	parser.add_option("-s", "--outfileprefix", dest="outputfilekeyname",
			help="out put file prefix")
	parser.add_option("-v", "--MaxDiff",dest="MaxMismatch",default=0.04,
			help="Max mismatch number [0.04]")
	parser.add_option("-x", "--MaxInsertSize",dest="MaxInSize",default=800,
			help="Max insert size [800]")
	parser.add_option("-l", type="int",dest="seed",default=35,
			help="seed size [35]")
	parser.add_option("-R", type="int", dest="mostHits",default=30,
			help="stop searching when there are >INT equally best hits [30]")
        parser.add_option("-t", type="int",dest="threads",default=4, 
			help="threads nubmer [4: total 4X2]")
	parser.add_option("-e", type="int", dest="MaxGapSize", default=-1,
			help="maximum number of gap extensions, -1 for disabling long gaps")
	parser.add_option("-i", type="int", dest="EndwithoutGap", default=5,
			help="do not put an indel within INT bp towards the ends")
	parser.add_option("-p", type="str", dest="BWA", default="/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/bwa-0.7.0/bwa",
			help="bwa. Default /ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/bwa-0.7.0/bwa")

	(options, args) = parser.parse_args(args=PEargs)
	if len(args) < 1:
		print "usage: PE -h, --help"
		sys.exit(0)

	if options.outputdir[-1] != '/':
		options.outputdir += '/'
	
	sai1= args[1].split('/')[-1]
	prefix=sai1
	sai1 = options.outputdir + sai1 + '.sai'
	sai2=args[2].split('/')[-1]
	sai2 = options.outputdir + sai2 + '.sai'
	prefix = options.outputdir +  options.outputfilekeyname


	cmd1 ='%s aln -m 200000 -o 1 -e %s -i %s -L -I -t %s -n %s -R %s %s %s -f %s' % (options.BWA, options.MaxGapSize, options.EndwithoutGap, options.threads, options.MaxMismatch, options.mostHits, args[0], args[1], sai1)
	cmd2 ='%s aln -m 200000 -o 1 -e %s -i %s -L -I -t %s -n %s -R %s %s %s -f %s' % (options.BWA, options.MaxGapSize, options.EndwithoutGap, options.threads, options.MaxMismatch, options.mostHits, args[0], args[2], sai2)

	
	cmd3='%s sampe -a %s %s %s %s %s %s | gzip >%s.sam.gz && rm %s %s' % (options.BWA, options.MaxInSize, args[0], sai1, sai2, args[1], args[2], prefix, sai1, sai2)
	
	# here we fork tow child of process to do 'align' 
	child = os.fork()

	if child == 0:
		work(cmd1)
		os._exit(0)
	work(cmd2)
	print os.getpid()


	while True:
		(id, stat) = os.waitpid(child, os.EX_OK)
		if stat == 0:
			break
		
	work(cmd3)

####################################################################################	
USAGE= sys.argv[0] + """ [command] <option>\n
   PE	Pairend alignment 
   SE	Singleend aligment
   MS	Merge by chromosome and sort the SAM result, support read sam.gz

"""	
####################################################################################
if __name__ == "__main__":

	
	if len(sys.argv) >1 and sys.argv[1] == "PE":
		
		bwaAlignPE(sys.argv[2:]) 
	else:
		print USAGE		
		sys.exit(0)
