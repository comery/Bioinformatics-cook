#!/usr/bin/python 
import os , sys , time , re , gzip
from optparse import OptionParser
import commands

def RunShell(Run_script):
        proRunID = {}
        working = 1
        while  working :
                if len (proRunID) >= 4 or len (Run_script) == 0 :
                        (pid , stat) = os.wait()
                        if stat % 255 :
                                print >>stderr , "run error"
                                exit(1)
                        proRunID.pop( pid )

                i = len(proRunID)
                while  i < int (parse.process) and len (Run_script) : ## add int
                        script = ""
                        script = Run_script.pop()
#                       print script 
                        child_pid = os.fork()
                        if child_pid == 0 :
                                ( stat , cmd ) = commands.getstatusoutput ( script ) #############
# os.system(script) 
                                sys.exit(0)
                        proRunID [child_pid ] = 1
                        i+=1
                if len(proRunID) == 0 and len(Run_script) == 0  :
                        working = 0
 
def GetLib(name):
        l = re.search("[A-Za-z0-9]+" , name)
        if l :
                return name[l.start()+4:l.end()]
        print >>sys.stderr , "bam file name error:"+ name +" \nplease keep fastq infomation !\n"
        sys.exit(1)


if __name__ == '__main__':
	usage = "Usage: sortmerge [options] 1.sam.gz 2.sam.gz .."
        opt = OptionParser(usage)
        opt.add_option("-o" , "--outdir" , dest = "outdir" , 
                        help = "[str] all chr output this dir default './'" , default = "./")
        opt.add_option("-p" , "--process" , dest = "process" , 
                        help = "[int] process to use for sort , rmdup " , default = 4)
        opt.add_option("-s" , "--samplename" , dest = "sname" , 
                        help = "[str] sample name ,default BWA" , default = "BWA")
        opt.add_option("-t" , "--samtools" , dest = "stools" , 
                        help = "[str] samtools soft ,default /ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/samtools-0.1.18/samtools" , 
                        default = "/ifs2/BC_GAG/Group/zhuangzhenhua/s_ware/samtools-0.1.18/samtools")

        (parse, args) = opt.parse_args()

	if(len(args)<1):
		print "please use -h or --help"
		sys.exit(0)

        print >>sys.stderr ,"sortmerge sam format :"
        print >>sys.stderr ,"\t" + "outdir     : "+parse.outdir
        print >>sys.stderr ,"\t" + "samplename : "+parse.sname
        print >>sys.stderr ,"\t" + "process    : "+str(parse.process)
        print >>sys.stderr ,"\t" + "samtools   : "+parse.stools
        samtools = parse.stools

        parse.outdir =  parse.outdir + "/" + parse.sname 
        rootOutDir = parse.outdir

        (stat ,cmd) =  commands.getstatusoutput ("mkdir -p  " + parse.outdir +"&& rm -f "+parse.outdir+"/"+"*unmap.sam*" )
        if stat%255 :
                print >>stderr , "error:"+cmd
                exit(1)
########get sam head chr infor ###################################################
        samHead = parse.outdir+"/head.sam"
        (stat ,cmd) =  commands.getstatusoutput ("rm -f "+samHead+"&&"+samtools + " view -H -S " + args[0] + "|tee -a "+samHead)
        #(stat ,cmd) =  commands.getstatusoutput (samtools + " view -H -S " + args[0] +" > "+samHead+" ")
        if stat % 255 : 
                print >>sys.stderr , "samtools error:" + cmd
                sys.exit(1)
        Head = cmd.split( "\n" )
        Headinfo = cmd
        lib={}
        chrinfo = {} 
        for samname in args:
                l = GetLib(samname)
                if lib.has_key(l): lib[ l ] += 1 
                else : lib [ l ] = 1
        for line in Head:
               if line.split("\t")[0] == "@SQ":
                        chrname = line.split( "\t" )[1].split("SN:")[1]
                        chrinfo [chrname] = line
                        (stat , cmd ) = commands.getstatusoutput ("mkdir -p "+parse.outdir+"/"+chrname + "&& rm -f " +parse.outdir +"/"+chrname + "/chr*.sam ")
                        if  stat % 255 :
                                print >>sys.stderr , "mkdir "+parse.outdir+"/"+chrname+" error:"+cmd+"\n"
                                sys.exit(1)
        awkSplitBam = {} 
        for samname in args :   
                libname = GetLib( samname )
                awkSplitBam.setdefault(libname , "" )
                if len(awkSplitBam[libname]) :
                                awkSplitBam[libname] += ";"
###add rm *unmap.sam* 2011.9.15 error   mv to mkdir outdir  2012.2.15
                awkSplitBam[libname] += "zcat -f "+ samname + "| awk '{if ($0 !~/^@/){if ($3 != \"*\"){ print $0 >> \""+parse.outdir+"/\"$3\"/\"$3\"."+libname+".sam\"}else {print $0 >> \"" + parse.outdir + "/"+os.path.basename(samname)+".PE.unmap.sam\" } }}';gzip "  + parse.outdir + "/"+os.path.basename(samname)+".PE.unmap.sam"
        print >>sys.stderr ,"[" + time.ctime() + "]split sam start" 
        RunShell(awkSplitBam.values()) ##################
        print >>sys.stderr ,"[" + time.ctime() + "]split sam OK!!" 
##############################################################
        chrScript={}
        for chrname in chrinfo :
                chrScript[chrname] = ""
                chrAllfile = ""
                chrdirName = parse.outdir +"/"+chrname+"/"+chrname
                finalBam = parse.outdir +"/"+chrname+"/"+parse.sname+"."+chrname
                for libname in lib :
                        Chrlib = chrdirName+"."+libname
                        chrScript[ chrname ] += "cat "+ samHead + " "+ Chrlib + ".sam|"+ samtools +" view -Sb - |"+samtools +" sort - "+Chrlib +" && "+samtools+ " rmdup "+Chrlib+".bam "+Chrlib+".rmdup.bam 2>"+Chrlib + ".rmdup.log ;rm " + Chrlib + ".sam ;"
                        chrAllfile += Chrlib+ ".rmdup.bam "
                if len(lib) > 1 :
                        chrScript[ chrname ] +=  samtools + " merge "+finalBam +".bam "+chrAllfile+ ";"+samtools+" mapstat "+finalBam +".bam > "+finalBam+".bam.stat "#;"rm " + chrAllfile ;
                else:
                        chrScript[ chrname ] +=  "mv "+chrAllfile+ " " +finalBam +".bam ;"+samtools+" mapstat "+finalBam +".bam > "+finalBam+".bam.stat";
                chrScript [chrname ] += ";rm "+chrdirName+"*bam"
                chrScript[ chrname ] += "; " +samtools +" depth " + " " + finalBam + ".bam 2> "+finalBam +".Dstat|gzip > "+ finalBam +".depth.gz"
        print >>sys.stderr ,"[" + time.ctime() + "]sort merge sam start" 
        RunShell(chrScript.values()) #####
        print >>sys.stderr ,"[" + time.ctime() + "]sort merge sam end" 
#exit(0)

        mapstat =[ 
    #     "total reads",
    #     "total bases",
     #    "read1",
     #    "read2",
         "mapped reads",
         "mapped bases",
         "PE mapped reads[with itself and mate mapped]",
         "PE bases",
         "SE mapped reads[singletons]",
         "SE bases",
         "properly paired reads",
         "properly paired bases",
         "trimed reads",
         "trimed bases",
         "insertion envent reads",
         "insertion bases",
         "deletion envent reads",
         "deletion bases",
         "with mate mapped to a different chr",
         "with mate mapped to a different chr (mapQ>=5)",
         "QC failure",
         "duplicates"
         ]
        duplicate = {} 
        infoStat = {}
        for chrname in chrinfo:
                finalBamStat = parse.outdir +"/"+chrname+"/"+parse.sname+"."+chrname+ ".bam.stat"
                statF = open (finalBamStat)
                i = 0
                try:
                        while True :
                                line = statF.readline()
                                if line :
                                        infoStat.setdefault(mapstat[i] , 0.0 )
                                        infoStat[mapstat[i]] += int( line.split()[0] )
                                        i+=1
                                else :        
                                       break 
                finally:
                        statF.close()
                for libname in lib:
                        Chrlib = parse.outdir +"/"+chrname+"/"+chrname+"."+libname  
                        dup = open(Chrlib+".rmdup.log")
                        duplicate.setdefault(parse.sname ,[0.0 , 0.0 ] )
                        try:
                                line = dup.readline()
				check = 0 #
				while line :
                                	#chrlibdup = dup.readline()
					chrlibdup = line 
					line = dup.readline()
					check = 1 # add 2011.11.23 
                                if chrlibdup and check:
                                        duplicate.setdefault(libname , [0.0 , 0.0 ] )
                                        chrlibdup = chrlibdup.split()
                                        func = lambda x , y : float(x)+ float(y)
                                        duplicate[libname] =map(func , duplicate[libname] , [chrlibdup[ 1] ,chrlibdup [3 ]]) #modif 2011.11.23 
                                        duplicate[parse.sname] = map(func , duplicate[parse.sname] ,[ chrlibdup[ 1 ], chrlibdup[3 ]])
                        finally:
                                dup.close()
        unmap = [0 , 0 , 0 , 0]
        for sname in args :
                unmapPE = gzip.open (parse.outdir + "/"+os.path.basename(samname)+".PE.unmap.sam.gz" , "rb")
                
                try:
                        while True:
                                line = unmapPE.readline()
                                if not line:
                                        break
                                if int(line.split()[1]) & 0x40 :
                                        unmap[0] += 1 
                                        unmap[1] += len(line.split()[9])
                                else:
                                        unmap[2] += 1
                                        unmap[3] += len(line.split()[9]) 
                finally:
                        unmapPE.close()
    #    infoStat["total reads"] += unmap[0] + unmap[2]
     #   infoStat["total bases"] += unmap[1] + unmap[3]
     #   infoStat["read1"] += unmap[0]
     #   infoStat["read2"] += unmap[2]

        allInfo = open(parse.outdir+"/mapping.stat" , "w")  
        try:
          #      allInfo.writelines( "%15d total reads\n" %(infoStat[mapstat[0]]))
          #      allInfo.writelines( "%15d total bases\n" %(infoStat[mapstat[1]]))
          #      allInfo.writelines( "%15d read1\n" %(infoStat[mapstat[2]]))
           #     allInfo.writelines( "%15d read2\n" %(infoStat[mapstat[3]]))
         #       allInfo.writelines( "%15d mapped reads(%.2f%%)\n" %(infoStat[mapstat[4]] ,infoStat[mapstat[4]]/infoStat[mapstat[0]]*100))
          #      allInfo.writelines( "%15d mapped bases(%.2f%%)\n" %(infoStat[mapstat[5]] ,infoStat[mapstat[5]]/infoStat[mapstat[1]]*100))
           #     allInfo.writelines( "%15d PE mapped reads[with itself and mate mapped](%.2f%%)\n" %(infoStat[mapstat[6]],infoStat[mapstat[6]]/infoStat[mapstat[0]]*100))
            #    allInfo.writelines( "%15d PE bases(%.2f%%)\n" %(infoStat[mapstat[7]],infoStat[mapstat[7]]/infoStat[mapstat[1]]*100))
            #    allInfo.writelines( "%15d SE mapped reads[singletons](%.2f%%)\n" %(infoStat[mapstat[8]],infoStat[mapstat[8]]/infoStat[mapstat[0]]*100))
             #   allInfo.writelines( "%15d SE bases(%.2f%%)\n" %(infoStat[mapstat[9]],infoStat[mapstat[9]]/infoStat[mapstat[1]]*100))
             #   allInfo.writelines( "%15d properly paired reads(%.2f%%)\n" %(infoStat[mapstat[10]],infoStat[mapstat[10]]/infoStat[mapstat[0]]*100))
             #   allInfo.writelines( "%15d properly paired bases(%.2f%%)\n" %(infoStat[mapstat[11]],infoStat[mapstat[11]]/infoStat[mapstat[1]]*100))
              #  allInfo.writelines( "%15d trimed reads(%.2f%%)\n" %(infoStat[mapstat[12]],infoStat[mapstat[12]]/infoStat[mapstat[0]]*100))
             #   allInfo.writelines( "%15d trimed bases(%.2f%%)\n" %(infoStat[mapstat[13]],infoStat[mapstat[13]]/infoStat[mapstat[1]]*100))
             #   allInfo.writelines( "%15d insertion envent reads(%.2f%%)\n" %(infoStat[mapstat[14]],infoStat[mapstat[14]]/infoStat[mapstat[0]]*100))
              #  allInfo.writelines( "%15d insertion bases(%.2f%%)\n" %(infoStat[mapstat[15]],infoStat[mapstat[15]]/infoStat[mapstat[1]]*100))
              #  allInfo.writelines( "%15d deletion envent reads(%.2f%%)\n" %(infoStat[mapstat[16]],infoStat[mapstat[16]]/infoStat[mapstat[0]]*100))
              #  allInfo.writelines( "%15d deletion bases(%.2f%%)\n" %(infoStat[mapstat[17]],infoStat[mapstat[17]]/infoStat[mapstat[1]]*100))
              #  allInfo.writelines( "%15d with mate mapped to a different chr\n" %(infoStat[mapstat[18]]))
              #  allInfo.writelines( "%15d with mate mapped to a different chr (mapQ>=5)\n" %(infoStat[mapstat[19]]))
               # allInfo.writelines( "%15d QC failure \n" %(infoStat[mapstat[20]]))
              #  allInfo.writelines( "%14.2f%% duplicates \n" %(duplicate[parse.sname][0]/duplicate[parse.sname][1]*100))
                duplicate.pop(parse.sname)
                for libname in lib: 
                        print  str(duplicate[libname][0])+"       "+str(duplicate[libname][1])
                        allInfo.writelines("%18s duplicates:%.2f%%\n" %(libname,duplicate[libname][0]/duplicate[libname][1]*100))        
        finally:
                allInfo.close()
