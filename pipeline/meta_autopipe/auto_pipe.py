#!/hwfssz1/ST_DIVERSITY/PUB/USER/yangchentao/software/bin/python3
import os
import sys
import pwd
import time
import json
import argparse
import subprocess
t = time.time()

description = ""
parser = argparse.ArgumentParser(description=description)

parser.add_argument('-l', dest='list', required=False, metavar='<str>',
                    help='sample and data list, no need when run in second time')

parser.add_argument('-o', dest='outdir', required=False, default="./",
                    metavar='<str>',
                    help='outdir, default = ./')

parser.add_argument('-c', dest='conf', required=False, metavar='<str>',
                    help='configure file')

parser.add_argument('-t', dest='interval', required=False, default=20,
                    metavar='<int>',
                    help='set interval time of checking by qstat, default 20 mins')

parser.add_argument('-d', dest='database', required=True, metavar='<str>',
                    default="TEST",
                    help='database name for saving the job work status, default=TEST')

args = parser.parse_args()

def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return 0
    else:
        return 1


def print_time(info, hd):
    print(info + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), file=hd)


def parser_conf(file):
    print("[INFO]: load configures from " + file)
    conf = {}
    source = {}
    with open(file, 'r') as cf:
        for i in cf:
            if i.startswith("#"):
                continue
            elif "|" in i:
                tmp = i.strip().split("|")
                qsub_conf = tmp[1].split(",")
                source[tmp[0]] = qsub_conf
            elif "=" in i:
                tmp = i.strip().split("=")
                conf[tmp[0]] = tmp[1]
            else:
                print("skipping lines")
    if len(conf.keys()) == 0:
        print("you gave empty software list!")
        exit()
    if len(source.keys()) == 0:
        print("you gave empty qsub parameters, so use default value")
    if 'filter' not in source.keys():
        source['filter'] = ['evo', '2g', '1']
    if 'merge' not in source.keys():
        source['merge'] = ['evo', '0.5g', '1']
    if 'rmcondam' not in source.keys():
        source['rmcondam'] = ['evo', '4g', '6']
    if 'assembly' not in source.keys():
        source['assembly'] = ['super', '300g', '10']
    if 'mappping' not in source.keys():
        source['mapping'] = ['evo', '4g', '6']

    return conf, source


def store_to_db(data, outfile):
    with open(outfile, 'w') as fw:
        json.dump(data,fw)

def load_from_db(handle):
    data = json.load(handle)
    return data


def submit(step, wkdir, shell, source, conf):
    shelldir = os.path.dirname(shell)
    script = os.path.basename(shell)
    # cd shell dir to submit
    os.chdir(shelldir)
    #print(shelldir)
    if step == 'filter':
        #print("filter")
        cmd = "{} {} {} {} {}".format(conf['sel-qsub'],
                                      source['filter'][0],
                                      source['filter'][1],
                                      source['filter'][2],
                                      script)
        subprocess.call(cmd, shell=True)
    elif step == 'merge':
        cmd = "{} {} {} {} {}".format(conf['sel-qsub'],
                                       source['merge'][0],
                                       source['merge'][1],
                                       source['merge'][2],
                                       script)
        subprocess.call(cmd, shell=True)
    elif step == 'rmcondam':
        cmd = "{} {} {} {} {}".format(conf['sel-qsub'],
                                     source['rmcondam'][0],
                                     source['rmcondam'][1],
                                     source['rmcondam'][2],
                                     script)
        subprocess.call(cmd, shell=True)
    elif step == 'assembly':
        cmd = "{} {} {} {} {}".format(conf['sel-qsub'],
                                      source['assembly'][0],
                                      source['assembly'][1],
                                      source['assembly'][2],
                                      script)
        subprocess.call(cmd, shell=True)
    elif step == 'mapping':
        cmd = "{} {} {} {} {}".format(conf['sel-qsub'],
                                      source['mapping'][0],
                                      source['mapping'][1],
                                      source['mapping'][2],
                                      script)
        subprocess.call(cmd, shell=True)
    else:
        print("Unknow step name")
        exit()
    # get back
    os.chdir(wkdir)


def get_fq_dir(path):
    dirs = os.listdir(path)
    fqs = []
    for i in dirs:
        if i.endswith("_1.fq.gz"):
            fqs.append(i.replace("_1.fq.gz", ""))
    return fqs


def run_filter(ID, fqdirs, outdir, source, conf):
    check_files = []
    wkdir = os.getcwd()
    outdir = os.path.abspath(outdir)
    sample_outdir = outdir + "/" + ID
    filter_outdir = sample_outdir + "/" + "filter"
    if os.path.exists(sample_outdir) == False:
        os.mkdir(sample_outdir)

    if os.path.exists(filter_outdir) == False:
        os.mkdir(filter_outdir)
    file_subfix = 0
    for fqdir in fqdirs:          # one sample maybe contain many fastq dir
        for i in get_fq_dir(fqdir):
            file_subfix += 1
            fq1 = fqdir + "/" + i + "_1.fq.gz"
            fq2 = fqdir + "/" + i + "_2.fq.gz"
            shell = filter_outdir + "/" + ID + "_" + str(file_subfix) + ".filter.sh"
            check_files.append(filter_outdir + "/" + ID + "_" + str(file_subfix) + ".filter.ok")
            with open(shell, 'w') as fh:
                fh.write("{} filter -1 {} -2 {} -l 10 -q 0.1 -m 20 -n 0.001 -G -Q 2 ".format(conf['SOAPnuke'], fq1, fq2))
                fh.write("-f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG\n")
                fh.write("touch " + filter_outdir + "/" + ID + "_" + str(file_subfix) + ".filter.ok")

            submit('filter', wkdir, shell, source, conf)
    return check_files


def run_merge(ID, outdir, source, conf):
    wkdir = os.getcwd()
    outdir = os.path.abspath(outdir)
    sample_outdir = outdir + "/" + ID
    filter_outdir = sample_outdir + "/" + "filter"
    merge_outdir = sample_outdir + "/" + "merge"
    if os.path.exists(sample_outdir) == False:
        print("can not find {} for sample {}".format(sample_outdir, ID))
        exit()

    if os.path.exists(merge_outdir) == False:
        os.mkdir(merge_outdir)

    shell = merge_outdir + "/" + ID + ".merge.sh"
    with open(shell, 'w') as fh:
        cmd = "zcat " + filter_outdir + "/*_1.fq.gz | gzip - >" + merge_outdir + "/" + ID + "_1.fq.gz" + "\n"
        cmd += "zcat " + filter_outdir + "/*_2.fq.gz | gzip - >" + merge_outdir + "/" + ID + "_2.fq.gz" + "\n"
        cmd += "touch " + merge_outdir + "/" + ID + ".merge.ok"
        fh.write(cmd + "\n")

    submit('merge', wkdir, shell, source, conf)


def run_rmcontam(ID, outdir, source, conf):
    wkdir = os.getcwd()
    outdir = os.path.abspath(outdir)
    sample_outdir = outdir + "/" + ID
    merge_outdir = sample_outdir + "/" + "merge"
    rmcondam_outdir = sample_outdir + "/" + "rmcondam"
    if os.path.exists(sample_outdir) == False:
        print("can not find {} for sample {}".format(sample_outdir, ID))
        exit()

    if os.path.exists(rmcondam_outdir) == False:
        os.mkdir(rmcondam_outdir)
    shell = rmcondam_outdir + "/" + ID + ".rmcondam.sh"
    with open(shell, 'w') as fh:
        fq1 = merge_outdir + "/" + ID + "_1.fq.gz"
        fq2 = merge_outdir + "/" + ID + "_2.fq.gz"
        cmd  = "{} -x {} -q -1 {} -2 {} --very-fast-local -p 6 --un-conc-gz {}_unconc.fq.gz | samtools view -bFS - > {}.bam\n".format(conf['bowtie2'], conf['ref'], fq1, fq2, ID, ID)

        cmd += "touch " + rmcondam_outdir + "/" + ID + ".rmcondam.ok"
        fh.write(cmd + "\n")

    submit('rmcondam', wkdir, shell, source, conf)


def run_assembly(ID, outdir, source, conf):
    wkdir = os.getcwd()
    outdir = os.path.abspath(outdir)
    sample_outdir = outdir + "/" + ID
    rmcondam_outdir = sample_outdir + "/" + "rmcondam"
    assembly_outdir = sample_outdir + "/" + "assembly"
    if os.path.exists(sample_outdir) == False:
        print("can not find {} for sample {}".format(sample_outdir, ID))
        exit()

    if os.path.exists(assembly_outdir) == False:
        os.mkdir(assembly_outdir)
    shell = assembly_outdir + "/" + ID + ".assembly.sh"
    with open(shell, 'w') as fh:
        fq1 = rmcondam_outdir + "/" + ID + "_unconc.fq.1.gz"
        fq2 = rmcondam_outdir + "/" + ID + "_unconc.fq.2.gz"
        cmd = "{} -o {} --meta --pe1-1 {} --pe1-2 {} -k auto -t 10 -m 800\n".format(conf['spades'], assembly_outdir, fq1, fq2)
        cmd += "touch " + assembly_outdir + "/" + ID + ".assembly.ok"
        fh.write(cmd + "\n")

    submit('assembly', wkdir, shell, source, conf)


def run_mapping(ID, outdir, source, conf):
    wkdir = os.getcwd()
    outdir = os.path.abspath(outdir)
    sample_outdir = outdir + "/" + ID
    rmcondam_outdir = sample_outdir + "/" + "rmcondam"
    assembly_outdir = sample_outdir + "/" + "assembly"
    mapping_outdir = sample_outdir + "/" + "mapping"

    if os.path.exists(sample_outdir) == False:
        print("can not find {} for sample {}".format(sample_outdir, ID))
        exit()

    if os.path.exists(mapping_outdir) == False:
        os.mkdir(mapping_outdir)
    shell = mapping_outdir + "/" + ID + ".mapping.sh"
    with open(shell, 'w') as fh:
        fq1 = rmcondam_outdir + "/" + ID + "_unconc.fq.1.gz"
        fq2 = rmcondam_outdir + "/" + ID + "_unconc.fq.2.gz"
        # select fasta 
        cmd = "{} -ap {} -lgt 200 -o {}/{}_scaffolds_200.fa {}/scaffolds.fasta\n".format(conf['fastaKit'],
                ID, mapping_outdir, ID, assembly_outdir)
        # statistic 
        cmd += "{} {}/{}_scaffolds_200.fa\n".format(conf['fastx-stats'], mapping_outdir, ID)
        # bowtie2-build index
        cmd += "{} {}/{}_scaffolds_200.fa {}_scaffolds \n".format(conf['bowtie2-build'],
                mapping_outdir, ID, ID)
        # bowtie2
        cmd += "{} -x {}_scaffolds -q -1 {} -2 {} --very-sensitive-local -p 6 --un-conc-gz {}_unconc_assem.fq.gz".format(conf['bowtie2'], ID, fq1, fq2, ID)

        cmd += " | {} view -bFS - > {}_assem.bam\n".format(conf['samtools'], ID)
        # samtools sort
        cmd += "{} sort {}_assem.bam {}_assem_sort\n".format(conf['samtools'], ID, ID)
        # calculate depth 
        cmd += "{} --outputDepth {}_depth.txt  {}_assem_sort.bam \n".format(conf['caldep'], ID, ID)
        # samtools index
        cmd += "{} index {}_assem_sort.bam\n".format(conf['samtools'], ID)
        # samtools idxstats
        cmd += "{} idxstats {}_assem_sort.bam > {}.dat\n".format(conf['samtools'], ID, ID)

        cmd += "touch " + mapping_outdir + "/" + ID + ".mapping.ok"
        fh.write(cmd + "\n")

    submit('mapping', wkdir, shell, source, conf)


def check_ok(ID, step, args, checklist=None):
    outdir = os.path.abspath(args.outdir)
    sample_outdir = outdir + "/" + ID
    if step == 'filter':
        sign = 1
        for i in checklist:
            if os.path.exists(i) == False:
                sign = 0
                break
        if sign == 0:
            return False
        else:
            return True
    else:
        sign_file = sample_outdir + "/" + step + "/" + ID + "." + step + ".ok"
        if os.path.exists(sign_file) == True:
            return True
        else:
            return False


def get_username():
    return pwd.getpwuid(os.getuid())[0]


def make_maintain_shell(outfile, this_script_name, args):
    user = get_username()
    wkdir = os.getcwd()
    default_bashrc = "/home/{}/.bashrc".format(user)
    if os.path.exists(default_bashrc) == False:
        default_bashrc = "/home/yangchentao/.bashrc"
    # run maintain shell
    cmd = "#！/usr/bin/bash\n"
    cmd += "cd " + wkdir + "\n"
    cmd += "source " + default_bashrc + "\n"
    cmd += "python3 {} -t {} -o {}  -d {}".format(this_script_name,
                                                      args.interval,
                                                      args.outdir,
                                                      args.database)

    with open(outfile, 'w') as fw:
        fw.write(cmd)


def main():

    """
    samples,filter,merge,rmcondam,assembly,mapping
    ID1,0,0,0,0,0
    ID2,0,0,0,0,0
    ...

    note: 0 means not start, 1 means job sumbmitted, 2 means job finished.
    """
    wkdir = os.getcwd()
    DB = "./" + args.database + ".json"
    # first time to run this script
    if os.path.exists(DB) == False:
        print("[INFO]: load sample list from " + args.list)
        database = {}
        if hasattr(args, 'conf') == False:
            print("-c is necessary for first time to run this script!")
            exit()
        elif hasattr(args, 'list') == False:
            print("-l is necessary for first time to run this script!")
            exit()
        else:
            # read configures
            conf, source = parser_conf(args.conf)
            database['conf'] = conf
            database['source'] = source

        """ when the first to run main(), read samples list and run fist step"""
        # read list
        check_files_filter = {}
        work_status = {}
        samples = {}
        with open(args.list, 'r') as fh, open(DB, 'w') as db:
            print("#Sample\tFilter\tMerge\tRmcondam\tAssembly\tMapping", file=db)
            for i in fh:
                if i.startswith("#"):
                    continue
                tmp = i.strip().split()
                sample = tmp[0]
                fqdir = tmp[1]
                if sample not in samples.keys():
                    samples[sample] = [fqdir,]
                else:
                    samples[sample].append(fqdir)

        for s in sorted(samples.keys()):
            check_files = run_filter(s, samples[s], args.outdir, source, conf)
            check_files_filter[s] = check_files
            work_status[s] = [1, 0, 0, 0, 0] # init work status, all steps for each sample are not fininshed.
        database['check_files_filter'] = check_files_filter
        database['work_status'] = work_status
        # now database contains all things we need
        store_to_db(database, DB)
        maintain_shell = args.database + ".maintain.sh"
        make_maintain_shell(maintain_shell, sys.argv[0], args)
        # make the maintain shell executable
        subprocess.call("chmod 755 " + args.database + ".maintain.sh", shell=True)
        # user need to set crontab profile
        print("Now all filter job have been submitted, the following issue matters!\n")
        print("crontab profile:\n" + "*/{} * * * * {}/{}.maintain.sh\n".format(args.interval,
                                                                                wkdir,
                                                                                args.database))
        print("please type \"crontab -e\" in your termnal and the save the last line in there")



    else:
        # regularly check work status
        TMPDB = DB + ".tmp"
        DONE = DB + ".done"
        LOG = DB + ".log"
        dh = open(DONE, 'a')
        with open(DB, 'r') as db, open(LOG, 'a') as log:
            # reload database from json file
            database = load_from_db(db)
            # unfold these dict from database
            work_status = database['work_status']
            check_files_filter = database['check_files_filter']
            conf = database['conf']
            source = database['source']
            # remove crontab file
            if len(work_status.keys()) == 0:
                print("all your samples are done!")
                subprocess.call("crontab -r", shell=True)
                exit()
            samples = sorted(work_status.keys())
            for sample in samples:
                signals = work_status[sample]
                status_sign = sum(signals)
                if status_sign == 1:
                    # check the first step result
                    if check_ok(sample, 'filter', args, check_files_filter[sample]) == True:
                        print_time(sample + " filter done", log)
                        # filter finished, do next step
                        run_merge(sample, args.outdir, source, conf)
                        new_status = [2, 1, 0, 0, 0]
                        work_status[sample] = new_status

                elif status_sign == 3:
                    # check the second step result
                    if check_ok(sample, 'merge', args) == True:
                        print_time(sample + " merge done", log)
                        # merge finished, do next step
                        run_rmcontam(sample, args.outdir, source, conf)
                        new_status = [2, 2, 1, 0, 0]
                        work_status[sample] = new_status
                elif status_sign == 5:
                    # check the third step result
                    if check_ok(sample, 'rmcondam', args) == True:
                        print_time(sample + " rmcondam done", log)
                        # rmcondam finished, do next step
                        run_assembly(sample, args.outdir, source, conf)
                        new_status = [2, 2, 2, 1, 0]
                        work_status[sample] = new_status
                elif status_sign == 7:
                    # check the fourth step result
                    if check_ok(sample, 'assembly', args) == True:
                        print_time(sample + " assembly done", log)
                        # assembly finished, do next step
                        run_mapping(sample, args.outdir, source, conf)
                        new_status = [2, 2, 2, 2, 1]
                        work_status[sample] = new_status
                elif status_sign == 9:
                    # check the fifth step result
                    if check_ok(sample, 'mapping', args) == True:
                        print_time(sample + " mapping done", log)

                        # mapping finished, remove it from list
                        dh.write(sample + "\n")
                        del work_status[sample]
                else:
                    # error
                    print("Unknow code {} for {}".format(status_sign, sample))

        dh.close()
        database['work_status'] = work_status
        # now update the database 
        store_to_db(database, TMPDB)
        cmd = "mv {} {}".format(TMPDB, DB)
        subprocess.call(cmd, shell=True)



if __name__ == '__main__':
    if os.path.exists(args.outdir) == False:
        os.mkdir(args.outdir)
    main()
