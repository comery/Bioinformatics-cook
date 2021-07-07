#/usr/bin/env python3
import sys
import os
import time
import subprocess

if len(sys.argv) < 2:
    sys.exit("python3 {} data.list [cp|ln|mv](optional)".format(sys.argv[0]))

def main():
    date = time.strftime("%Y-%m-%d", time.localtime())
    log_file = sys.argv[1] + "_" + date + ".log"
    log = open(log_file, 'w')

    if len(sys.argv) > 2:
        access_method = sys.argv[2]
        if access_method not in ['ln', 'cp', 'mv']:
            sys.exit("can not get data by {} way!".format(access_method))
    else:
        access_method = "ln -s "
    print("[INFO]: transfer method: " + access_method)
    workdir = os.getcwd()

    with open(sys.argv[1], 'r') as fh:
        for i in fh:
            if i.startswith("#"):
                continue
            tmp = i.strip().split()
            sample = tmp[0]
            barcodes = tmp[1].split(",")
            raw_path = tmp[2]
            none, dirname = os.path.split(raw_path)
            if len(tmp) >3 and os.path.isdir(tmp[3]) == True:
                raw_path = tmp[3] # splited data path instead of original path
                print("[NOTE]: {} it will get access data from {}".format(sample, raw_path), file=log)

            if os.path.exists(sample) == False:
                os.mkdir(sample)
            for b in barcodes:
                read1 = "{}/{}_{}_1.fq.gz".format(raw_path, dirname, b)
                read2 = "{}/{}_{}_2.fq.gz".format(raw_path, dirname, b)
                if os.path.exists(read1) == False or os.path.exists(read2) == False:
                    print("barcode {} of sample {} is not valid!".format(b, sample), file=log)
                else:
                    cmd = "{} {}/{}_{}_*.fq.gz {}/{}/".format(access_method, raw_path, dirname, b, workdir, sample)
                    subprocess.call(cmd, shell=True)
    log.write("All done")
    log.close()

if __name__ == '__main__':
    main()


