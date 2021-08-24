#!/usr/bin/env python3
import os
import sys
import glob
from icecream import ic


def get_jobs(cycles):
    jobs = []
    for c in cycles:
        qid = c.split(".")[-1].replace("o", "")
        jobs.append(qid)
    return jobs

def check_status(f, sign):
    with open(f, 'r') as fh:
        info = fh.read()
    if sign in info:
        return True
    else:
        return False

def testcmd(command):
    import subprocess
    ret = subprocess.run(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding="utf-8",timeout=1)
    if ret.returncode == 0:
        return True
    else:
        return False


def jobsrun(job):
    #user =  os.system('whoami')
    command = 'qstat -j ' + str(job)
    result = testcmd(command)
    return result


def main():
    regx = sys.argv[1]
    sign = sys.argv[2]
    files = glob.glob(regx)

    for f in files:
        cycles = glob.glob(f + ".o*")
        # not run yet
        if len(cycles) < 1:
           print(f"perhaps {f} is not run")
           os.system(f"sel-qsub evo 5g 1 {f}")
        else:
            jobs = get_jobs(cycles)
            running = []
            done = 0
            for j in jobs:
                stdout = f + ".o" + j
                if check_status(stdout, sign): # job runs good and finished
                    #print(f"{f} is done by qid {j}")
                    done += 1
                else:
                    # job is not done
                    if jobsrun(j): # still running
                        #print(f"{f} is running by qid {j}")
                        running.append(j)
                    else: # not running, and no good sign
                        print(f"{f} run failed in cycle {j}")
                        #os.system("rm {f}.e{j} {f}.o{j}")

            if len(running) > 1:
                print(f"BAD, shell {f} is running under multi jobs: {running}")

            if len(running) == 0 and done == 0:
                print(f"{f} run failed")
                #os.system("qsub-sel evo 5g 1 {f}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(f"python3 {sys.argv[0]} regx|\*.sh sign|done")
    else:
        main()
