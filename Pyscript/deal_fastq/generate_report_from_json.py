import sys
import json
from icecream import ic

def load_from_db(handle):
    data = json.load(handle)
    return data

sample = sys.argv[1].replace(".json", "").split("-")[0]

raw_reads, raw_bases, clean_reads, clean_bases, dup_rate = 0, 0, 0, 0, 0
number = 0

for f in sys.argv[1:]:
    number += 1
    with open(f, 'r') as fh:
        json_dict = load_from_db(fh)
        raw_reads += int(json_dict['summary']['before_filtering']['total_reads'])
        raw_bases += int(json_dict['summary']['before_filtering']['total_bases'])
        clean_reads += int(json_dict['summary']['after_filtering']['total_reads'])
        clean_bases += int(json_dict['summary']['after_filtering']['total_bases'])
        #clean_q20 += int(json_dict['summary']['after_filtering']['q20_rate'])
        #clean_q30 +=  int(json_dict['summary']['after_filtering']['q30_rate'])
        dup_rate += float(json_dict['duplication']['rate'])

average_duprate = "{:.4f}".format(dup_rate/number)
filtered_rate = "{:.4f}".format(1-(clean_reads/raw_reads))
print(f"{sample}\t{raw_reads}\t{raw_bases}\t{clean_reads}\t{clean_bases}\t{filtered_rate}\t{average_duprate}")
