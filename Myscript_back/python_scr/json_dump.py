import json
def store_to_db(data, outfile):
    with open(outfile, 'w') as fw:
        json.dump(data,fw, indent=4)


data = {'a':1,
        'b':2}

store_to_db(data, "test.json")
