with open("test.fq",'r') as fh:
    id = fh.readline().strip()
    while id:
        seq = fh.readline().strip()
        fh.readline().strip()
        qul = fh.readline().strip()
        print("id= " + id + "\n" +\
                "seq= " + seq + "\n" +\
                "qul= " + qul + "\n")
        id = fh.readline().strip()

