



class sequence:
    def __init__(self):
        self.seq=""
        self.desc=""
	self.name=""
    
class aln:
    def __init__(self):
        self.count=0
        self.len=0
        self.seqs={}
    
def readphyaln(infile):
    ali=aln()
    head=infile.readline()
    fields=head.split()
    ali.count=int(fields[0])
    ali.len=int(fields[1])
    for n in range(ali.count):
        line=infile.readline()
        fields=line.split()
        name=fields[0]
        ali.seqs[name]=sequence()
        ali.seqs[name].seq=fields[1]
        for f in range(2,len(fields)):
            ali.seqs[name].seq=ali.seqs[name].seq+fields[f]
        #print(name+'\t'+ali.seqs[name].seq)
    
    i=len(ali.seqs[name].seq)
    while i<ali.len:
        line=infile.readline()
        for name in ali.seqs:
            line=infile.readline()
            fields=line.split()
            for f in range(len(fields)):
                ali.seqs[name].seq=ali.seqs[name].seq+fields[f]
        i=len(ali.seqs[name].seq)
    return(ali)


        
def readphyalns(phyalnfile):
    alns=[]
    eof=False
    while not eof:
        try:
            alns.append(readphyaln(phyalnfile))
        except:
            eof=True
        #for name in alns[0].seqs:
         #   print(str(len(alns))+'\t'+name+'\t'+alns[0].seqs[name].seq)
    return(alns)



def readfsa(file):
    seqdict = {}
    inorder = []
    names   = {}
    i       = 0
    for line in file:
        if line[0:1] ==">":
            i += 1
            entry = "Entry_%d" % i
            seq   = ""
            desc  = ""
            name  = ""
            fields = line.split()
            if fields[0]==">":
                name = fields[1]
                desc = fields[2:]
                if name in names:
                    names[name] += 1
                    name = name + desc[0]
                    name = name[:10]
                    if name in names:
                        name = name + str(names[name])
                else:
                    names[name] = 0
            else:
                name=fields[0][1:]
                desc=fields[1:]
                if len(desc)>0:
                    name += desc[0]
                    name = name[:10]
            inorder.append(entry)
            seqdict[entry]      = sequence()
            seqdict[entry].desc = desc
	    seqdict[entry].name = name
        else:
            if name != "":
                fields=line.split()
                if len(fields)>0:
                    seqdict[entry].seq = seqdict[entry].seq+fields[0]
    return(seqdict, inorder)


def readfsa_mn(file):
    seqdict = {}
    inorder = []
    names   = {}
    i       = 0
    for line in file:
        if line[0:1] ==">":
            i += 1
            entry = "Entry_%d" % i
            seq   = ""
            desc  = ""
            name  = ""
            fields = line.split()
            if fields[0]==">":
                name = fields[1]
                desc = fields[2:]
            else:
                name=fields[0][1:]
                desc=fields[1:]

            if len(desc)>0:
                name = name + "_" + desc[0]
                name = name[:15]

	    if name in names:
                names[name] += 1
                name = name[:13] + "_" + str(names[name])
            else:
                names[name] = 0

            inorder.append(entry)
            seqdict[entry]      = sequence()
            seqdict[entry].desc = desc
	    seqdict[entry].name = name
        else:
            if name != "":
                fields=line.split()
                if len(fields)>0:
                    seqdict[entry].seq = seqdict[entry].seq+fields[0]

    return(seqdict, inorder)

