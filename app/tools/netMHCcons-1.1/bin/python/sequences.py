
import sys


class sequence:
    def __init__(self):
        self.seq=""
        self.desc=""
        self.name=""
        self.asgn=""
    
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


#########  FASTA ###################


def readfsa(file):
    seqdict={}
    #seqdict['order']=[]
    name = ""
    for line in file:
        if line[0:1] ==">":
            seq=""
            desc=""
            name=""
            asgn=""
            fields = line.split()
            if fields[0]==">":
                name=fields[1]
                desc=fields[2:]
            else:
                name=fields[0][1:]
                desc=fields[1:]
            seqdict[name]=sequence()
            seqdict[name].desc=desc
            seqdict[name].name=name
        else:
            if name != "":
                fields=line.split()
                if len(fields)>0:
                    seqdict[name].seq=seqdict[name].seq+fields[0]
    
    return(seqdict)

def writefsa(outfile, entry):
    header = ">"+entry.name
    for desc in entry.desc:
        header += " "+desc
        
    outfile.write(header+"\n")
    splits = range(0, len(entry.seq), 80)
    i=0
    for aa in entry.seq:
        outfile.write(aa)
        i+=1
        if i in splits:
            outfile.write("\n")
    if not i in splits:
        outfile.write("\n")

def readfsaII(file):
    seqdict={}
    #seqdict['order']=[]
    name = ""
    for line in file:
        if line[0:1] ==">":
            seq=""
            desc=""
            name=""
            fields = line.split()
            if fields[0]==">":
                name=fields[1]
                desc=fields[2:]
            else:
                name=fields[0][1:]
                desc=fields[1:]
            seqdict[name]=sequence()
            seqdict[name].desc=desc
        else:
            if name != "":
                fields=line.split()
                if len(fields)==1:
                    seqdict[name].seq=seqdict[name].seq+fields[0]
    return(seqdict)


############ HOW ##########################

def readhow(infile):
     seqdict={}
     inentry = False
     inass   = False
     inseq   = False
     for line in infile:
        if not inentry:
            if len(line)>1 and line[0]==" ":
                inentry = True
                fields  = line.split()
		fields.reverse()
                howlen  = fields.pop()
         	howlen  = int(howlen)
                name    = fields.pop()
                desc    = fields
                seq     = ""
                inseq   = True
        elif inseq:
            fields = line.split()
            seq   += fields[0]
            if len(seq)==howlen:
                inass = True
                inseq = False
                ass   = ""
            elif len(seq)>howlen or not len(fields[0]) == 80:
                sys.stderr.write("ERROR !  SEQ LEN != DEF LEN (%d:%d)\n" %(len(seq), howlen))
                sys.exit()
        elif inass:
            fields = line.split()
            ass   += fields[0]
            if len(ass)==howlen:
                inass   = False
                inentry = False
                seqdict[name]      = sequence()
                seqdict[name].desc = desc
                seqdict[name].name = name
		seqdict[name].seq  = seq
                seqdict[name].asgn = ass
            elif len(ass)>howlen or not len(fields[0]) == 80:
                sys.stderr.write("ERROR !  ASSIGN LEN != DEF LEN (%d:%d)\n" %(len(ass), howlen))
                sys.exit()
     return(seqdict)

                
                
                
     

def writehow(file, seq):
    
    l          = len(seq.seq)
    linebreaks = range(79, l, 80)
    
    file.write("%6d %s\t" % (l, seq.name))
    for word in seq.desc:
        file.write("%s " % word)
    file.write("\n")
               
    
    for i in range(l):
        file.write("%s"% seq.seq[i])
        if i in linebreaks:
            file.write("\n")
    if not (l-1) in linebreaks:
        file.write("\n")
        
    for i in range(l):
        try:
            file.write("%s"% seq.asgn[i])
        except:
            file.write("-")
        if i in linebreaks:
            file.write("\n")
    if not (l-1) in linebreaks:
        file.write("\n")
    return()



################## TAB ###############################        

def readtab(file):
    seqlist = []
    for line in file:
        fields       = line.split('\t')
        thisseq      = sequence()
        thisseq.name = fields[0].strip()
        thisseq.seq  = fields[1].strip()
        thisseq.desc = fields[3].strip()
        seqlist.append(thisseq)
    file.close()
    return(seqlist)
