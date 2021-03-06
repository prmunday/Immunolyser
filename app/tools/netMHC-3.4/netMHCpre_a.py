#! /usr/local/python/bin/python

# NetMHC 3.4 main script
# 
# Morten Nielsen mniel@cbs.dtu.dk, March 2013
# Center for Biological Sequence Analysis


###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

###############################################################################
#               NOTHING SHOULD NEED CHANGING BELOW THIS LINE!
###############################################################################


import sys, os, math, tempfile, datetime, time, copy

from optparse import OptionParser

# Path to the 'netMHC-3.4' directory (mandatory)
#preddir = "/usr/cbs/bio/src/netMHC-3.4";
preddir = os.environ["NMHOME"] 
#print os.environ["NMHOME"]
#print preddir

version = "3.4"

#tmpdir        = tempfile.mkdtemp(dir='/var/tmp')
#tmpdir         = tempfile.mkdtemp(dir='/scratch')
tmpdir         = tempfile.mkdtemp(dir=os.environ["TMPDIR"])

# Build commandline parser
parser = OptionParser(usage="usage: %prog [options] file", version=version)

parser.add_option("-a", "--mhc", type="string", dest="alleles", metavar="STR",
                  help="Allele names ( ',' -separated)")

parser.add_option("-l", "--peplen", type="int", dest="peplen", metavar="NUM",
                  help="Length of subpeptides to predict")

parser.add_option("-x", "--xls", type="string", dest="xlsname", metavar="STR",
                  help="Name of tab separated output file")

parser.add_option("-s", "--sort", action="store_true", dest="sort",
                  help="Sort output on descending affinity")

parser.add_option("-p", "--peptide", action="store_true", dest="pepfile",
                  help="infile is in peptide format")
		  
parser.add_option("-n", "--nodirect", action="store_true", dest="nodirect",
                  help="Do not use direct prediction (use 9mer aproximation)")

parser.add_option("-b", "--noblacklist", action="store_true", dest="noblacklist",
                  help="Do not use blacklist")

parser.add_option("-A", "--Alleles", action="store_true", dest="Showall",
                  help="Show available alleles and exit")



parser.set_defaults(alleles='', xlsname='', peplen=9, sort=False, pepfile=False, Showall=False, nodirect=False, noblacklist=False)

(opts, args) = parser.parse_args()

if not opts.Showall:
	infilename=args[0]

infohandle = os.popen('uname -s')
platform   = infohandle.readline().rstrip()
infohandle.close()
infohandle = os.popen('uname -m')
platform   = platform + "_" + infohandle.readline().rstrip()
infohandle.close()

#preddir = os.environ.get("NETMHCIII")
os.environ["NETMHCIII"] = preddir
#wwwrun                  = os.environ.get("USER") == None
wwwrun                  = os.environ.get("USER") == "www"

#sys.stderr.write("### %s\n" % wwwrun)
#print os.environ.get("user")
#print os.environ.get("USER")
#sys.exit()

sys.path.append(preddir+'/lib/python/') ### add local lib to lib path 

import sequences ### use local module sequences

NN_listname   = preddir+'/etc/NN.list'
pydate        = datetime.date.today()
date          = pydate.strftime("%A %B %e %Y")
outtime       = time.strftime("%R")
    
whitelistfilename = "%s/etc/whitelist" % preddir
blacklistfilename = "%s/etc/blacklist" % preddir

if wwwrun:
	outputdir           = tempfile.mkdtemp(dir='/usr/opt/www/pub/CBS/services/NetMHC-3.4/tmp/')
	dummy, xlfilename   = tempfile.mkstemp(dir=outputdir, suffix='.xls')
	#os.system('chmod a+rwx %s' % ( outputdir ))
else:
	xlfilename        = opts.xlsname
	
dummy, tmpfsafilename  = tempfile.mkstemp(dir=tmpdir)
dummy, pepfilename     = tempfile.mkstemp(dir=tmpdir)
dummy, predfilename    = tempfile.mkstemp(dir=tmpdir)

os.system('chmod a+rwx %s' % ( tmpdir ))
os.system('chmod g+r %s && chmod g+r %s' % (pepfilename, predfilename))

def showalleles(NNfile, tdir):
    allist = []
    print "### Alleles with ANN predictors:"
    for line in NNfile:
        fields = line.split()
        allist.append(fields[0])
    NNfile.close()
    allist.sort()
    for al in allist:
        print al
        
    os.system('rm -rf %s/*' % tdir)
    os.system('rm -rf %s/.??*' % tdir)
    try:
        os.rmdir(tdir)
    except:
        sys.stderr.write("COULD NOT REMOVE TEMPORARY DIRECTORY %s\nPROBABLY ON MOUNTED FILE SYSTEM!\n\n" % tdir) 
    sys.exit()
    return()
        

def check_alleles(NNfile, allist):
    newlist = []
    NN_alleles = []
    for line in NNfile:
        fields = line.split()
        NN_alleles.append(fields[0])
    NNfile.close()
    for al in allist:
        if al in NN_alleles:
            newlist.append(al)
        else:
            sys.stderr.write("ERROR!!\nNO PREDICTION METHODS AVAILABLE FOR ALLELE %s\nFOR A LIST OF AVAILABLE ALLELES USE OPTION -A\n" % al)
            sys.exit()
    return(newlist)


def build_whitelist(filename, nodirect):
    try:
        whitelistfile = open(filename, 'r')
    except:
        print "COULD NOT OPEN %\n" % filename
        sys.exit()

    whitelist = {}
    if not nodirect:
     for line in whitelistfile:
         fields = line.split()
	 mhc = fields[0]
         if not whitelist.has_key(mhc):
             whitelist[mhc] = []
         whitelist[fields[0]].append(int(fields[1]))
    whitelistfile.close()
    return(whitelist)


def pep2fsa(pepfile, fastafile):
    seq = ""
    for line in pepfile:
        fields = line.split()
        seq += fields[0]
	plen = len(fields[0])
    pepfile.close()
    fastafile.write(">Peptides\n%s\n" % seq)
    fastafile.close()
    return(plen)

def get_thresholds(file):
    thrs = {}
    for line in file:
        fields = line.split()
        thrs[fields[0]] = fields
    file.close()
    return(thrs)

def readseqs(filename):
    file=open(filename,'r')
    entries, order = sequences.readfsa(file)
    file.close()
    return(entries, order)

def seq2pepfile(seq, filename,plen,pepinp):
    pepfile = open(filename, 'w')
    if pepinp:
        for i in  range(0, len(seq)-1, plen):
            pepfile.write('%s\n' % seq[i:i+plen])
    else:
        for pos in range(len(seq)-(plen-1)):
            pepfile.write('%s\n' % seq[pos:pos+plen])
    pepfile.close()
    return()

def predict(NN, infilename, tdir, pdir, platform, mhc, outfilename, plen, wl):

    usedirect = False
    
    if wl.has_key(mhc):
        if plen in wl[mhc]:
            usedirect = True
            mhc = mhc+"_"+str(plen)+"mer"

    if plen>9:
        if not usedirect:
            w = plen - 9
            dummy, newpepfilename = tempfile.mkstemp(dir=tdir)
            npepfile = open(newpepfilename, 'w')
            infile = open(infilename, 'r')
            n = 0
            for line in infile:
                pep = line.rstrip()
                for i in range(6):
                    npep = pep[:3+i]+pep[3+i+w:]
                    npepfile.write("%s %s %d\n" % (npep, pep, n))
                n+=1
            npepfile.close()
            infilename = newpepfilename
            
    elif plen==8:
        if not usedirect:
            dummy, newpepfilename = tempfile.mkstemp(dir=tdir)
            npepfile = open(newpepfilename, 'w')
            infile = open(infilename, 'r')
            n = 0
            for line in infile:
                pep = line.rstrip()
                for i in range(3,8):
                    npep = pep[:i]+"X"+pep[i:]
                    npepfile.write("%s %s %d\n" % (npep, pep, n))
                n+=1
            npepfile.close()
            infilename = newpepfilename
    else:
        usedirect = True

    if NN:
    	#print "### NN ###"
	#sys.exit()
        bindir      = "%s/bin/%s/" % (pdir, platform)
        syndir      = pdir+"/etc/net/" + mhc
        blsynlist   = syndir + "/bl50/synlist"
        spsynlist   = syndir + "/sparse/synlist"
         
        gawkcmd = "gawk '{if($1 == %s#%s) {print($0)} else{ printf(%s%s%s%s%s%s,$0,exp(log(50000)-($NF*log(50000))))}}'" % ('"', '"', '"', '%s', '\\t', '%d', '\\n', '"')

	nihpred = "%s/epipred3.0 -dirty -v -bdir %s -tdir %s -seq2inp seq2inp -blsyn %s -spsyn %s  %s | %s" % ( bindir, bindir, tdir, blsynlist, spsynlist, infilename, gawkcmd) 

        cmd = '%s | grep -v "#" > %s' % (nihpred, outfilename)

    try:
        os.system(cmd)
    except:
        print 'Could not execute command %s\n' % cmd
        sys.exit(1)
	#os.system('rm -f %s' % infilename)
    return(mhc, usedirect)

def predfile2predlist(filename, NN, plen, apprx, ave):
    list=[]
    predfile  = open(filename, 'r')
    predfile2 = predfile
    if ave:
	    predfile = open(filename+"_2", 'r')
    nl = 0
    s = 0.0
    for line in predfile:
        predorder = []
        fields=line.split()
        if NN:
            if plen>9 and apprx:
                s = s + float(fields[4])
                nl += 1
                if nl == 6:
		    s = s/6.0
		    if ave:
			    myline = predfile2.readline()
			    myfields = myline.split()
			    s = (s + float(myfields[3]))/2.0
                    aff = math.exp(math.log(50000.0)*(1.0-s))
                    predorder = [s, fields[2], fields[3], aff]
                    nl = 0
                    s = 0.0
            elif (plen == 8) and apprx:
                s = s + float(fields[4])
                nl += 1
                if nl == 5 :
		    s = s/5.0
		    if ave:
			    myline = predfile2.readline()
			    myfields = myline.split()
			    s = (s + float(myfields[3]))/2.0
                    aff = math.exp(math.log(50000.0)*(1.0-s))
                    predorder = [s, fields[2], fields[3], aff]
                    nl = 0
                    s = 0.0
            else:
                predorder = [fields[3], fields[1], fields[0], fields[4]] # Predorder 0:log50k, 1:peptide, 2:position, 3:affinity(nM)
        else:
            if plen>9:
                s = s + float(fields[3])
                nl += 1
                if nl == 6 :
                    s = s/6.0
                    predorder = [s, fields[2], fields[0]]
                    nl = 0
                    s = 0.0
            elif plen == 8:
                s = s + float(fields[3])
                nl += 1
                if nl == 5 :
                    s = s/5.0
                    predorder = [s, fields[2], fields[0]]
                    nl = 0
                    s = 0.0
            else:
                predorder = [float(fields[2]), fields[1], fields[0]] # Predorder 0:score, 1:peptide, 2:position
        if len(predorder)>0:
            list.append(predorder)
    predfile.close()
    predfile2.close()
    return(list)

def set_thr(NN, thresholdlist, mhc):
    if NN:
        thr1 = 50
        thr2 = 500
    else:
        thresfile = open(thresholdlist, 'r')
        thresholds = get_thresholds(thresfile)
        dummy, thr1, thr2 = thresholds[mhc]
    return(thr1, thr2)

def printresults(NN, order, preds, allele, methodlist, plen, www):
    thr1 = 50
    thr2 = 500
    thr1string = str(thr1)+" nM"
    thr2string = str(thr2)+" nM"
    method = 'Artificial Neural Networks'
    if www:    
        downloadstring = '<a href="%s">Download output sheet</a>' % xlfilename[20:]
    else:
        downloadstring = ''
    xlsdict       = {}
    webheader         = "%4s %10s %13s %12s %10s %15s %7s" % ('pos', 'peptide', 'logscore', 'affinity(nM)', 'Bind Level', 'Protein Name', 'Allele')


    
    method = method + " - " + methodlist
    print downloadstring
    print "\nNetMHC version %s. %smer predictions using %s. Allele %s. \nStrong binder threshold %6s. Weak binder threshold score %6s\n\n%s\n" % (version, plen, method, allele, thr1string, thr2string, downloadstring)
    
    print "----------------------------------------------------------------------------------------------------"
    print "%s" % webheader
    print "----------------------------------------------------------------------------------------------------"

    for entry in order:
        xlsdict[entry] = []
        if NN:
            for  (score, pep, pos, aff) in preds[entry]:
                bindlevel = ''
                if float(aff)<thr2:
                    bindlevel = 'WB'
                    if float(aff)<thr1:
                        bindlevel = 'SB'
                try:
			xlsdict[entry].append([int(pos), pep, aff, float(score)/0.64])
		except:
			print pos
			print pep
			print aff
			print score
			sys.exit()
                print "%4i %10s %13.3f %12d %10s %15s %7s" % (int(pos), pep, float(score), int(aff), bindlevel, entry, allele)
        xlsdict[entry].sort()
        print "--------------------------------------------------------------------------------------------------"
    print downloadstring
    return(xlsdict)
    
        

def printxls(fulldict, allist, xlout, methlst, plen):
    xlout.write("NetMHC version %s\tLength of predicted peptides:%d\tDate:\t%s\ttime:\t%s\n\n" % (version, plen, date, outtime))
    #xlout.write("Analyzed File:\t%s\n\n" % opts.infilename)
    xlout.write("Protein\tPosition\tPeptide")
    averages = {}
    for hla in allist:
        method = methlst[hla]
        xlout.write("\t%s ANN/Mat %s predicted affinity (Kd, nM)/Matscore" % (hla, method))
    #xlout.write("\tAverage score (higher score = stronger affinity)\n")
    xlout.write("\n")
    for entry in fulldict[allist[0]]:
        #print entry
        averages[entry] = []
        for pepline in fulldict[allist[0]][entry]:
            #print pepline
            pos = int(pepline[0])
            pep = pepline[1]
	    NN  = pepline[-1]
            xlout.write("%s\t%d\t%s" % (entry, pos, pep))
            
            for hla in allist:
                affscore = float(fulldict[hla][entry][pos][2])
                score    = float(fulldict[hla][entry][pos][3])
                if hla == allist[0]:
                    averages[entry].append(score)
                else:
                    averages[entry][pos] += score
                    
                if NN:
                    xlout.write("\t%d" % (int(affscore)))
                else:
                    xlout.write("\t%.2f" % (float(affscore)))
            #xlout.write("\t%.3f\n" % (float(averages[entry][pos]/float(len(allist)))))
            xlout.write("\n")
    xlout.write("\n")
    xlout.close()
    return()

    
######### Main #################
def main():

    if opts.Showall:
        showalleles(open(NN_listname, 'r'), tmpdir)
    if infilename == '':
        sys.stderr.write("ERROR!\n\nMISSING INFILE (option -f filename)\n\n")
        sys.exit()
    else:
        try:
            dummy = open(infilename, 'r')
            dummy.close()
        except:
            sys.stderr.write("ERROR!\n\nCOULD NOT OPEN INPUT FILE: %s\n\n" % infilename)
            sys.exit()
        
    alleles = []
    if not opts.alleles == '':
        alleles = opts.alleles.split(",")
    if len(alleles)==0:
            sys.stderr.write("ERROR!\n\nYOU MUST SELECT AT LEAST ONE ALLELE!\n\n")        
            sys.exit(1) 
    alleles    = check_alleles(open(NN_listname, 'r'), alleles)
    
    peplen     = opts.peplen
    preds      = {}
    badchars   = "BJOUZ"
    pinp       = opts.pepfile # peptide input
    
    whitelist  = build_whitelist(whitelistfilename, opts.nodirect)
    blacklist  = build_whitelist(blacklistfilename, opts.noblacklist)
    
    if pinp:    
        if not wwwrun:
	   peplen = pep2fsa(open(infilename,'r'), open(tmpfsafilename, 'w'))
	   fastafilename = tmpfsafilename
	else:
	   fastafilename = infilename
    else:
        fastafilename = infilename
    
    seqs, inorder = readseqs(fastafilename)
    fulldict={}
    predtypes={}
    print "%s %s" % (date, outtime)
    for entry in seqs:
        seq=seqs[entry].seq.upper()
        for char in badchars:
            if char in seq:
                sys.stderr.write("\nERROR: Illegal character '%s' in sequence %s \n\n" % (char, entry))
                sys.exit()

    methodlist = {}
    for allele in alleles:
        if allele[-3:]=='mat':
            NNpred = False
        else:
            NNpred = True

        useapprox = True

        if blacklist.has_key(allele) and peplen in blacklist[allele]:
            useapprox = False
            
        for entry in seqs:
            seq=seqs[entry].seq.upper()
            seq2pepfile(seq, pepfilename, peplen, pinp)
            allele2, direct = predict(NNpred, pepfilename, tmpdir, preddir, platform, allele, predfilename, peplen, whitelist )
	    ave = False
	    if (not peplen == 9) and direct and useapprox:
		    direct = False
		    ave    = True
		    mypredfilename = predfilename + "_2"
		    dummy1, dummy2 = predict(NNpred, pepfilename, tmpdir, preddir, platform, allele, mypredfilename, peplen, {})
		
            approxpred = (not direct)
            if approxpred:
                if ave:
                    methodlist[allele] = "Approx+Direct"
                else:
                    methodlist[allele] = "Approximation" 
            else:
                 methodlist[allele] = "Direct"

            predtypes[allele] = NNpred
            preds[entry] = predfile2predlist(predfilename, NNpred, peplen, approxpred, ave)
            if opts.sort:
                preds[entry].sort()
                preds[entry].reverse()

        fulldict[allele] = printresults(NNpred, inorder, preds, allele, methodlist[allele], peplen, wwwrun)
        os.system('rm -f %s %s' % (pepfilename, predfilename))
    if not xlfilename == '':
        #sys.stderr.write("### %s\n" % xlfilename)
        printxls(fulldict, alleles, open(xlfilename, 'w'), methodlist, peplen)
    return()
        
            
main()


os.system("rm -rf %s/*" % tmpdir)
os.system("rm -rf %s/.??*" % tmpdir)
try:
    os.rmdir(tmpdir)
except:
    sys.stderr.write("COULD NOT REMOVE TEMPORARY DIRECTORY %s\nPROBABLY ON MOUNTED FILE SYSTEM!\n\n" % tmpdir) 
