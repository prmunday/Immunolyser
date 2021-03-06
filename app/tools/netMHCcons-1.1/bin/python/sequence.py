#
# Define sequence class
#
import re
class sequence:
  def __init__(self):
    self.name = ""
    self.seq = ""
    self.assign = ""
    self.desc = ""
  #
  # Length
  #    
  def length(self):
    return(len(self.seq))
  #
  # fasta
  #
  def readfasta(self,file):
    for line in file:
      fields=line.split()
      if len(line)>1:
        if fields[0][0] == ">":
          self.name = fields[0][1:] # skip >
          self.desc = re.sub(r"^[^\s]+\s","",line.strip())
        else:
          self.seq = self.seq  + fields[0]
  def writefasta(self,file):
    file.write(">%s %s\n" % (self.name, self.desc))
    start=0
    while start < len(self.seq):
      file.write("%s\n" % (self.seq[start:start+60]))
      start +=60
  #
  # how
  #
  def readhow(self,file):
    line = file.readline()
    fields=line.split()
    if len(line)>1:
      if line[0][0] == " ":
        length = int(fields[0])
        self.name = fields[1]
	self.desc = re.sub(r"^[^\s]+\s","",line.rstrip()[7:])
      while len(self.seq) < length:
        line = file.readline()
	line.strip()
        fields=line.split()
        self.seq = self.seq  + fields[0]
      while len(self.assign) < length:
        line = file.readline()
	line.strip()
        fields=line.split()
        self.assign = self.assign  + fields[0]
  def writehow(self,file):
    file.write("%6d %s %s\n" % (len(self.seq),self.name, self.desc))
    start=0
    while start < len(self.seq):
      file.write("%s\n" % (self.seq[start:start+60]))
      start +=60
    start=0
    if len(self.assign) == 0:
      self.assign = "."*len(self.seq)
    while start < len(self.seq):
      file.write("%s\n" % (self.assign[start:start+60]))
      start +=60
