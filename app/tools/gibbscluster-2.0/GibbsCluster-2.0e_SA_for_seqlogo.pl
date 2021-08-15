#!/usr/bin/perl -w
#Stand-alone version
#Massimo Andreatta - Nov 2016
#
# 2.0 version implements insertions and deletions in the alignment
# NOTE: slow mode has been discontinued from this version

use strict;
use Getopt::Std;			# Command line options
use FileHandle;
use File::Spec;
my $bin;
BEGIN {
        $bin = File::Spec->rel2abs($0); #absolute path of caller
        $bin = readlink($bin) if -l $bin;  #called from symlink?
        $bin =~ s/(.*)(\/\S+)/$1/;
}
use lib "$bin/lib";   #add library to path

use Parallel::ChildManager;

my $call = join " ", $0, @ARGV;

########################
# Get command line options. Input file is mandatory (in format Pep || Pep TAB Annotation)
use vars qw($opt_b $opt_c $opt_f $opt_g $opt_h $opt_i $opt_j $opt_k $opt_l $opt_n $opt_p $opt_q $opt_r $opt_s $opt_t $opt_u $opt_z $opt_B $opt_C $opt_D $opt_G $opt_H $opt_I $opt_P $opt_R $opt_S $opt_T $opt_X $opt_W);
getopts('b:c:f:g:hi:j:k:l:n:pq:r:s:t:u:z:B:CD:G:H:I:P:R:S:TXW:');

unless (defined($opt_f)) {
    usage();
    
    print "ERROR. No training data uploaded\n<br>";
	exit;
}
usage() if $opt_h;

my $bdir="$bin/bin"; #binary directory
my @warnings; ##array of warning messages

################################################################
### Platform type
my $platform = qx(uname -sm);
chomp $platform;

my $subbin;

if ($platform eq "Linux x86_64") {
    $subbin = "Linux_x86_64";
} elsif ($platform eq "Darwin x86_64") {
    $subbin = "Darwin_x86_64";
} else {
    print "Platform: $platform not supported.";
    exit(0);
}

###############################
## Results directory (assign a unique name using the variable $$) 

my $pwd = qx(pwd); 
chomp $pwd;
my $resdir=$pwd;     ##results go to current directory, unless specified
$resdir="$opt_R" if defined($opt_R);

# my $sID=$$;  #session ID
# $sID++ while(-e "$resdir/$sID");  ##find a non-used session ID
# my $prefix=$sID;
# if (defined($opt_P)) {
# 	my $P=$opt_P;
# 	$P =~ s/\s+/_/g;  ##remove problematic characters
# 	$P =~ s/\\/_/g;
# 	$P =~ s/\//_/g;
# 	$P =~ s/[><]/_/g;
# 	$prefix = $P . "_$sID";  #use custom prefix for files
# }
# $resdir .= "/$prefix";

#my $pathbase = $resdir;
my $pathbase = "."; ##use relative paths for html display


my $c="mkdir -p $resdir/logos; mkdir $resdir/images; mkdir $resdir/cores; mkdir $resdir/matrices; ";
$c .= "mkdir $resdir/data; mkdir $resdir/res";
# system($c)==0 or die "Failed to execute $c: $?\n";  #failed to create results directory

my $timeanddate = localtime();

print "#Call: $call\n";
print "#$timeanddate\n";
# print "#Session ID: $sID\n";
# print "#Run name: $prefix\n";

#### options for Clustering algorithm

my $cluster_mode=0;
$cluster_mode = 1 if defined($opt_C);
my $lgt = 9;
$lgt = $opt_l if defined($opt_l);
my $dlen = 0;
$dlen = $opt_D if defined($opt_D);
my $ilen = 0;
$ilen = $opt_I if defined($opt_I);
my $iter = 10;
$iter = $opt_i if defined($opt_i);
my $seed = 1;
$seed = $opt_S if defined($opt_S);
my $ts = 1.5;
$ts = $opt_t if defined($opt_t);
my $nt = 20;
$nt = $opt_n if defined($opt_n);
my $fs = 100;
$fs = $opt_s if defined($opt_s);
my $fr = 20;
$fr = $opt_r if defined($opt_r);
my $findels = 10;
$findels = $opt_u if defined($opt_u);
my $p1a = 0;
$p1a = 1 if defined($opt_p);
my $lambda = 0.8;
$lambda = $opt_b if defined($opt_b);
my $swt = 0;
$swt=$opt_c if defined($opt_c);
my $blosum = "62";
$blosum = $opt_B if defined($opt_B);
my $wlc = 200;
$wlc = $opt_W if defined($opt_W);
my $bg=1;
$bg = $opt_z if defined($opt_z);
my $baseline=0;
$baseline=$opt_j if defined($opt_j);
my $sigma = 5;
$sigma = $opt_q if defined($opt_q);
my $trash = 0;
$trash = 1 if defined($opt_T);
my $procs = 1;
$procs = $opt_k if defined($opt_k);

my $indels = ($ilen>0 || $dlen>0) ? 1 : 0; 

my $ming = 1;
my $maxg = 5;
if (defined ($opt_g)) {
    if ($opt_g =~ m/(\d+)-(\d+)/) {
	if ($2>$1) {
	    $maxg=$2;
	    $ming=$1;
	} else {
	    $maxg=$1;
	    $ming=$2;
	}
    } elsif ($opt_g =~ m/(\d+)/) {
	$maxg=$1;
	$ming=$1;
    } else {
  
	&error_die("Wrong format. Number of clusters must be a number (e.g. 3) or an interval (e.g. 1-7)");
    }
}

###Some error checks
unless (isIntPos($iter)) {
	&error_die("Option -i must be an integer positive number");
}
unless (isIntPos($lgt)) {
	&error_die("Option -l must be an integer positive number");
}
unless (isIntPos($seed)) {
	&error_die("Option -S must be an integer positive number");
}
unless (isInt($wlc)) {
	&error_die("Option -W must be an integer number");
}
unless (isIntPos($procs)) {
	&error_die("Option -k must be an integer positive number");
}
unless (isIntPos($nt)) {
	&error_die("Option -n must be an integer positive number");
}
unless (isInt($ilen)) {
	&error_die("Option -I must be an integer number");
}
unless (isInt($dlen)) {
	&error_die("Option -D must be an integer number");
}
unless (isIntPos($fs)) {
	&error_die("Option -s must be an integer positive number");
}
unless (isIntPos($fr)) {
	&error_die("Option -r must be an integer positive number");
}
unless (isIntPos($findels)) {
	&error_die("Option -u must be an integer positive number");
}
unless ($bg eq "0" || $bg eq "1" || $bg eq "2") {
	&error_die("Option -z must be either 0, 1 or 2");	
}
unless ($swt eq "0" || $swt eq "1" || $swt eq "2") {
	&error_die("Option -c must be either 0, 1 or 2");	
}
unless (isAnumber($ts)) {
	&error_die("Option -t must be a number");
}
unless (isAnumber($lambda)) {
	&error_die("Option -b must be a number");
}
unless (isAnumber($sigma)) {
	&error_die("Option -q must be a number");
}
unless (isAnumber($baseline)) {
	&error_die("Option -j must be a number");
}

my $bgmode = "Uniprot pre-calculated";
$bgmode = "Flat" if $bg==0;
$bgmode = "From data" if $bg==2; 

my $totjobs = $seed*($maxg-$ming+1);

## R and seq2logo
my $R;
my $seq2logo;

$seq2logo = $opt_G if defined($opt_G);
$R = $opt_H if defined($opt_H);

### Data directory
my $ddir="$bin/data";

### Executables etc.

my $gibbsC = "$bdir/gibbs_cluster_1.2.2_DB." . $subbin;

my $Rmakebarplot = "$bdir/R-KLDvsClust_barplot.R";

## child manager
my $cm=new ChildManager($procs) || die "Died. Cannot create ChildManager object\n";

##############################################
## read input file, save to convenient format into local file

my $rawdata_unsorted = "$resdir/data/pepdata_us.txt";
my $rawdata = "$resdir/data/pepdata.txt";
my $warnflag = 0;
my $ln=0;
my $unique=0;
my $tooshort=0;
my $noP1rem=0;
my %flag;
# open (IN,'<',$opt_f) or die "Cant open file $opt_f: $!";
# open (OUT,'>',$rawdata_unsorted) or die "Canot create file $rawdata_unsorted: $!";
# while (defined(my $l=<IN>)) {
#    chomp $l;
#    next if length($l)<=1;  #don't kill on empty lines

#    my $pep;
#    my $annotation;
   
#    if ($l =~ /^>/) {  ##skip Fasta headers
#       next;
#    }
   
#    if ($l =~ /(\S+)\s+(\S+)/) {
#        $pep=$1;
#        $annotation=$2;
#    }
#    elsif ($l =~ /(\S+)/) {
#        $pep=$1;
#        $annotation="Peplist";
#    }
   
#     $pep =~ tr/a-z/A-Z/;
#     if ($pep=~m/[^ACDEFGHIKLMNPQRSTVWY]/) {
#          &error_die("Wrong format at line $ln (unknown amino acid):\n$l"); 
#     }

#     $ln++;
    
#     if (length($pep)+$ilen<$lgt) { #too short
# 		$tooshort++;
# 		next;
#     }
#     next if exists($flag{$pep});
#     $flag{$pep}=1;
#     $unique++;

#    if ($cluster_mode==1 && $indels==0) {  #pre-aligned and no indels
# 	    my $tmp = length($pep); 
# 	    if ($tmp != $lgt) {		 
# 		 	$warnflag = 1; 
# 	    }
#    }
#     print OUT "$pep\t$annotation\n";

# }
# close IN;
# close OUT;

## Sort alphabetically (for reproducibility)
# system("cat $rawdata_unsorted | sort -k1 > $rawdata; rm $rawdata_unsorted");

%flag=();

# &error_die("No sequences submitted") if $ln==0;
#&error_die("Too few input sequences. Give at least $Wminseq") if $ln<$Wminseq;
# &error_die("All sequences are shorter than the motif length. Select a shorter length") if $unique==0;
#&error_die("Not enough sequences longer than the motif length. Select a shorter length") if ($ln>=$Wminseq && $unique<$Wminseq);

if ($warnflag==1) {
	push(@warnings, "Running in pre-aligned mode, but sequence length is different from specified core length (option -l). You probably want to align the sequences (-C flag off) or allow insertions or deletions (-D and/or -I > 0)");
}

# print "#Read $unique unique sequences";
# print $tooshort>0 ? " (removed $tooshort shorter than motif length)\n" : "\n";
# print "##### Settings: #####\n";

# print $cluster_mode==1 ? "#No shift moves, cluster move at every iteration\n" : "#Shift moves and cluster moves are activated\n";
print "#Number of clusters: ";
print $ming!=$maxg ? "$ming - $maxg\n" : "$maxg\n";
# print "#Motif length: $lgt\n";
# print "#Initial MC temperture: $ts\n";
# print "#Number of temperature steps: $nt\n";
# print "#Number of iterations x Sequence x Tstep: $iter\n";
# print "#Max insertion length: $ilen\n";
# print "#Max deletion length: $dlen\n";
# print "#Interval between Indel moves: $findels\n" if ($ilen>0 || $dlen>0);
# print "#Interval between Single Peptide moves: $fr\n" if ($cluster_mode==0);
# print "#Interval between Phase Shift moves: $fs\n" if ($cluster_mode==0);
# print "#Number of initial seeds: $seed\n";
# print "#Penalty lambda: $lambda\n";
# print "#Weight on small clusters: $sigma\n"; 
# print "#Preference for hydrophobic P1: $p1a\n" if ($cluster_mode==0);
# print "#(removed $noP1rem seq not containing allowed P1 aa)\n" if $noP1rem>0;
# print "#Sequence weighting type: $swt\n";
# print "#Use trash cluster to remove outliers: $trash\n";
# print "#Threshold for trash cluster: $baseline\n" if $trash==1;

#####################################################
## run CLUSTERING ALGORITHM

# print "#### Running Gibbs clustering...\n";
my $cmd;
# for (my $g=$ming; $g<=$maxg; $g++) {
#     for (my $s=1; $s<=$seed; $s++) {
    	 
# 		 my $clust_outfile = "$resdir/res/gibbs.${g}g.${s}s.out";
# 		 my $corefile = "$resdir/cores/gibbs.${g}g.${s}s.cores";
# 		 my $alnfile = "$resdir/aln/gibbs.${g}g.${s}s.aln";

# 		 $cmd = "$gibbsC ";
# 		 $cmd .= "-mhc1 " if $cluster_mode == 1;
# 		 $cmd .= "-p1a " if $p1a == 1;
# 		 $cmd .= "-trash " if $trash == 1;
# #		 $cmd .= "-fast " if $fast_mode == 1;
# 		 $cmd .= "-l $lgt -g $g -ts $ts -blf $ddir/blosum%i.freq_rownorm -blm $blosum -base $baseline ";
# 		 $cmd .= "-bg $bg -i $iter -s $s -nt $nt -fr $fr -fs $fs -wlc $wlc -swt $swt -lambda $lambda ";
# 		 $cmd .= "-dlen $dlen -ilen $ilen -findels $findels ";
# #		 $cmd .= "-sigma $sigma -pc -c $corefile -a $alnfile $rawdata > $clust_outfile";
# 		 $cmd .= "-sigma $sigma $rawdata > $clust_outfile";
# 		 $cm->start($cmd);
		 
# #		 print "###Gibbs call: $cmd\n";
		 
# 		 print "#Clustering with $g groups, seed number $s\n";
#     }
# }
# $cm->wait_all_children;

# print "#Clustering complete!\n";
# print "### Determining seeds with highest KLD...\n";

### find best seed for each length
# for (my $g=$ming; $g<=$maxg; $g++) {
#     my $bkldS=-999;
#     my $bS=1;
#     for (my $s=1; $s<=$seed; $s++)
#     {
# 		  my $clust_outfile = "$resdir/res/gibbs.${g}g.${s}s.out";
# 		  open (RES,'<',$clust_outfile) or die "Cannot open file $clust_outfile: $!";
# 		  while (defined(my $l=<RES>)) {
# 			  chomp $l;
# 			  if ($l =~ m/iterations: (\S+)/) {
# 			      if ($1>$bkldS) {
# 				     $bkldS=$1;
# 				     $bS=$s;
# 			  	  }
# 	          }
# 	      }
# 	      close RES;
#      }

# 	 print "#Best $g groups, seed $bS\n";     
#      #only keep the best seed for each $g
     
#      $cmd = "cp $resdir/res/gibbs.${g}g.${bS}s.out $resdir/res/gibbs.${g}g.out; ";
# #     $cmd .= "cp $resdir/cores/gibbs.${g}g.${bS}s.cores $resdir/cores/gibbs.${g}g.cores; ";   
# #    $cmd .= "cp $resdir/aln/gibbs.${g}g.${bS}s.aln $resdir/aln/gibbs.${g}g.aln";
#      system($cmd);
# }

####################################
## extract relevant Info from the output files
my (%kld,%size,%totseq,%kldAVG,%nonemptyG,%trashsize);
my (%kldprint,%kldAVGprint);
my ($bestG,$bestKLD)=(1,-999);
my %cores;

# print "###Parsing result files...\n";

# for (my $g=$ming; $g<=$maxg; $g++) {
#     my $clust_outfile = "$resdir/res/gibbs.${g}g.out";
#     my $clust_outfile_display = "$resdir/res/gibbs.${g}g.ds.out";
#     $totseq{$g}=0;
#     $nonemptyG{$g}=0;
#     open (RES,'<',$clust_outfile) or die "Cannot open file $clust_outfile: $!";
#     open (DIS, '>', $clust_outfile_display) or die "Cannot make file $clust_outfile_display: $!";
#     my $header="G Gn  Num           Sequence       Core o of ip  IP il  IL dp  DP dl  DL Annotation sS   Self bgG bgG bgS bgScor cS cScore\n";
# 	print DIS $header;
    
#     while (defined(my $l=<RES>)) {
		
# 		chomp $l;
	
# 		if ($l =~ m/^G\s+\d+\s+\d+\s+.*\s+ip\s+/) {
# 			my @e = split(' ', $l);
# 			my $group = $e[1];
# 			my $core = $e[4];
# 			push(@{$cores{$g}{$group}}, $core);
# 			print DIS "$l\n";
# 		} elsif ($l =~ m/^#Trash/) {
# #			print DIS "$l\n";
# 		}
# 		if ($l =~ m/group (\d+), size (\d+).+epochs - (\S+)/) {
# 			$kld{$g}{$1+1}=$3;
# 			$size{$g}{$1+1}=$2;
# 			$totseq{$g}+=$2;
# 			$kldprint{$g}{$1+1}=sprintf("%.3f",$3); 
# 		} elsif ($l =~ m/iterations: (\S+)/) {
# 			$kldAVG{$g}=$1;
# 			$kldAVGprint{$g}=sprintf("%.3f",$1);
# 			if ($kldAVG{$g}>$bestKLD) {
# 			$bestKLD=$kldAVG{$g};
# 			$bestG=$g;
# 			}
# 		} elsif ($l =~ m/removed (\d+) outliers/) {
# 			$trashsize{$g}=$1;
# 		} elsif ($g==$ming && $l =~ m/not contain allowed P1/) {  ##only on first file
# 			$noP1rem++;
# 		}
#     }
#     close RES;
#     close DIS;
# }

########################################################
# Save KLD vs. number of groups to file, for a barplot of cluster quality
# my $data4barplot="$resdir/images/gibbs.KLDvsClusters.tab";
# open (TAB,'>',$data4barplot) or die "Cannot create file $data4barplot: $!\n";
# print TAB "1";
# for (my $g=2; $g<=$maxg; $g++) {
#     print TAB "\t$g";
# }
# print TAB "\n";
# for (my $g=$ming; $g<=$maxg; $g++) {
#     my $kk = $kldAVG{$g};
#     print TAB "$g";
#     for (my $p=1; $p<=$maxg; $p++) {
# 	if ($p<=$g) {
# 	   # my $kld = exists($kld{$g}{$p}) ? $kld{$g}{$p} : 0;
# 	    my $size = exists($size{$g}{$p}) ? $size{$g}{$p} : 0;
# 	    $nonemptyG{$g}++ if $size>0;
# 	    my $val=$kk*$size/$totseq{$g};  #relative cluster weight
# 	    $val = 0 if $val<0;
# 	    print TAB "\t$val";
#         } else {
# 	    print TAB "\t0";
# 	}
#     }
#     print TAB "\n";
# }
# close TAB;

# make a barplot
# my $barplot = "$resdir/images/$prefix.gibbs.KLDvsCluster.barplot.pdf";
# my $barplot = "$resdir/images/gibbs.KLDvsCluster.barplot.JPG";

# if (defined($R)) {  ##only if R is installed
# 	$cmd = "cat $Rmakebarplot | $R --vanilla --args $data4barplot $barplot > $resdir/images/R.log";
# 	system($cmd);
# } 

###################################################################
#### put cores into individual files, for logo-generation

# for (my $g=$ming; $g<=$maxg; $g++) {
# open (C,'>',"$resdir/cores/gibbs.${bestG}g.cores") or die "Cannot create core file: $!\n";
# if (exists($cores{$bestG})) {
# 	my %levels = %{$cores{$bestG}};
# 	for (my $lev=0; $lev<=$maxg; $lev++) {
# 		if (exists($levels{$lev})) {
# 			my @seqs = @{$levels{$lev}};
# 			my $cg = $lev+1;
# 			print C "## Alignment cores for group $cg\n";	
# 			open (OUT,'>',"$resdir/cores/gibbs.${cg}of${bestG}.core") or die "Cannot create core file: $!\n";
# 			for (my $s=0; $s<=$#seqs; $s++) {
# 				my $core = $seqs[$s];
# 				print OUT "$core\n";
# 				print C "$core\n";
# 			}
# 			close OUT;
# 		}
# 	}
# }
# close C;
# }

#######################################
## seq2logo
if (defined($seq2logo)) {
#    for (my $g=$ming; $g<=$maxg; $g++) {
       for (my $p=1; $p<=$maxg; $p++) {
			my $logofile = "$resdir/logos/gibbs_logos_${p}of${maxg}";
			my $corefile = "$resdir/cores/gibbs.${p}of${maxg}.core";
			# my $matfile = "$resdir/matrices/gibbs.${p}of${bestG}.mat";
			my $title = quotemeta("Group $p of $maxg");

	   		next unless (-e $corefile);

			$cmd = "cd $resdir/logos; ";  #seq2logo has very short path limitation for the output file.
			
	   		$cmd .= "$seq2logo -f $corefile -o $logofile -I 2 --format [JPEG] -b $wlc -C 2 -t $title &>/dev/null; ";

	   		$cmd .= "cd $pwd; sleep 2";
			   
	   		system($cmd);
	   		
	   		##copy the Log-odds matrix made by Seq2Logo
	   		# $cmd = "cat $logofile.txt | grep -v '#' > $matfile";
			
	   		# system($cmd);
       }
#    }
}

######################################
### Show cluster statistics in plain text

# for (my $g=$ming; $g<=$maxg; $g++)
# {
#    print "#RESULTS for $g CLUSTERS\n";
#    print "#$g Final Average KLD: $kldAVG{$g}\n";
#    for (1..$g) {
# 	   my $p=$_;
# 	   print "#$g\t$p\t$size{$g}{$p}\t$kldprint{$g}{$p}\n";
#    }
#    if ($trash==1) {
#    		print "#Outliers: $trashsize{$g}\n#\n";
#    }
# }

#########################################
### Show warnings to log

if ($#warnings>-1) {
	print "###THERE ARE WARNINGS:\n";
	
	for (0..$#warnings) {
		my $w = $warnings[$_];
		print "#$_" . ": $w\n";
	}	
}

###############################
##### make HTML output

# my $html = "$resdir/$prefix" . "_report.html";

# my $courier = "<font face=\"COURIER\">";

# open (H,'>',$html) or die "Cannot create html file $html: $!"; 

# print H "<!DOCTYPE html>\n";
# print H "<html>\n<head><title>Gibbs Cluster Report</title>\n";

# print H "<style type=\"text/css\">\n";

# print H "body { \nfont-family: Arial, Helvetica, sans-serif;\ncolor: #000000;background-color: #FFFFFF;\n";
# print H "font-size: 15px; border: 0;\nmargin: 0;\npadding: 0;\n}\n";

# print H "div.header {\ncolor: #FFFFFF;\nbackground-color: #FF3399;\nfont-size: 200%;\nfont-weight: bold;\n";
# print H "border:0;\nmargin:0;\npadding: 0.5em;\n}\n";

# print H "div.footer {\ncolor: #FFFFFF;\nbackground-color: #FF3399;\nfont-size: 130%;\nfont-weight: bold;\n";
# print H "border:0;\nmargin:0;\npadding: 0.5em;\n}\n";

# print H "div.main {\nmargin-left: 2em;\nmargin-right: 2em;\n}\n";

# print H "a:link {\ncolor: #FF3399\n}\n";
# print H "a:visited {\ncolor: #FF3399\n}\n";
# print H "a:hover {\ncolor: #642EFE\n}\n";
# print H "a:active {\ncolor: #642EFE\n}\n";

# print H "</style>\n</head>\n";
# print H "<body>";

# print H "<div class=\"header\"> GibbsCluster report</div>\n";
# print H "<div class=\"main\">\n";

# print H "<br><b>Version:</b> 2.0\n<br>";
# print H "<b>Run ID:</b> $sID\n<br>";
# print H "<b>Run name:</b> $prefix\n<br>" if $sID ne $prefix;
# print H "<b>Platform:</b> $platform\n<br>";

# print H "<br><br>Read <b>$unique unique sequences</b> from file";
# print H $tooshort>0 ? " (removed $tooshort shorter than motif length)<br><br>\n" : "<br><br>\n";
# print H "<b>Settings:</b><br><small>\n";

# ##show options
# print H $cluster_mode==1 ? "<b>No shift moves, cluster moves at every iteration</b><br>\n" : "<b>Shift moves and cluster moves activated</b><br>\n";
# print H "Number of clusters: ";
# print H $ming!=$maxg ? "$ming - $maxg<br>\n" : "$maxg<br>\n";

# print H "Motif length: $lgt<br>\n";
# print H "Initial MC temperture: $ts<br>\n";
# print H "Number of temperature steps: $nt<br>\n";
# print H "Number of iterations x Sequence x Tstep: $iter<br>\n";
# print H "Max insertion length: $ilen<br>\n";
# print H "Max deletion length: $dlen<br>\n";
# print H "Interval between Indel moves: $findels<br>\n" if ($ilen>0 || $dlen>0);
# print H "Interval between Single Peptide moves: $fr<br>\n" if ($cluster_mode==0);
# print H "Interval between Phase Shift moves: $fs<br>\n" if ($cluster_mode==0);
# print H "Number of initial seeds: $seed<br>\n";
# print H "Penalty lambda: $lambda<br>\n";
# print H "Weight on small clusters: $sigma<br>\n";
# print H "Preference for hydrophobic P1: $p1a<br>" if ($cluster_mode==0);
# print H " (removed $noP1rem seq not containing allowed P1 aa)<br>" if $noP1rem>0;
# print H "Sequence weighting type: $swt<br>\n";
# print H "Background model: $bgmode<br>\n";
# print H "Use trash cluster to remove outliers: $trash<br>\n"; 
# print H "Threshold for trash cluster: $baseline<br>\n" if $trash==1;

# print H "</small>\n";

# print H "<b>Install Seq2Logo to visualize the sequence motifs!</b><br>\n" unless defined($seq2logo);

# print H "</font><hr>";
# ## show results on two columns
#  print H "<table width=\"100%\" border=\"0\">\n";
#   print H "<tr valign=\"top\"><td width:50%;text-align:center>\n<div style=\"padding-left: 2em\">\n";

# print H defined($R) ? "<b>KLD vs. Number of clusters with &#955; = $lambda</b></div>\n" : "</div>\n";

#  print H "</td> <td width:50%;text-align:center>\n<div style=\"padding-left: 2em\">";
# print H defined($seq2logo) ? "<b>Identified $nonemptyG{$bestG} sequence motifs</b>\n</div>\n" : "</div>";
# print H "</td></tr>";

# print H "<tr valign=\"middle\"><td width:50%;text-align:center>\n";
# ##barplot

# if (defined($R)) {
# 	print H "<img src=\"$pathbase/images/$prefix.gibbs.KLDvsCluster.barplot.png\" alt=\"\" title=\"Barplot of combined Kullbach Leibler Distance versus the Number of Partitions applied to the dataset. \" width=\"400\" />\n<br>";
# 	print H "<span style=\"padding-left:5em\"><font size=\"2\">View the <a href=\"$pathbase/images/$prefix.gibbs.KLDvsCluster.barplot.png\" TARGET='_blank'>barplot</a> in full size</font></span>\n<br>";
# }

# ##Logos with the optimal number of clusters
# print H "</td> <td width:50%;text-align:left;vertical-align:middle>\n";

# my $width=200;
# if ($nonemptyG{$bestG}==1) {
# 	$width=300;
# } elsif ($nonemptyG{$bestG}==2) {
# 	$width=250;
# }
# my $cnt=0;
# for (my $nlogo=1; $nlogo<=$bestG; $nlogo++) {
# 	if (-e "$resdir/logos/gibbs_logos_${nlogo}of${bestG}-001.png") {
# 			$cnt++;
# 			print H "<img src=\"$pathbase/logos/gibbs_logos_${nlogo}of${bestG}-001.png\" alt='' title=\"LOGO $nlogo of $bestG - KLD=$kldprint{$bestG}{$nlogo}\" width=\"$width\" />";
# 			print H "<br>\n" if ($cnt%3==0 && $nlogo<$bestG);  #many logos. go to a new line every 3
# 	}
# }	
# print H "</td></tr></table>\n<br>\n<hr>\n";

# print H "<table width=\"100%\" border='0'>\n";

# for (my $g=$ming; $g<=$maxg; $g++)
# {
#    print H "<tr valign=\"top\"><td colspan='2'>\n";

#    print H "<h3>RESULTS for $g CLUSTERS</h3>\n</td></tr>";

#    print H "<tr valign='top'> <td> $courier\n";
#    print H " Final Average KLD: <b>$kldAVG{$g}\n</b></font><br><br>\n";

#    print H "<table align='left'>\n";
#    print H "$courier<tr valign='top'> <th>&nbsp;Group&nbsp;</th> <th>&nbsp;Size&nbsp;</th> <th>&nbsp;KLD&nbsp;</th> <th>&nbsp;Seq2Logo&nbsp;</th> <th>&nbsp;Matrix&nbsp;</th> </tr>\n";

#    for (1..$g)
#    {
# 	   my $p=$_;
# 	   print H "<tr valign='top'> <td>&nbsp;$p&nbsp;</td> <td>&nbsp;$size{$g}{$p}&nbsp;</td> <td>&nbsp;$kldprint{$g}{$p}&nbsp;</td>\n<td>";
# 	##transfer peptides to Seq2logo Webserver
		
# 		my $cores = "$resdir/cores/gibbs.${p}of${g}.core";
# 		if (-e $cores)  
# 		{
# 		   print H "<div style='display: inline;'>";
# 		   print H "<form style='display:inline; margin:0;' action='http://www.cbs.dtu.dk/biotools/Seq2Logo/index.php' method='POST' target='_blank'> <input type='hidden' name='INPASTE' value='";

	   
# 		   open( C, "<", $cores ) || die "Can't open $cores: $!";	   
# 		   my $text = join('', <C>);
# 		   close C;
# 		   print H $text;

# 		   print H "'/> <input type='submit' value='LOGO' /></form></div>\n";
# 	   }
# 	   ## Link to matrix file
# 	   if (-e "$resdir/matrices/gibbs.${p}of${g}.mat") {
# 		   print H "<td>&nbsp; <a href=\"$pathbase/matrices/gibbs.${p}of${g}.mat\" TARGET='_blank'>Mat_$p.$g</a> &nbsp;</td>\n";
# 	   }
   
#    		print H " </td>\n</tr>";
#    }
#    if ($trash==1) {
# 	 print H "<tr valign='top'> <td>&nbsp;Outliers&nbsp;</td> <td>&nbsp;$trashsize{$g}&nbsp;</td> <td>&nbsp; &nbsp;</td>\n <td>&nbsp; &nbsp;</td></tr>";
   
#    }
#    print H "</font>";
#    print H "<tr valign='top'> <td colspan='4'><br>Raw <a href=\"$pathbase/res/gibbs.${g}g.out\" TARGET='_blank'>Clustering Report</a>\n";
#    print H "<br>Formatted <a href=\"$pathbase/res/gibbs.${g}g.ds.out\" TARGET='_blank'>Clustering Solution</a>\n";
#    print H "<br>Clustered <a href=\"$pathbase/cores/gibbs.${g}g.cores\" TARGET='_blank'>Alignment Cores</a></td></tr>\n";

# 	print H "</table>\n";

# 	print H "</td> <td width:50%;text-align:left;vertical_align:middle>\n";

# 	my $cnt=0;
# 	 for (my $nlogo=1; $nlogo<=$g; $nlogo++) {
# 		 if (-e "$resdir/logos/gibbs_logos_${nlogo}of${g}-001.png") {
# 			 $cnt++;
# 			 print H "<img src=\"$pathbase/logos/gibbs_logos_${nlogo}of${g}-001.png\" alt='' title=\"LOGO $nlogo of $g - KLD=$kldprint{$g}{$nlogo}\" width=\"150\" />";
# 			 print H "<br>\n" if ($cnt%3==0 && $nlogo<$g);  #many logos. go to a new line every 3
# 		 }
# 	 }
# 	 print H "</td></tr>\n<tr><td colspan=\"2\"></td></tr>\n";
# }
# print H "</table>";

# print H "\n<br><hr><br>Explain the <a href=\"http://www.cbs.dtu.dk/services/GibbsCluster/output.php\" target=\"_blank\">output</a>\n<br><br>\n";
# print H "Read the <a href=\"http://www.cbs.dtu.dk/services/GibbsCluster/instructions.php\" target=\"_blank\">instructions</a></div><br><br>";
# print H "\n<br><div class=\"footer\"> Made with GibbsCluster version 2.0</div>\n";
# print H "</body></html>\n";

# close H;

exit(0);

sub usage{
  print("$0\n");
  print(".\n");
  print("\n");
  print("Usage: $0 -f peptidelist [-g number of clusters] [-h] [more options]\n");
  print("Command line options:\n");
  print("  -f       upload training set (Peptide || Peptide TAB Annotation)\n");
  print("  -P       name for this run (no spaces)\n");
  print("  -g (1-5) number of clusters (as single number [3] or interval [1-5])\n");
  print("  -l (9)   motif length\n");
  print("  -i (10)  number of iterations x sequence x T step\n");
  print("  -t (1.5) start MC temperature\n");
  print("  -n (20)  number of temperature steps\n");
  print("  -b (0.8) lambda penalty on inter-cluster similarity\n");
  print("  -q (5)   sigma weight on small alignments\n");
  print("  -T       use trash cluster to remove outliers (switch)\n");
  print("  -j (0)   threshold for trash cluster\n");
  print("  -C       perform Single sequence moves at every iteration (default is every -r iterations)\n");
  print("  -u (10)  interval between Indel moves\n");
  print("  -r (20)  interval between Single sequence moves\n");
  print("  -s (100) interval betweem Phase shift moves\n");
  print("  -D (0)   max length of deletions\n");
  print("  -I (0)   max length of insertions\n");
  print("  -p       turn on prefence for hydrophobic at P1 (switch)\n");
  print("  -S (1)   number of seeds for initial starting conditions\n");
  print("  -c (0)   sequence weighting: [0] 1/ns [1] clustering [2] none\n");
  print("  -z (1)   background frequencies: [0] flat [1] pre-calculated [2] from data\n");
  print("  -W (200)  weigth on low count\n");
  print("  -G       path to Seq2Logo\n");
  print("  -H       path to R\n");
  print("  -k (1)   number of processes for multi-threading\n");
  print("  -R       results directory (defaults to current directory)\n");
  print("  -h       print this message\n");

  exit;
}

## error message
sub error_die {
    my $message = $_[0];

    print "ERROR\n$message\n"."Died.\n";
    exit(0);
}

## check format of a number
sub isAnumber {
    my $test = shift;

    eval {
        local $SIG{__WARN__} = sub {die $_[0]};
        $test += 0;
    };
    if ($@) {
	return 0;}
    else {
	return 1;} 
}

## check for integer
sub isInt {
    my $test=shift;
    
    if ($test =~ m/^\d+$/) {
	return 1; }
    else {
	return 0; }
}

## check for positive integer
sub isIntPos {
    my $test=shift;
    
    if ($test =~ m/^\d+$/ && $test>0) {
	return 1; }
    else {
	return 0; }
}
