#! /usr/bin/perl -w

# Author: Edita Karosiene, edita@cbs.dtu.dk, 03 Sept 2014

use strict;

use Getopt::Long;
use Env;
use Sys::Hostname;

# specify version of the method
my $version = "1.1";

# variables used to save options
my ($method, $rdir, $NetMHC, $NetMHCpan, $PickPocket, $n_pep, $n_bind, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W, $dist, $inptype, $length, $tdir, $alleles, $a1, $a2, $a3, $a4, $training_MHC, $training_Pan_Pick, $file, $dirty, $hlaseq, $sort, $xls, $xlsfile, $filter, $rank_f, $aff_f, $w, $help, $v);

# associating option values with variables

&GetOptions (
'rdir:s' => \$rdir,
'mhc:s' =>   \$NetMHC,
'pan:s'=>   \$NetMHCpan,
'pick:s'=>   \$PickPocket,
'method:s'=>\$method,
'np=s' => \$n_pep,
'nb=s' => \$n_bind,
'affS=s' => \$thr_aff_S,
'affW=s' => \$thr_aff_W,
'rankS=s' => \$thr_rank_S,
'rankW=s' => \$thr_rank_W,
'dist=s' => \$dist,
'inptype=i' => \$inptype,
'length:s' => \$length,
'xlsfile:s' => \$xlsfile,
'tdir:s' => \$tdir,
'trMHC:s' => \$training_MHC,
'trPanPick:s' => \$training_Pan_Pick,
'filter:s' => \$filter,
'rankF=s' => \$rank_f,
'affF=s' => \$aff_f,
'a:s' => \$alleles,
'a1:s' => \$a1,
'a2:s' => \$a2,
'a3:s' => \$a3,
'a4:s' => \$a4,
'f:s'=>   \$file,
'hlaseq:s' => \$hlaseq,
's' => \$sort,
'xls' => \$xls,
'w' => \$w,
'dirty' => \$dirty,
'v' => \$v,
'h' => \$help
);

# define platform
my $PLATFORM = $ENV{'PLATFORM'};

#define tmp directory
my $TMPDIR = $ENV{'TMPDIR'};

#set environment variable for $NETMHC
my $NETMHC_env = $ENV{'NETMHC_env'};

#set environment variable for $NETMHCpan
my $NETMHCpan_env = $ENV{'NETMHCpan_env'};

#set environment variable for $PickPocket
my $PICKPOCKET_env = $ENV{'PICKPOCKET_env'};

# Define session Id
my $sId = $$;

################################################################################################################
# Define default values for optional variables and check if provided options and files are in the right format #
################################################################################################################

# Which mehtod to use
if (!defined $method) {
	$method = "NetMHCcons";
}
elsif (defined $method and $method eq "") {
	print "ERROR:insufficient arguments for option -method\nCheck ussage of the program using -h option\n";
	exit;
}
if (($method ne "NetMHCcons") and ($method ne "NetMHCpan") and ($method ne "NetMHC") and ($method ne "PickPocket")) {
	print "Wrong name of the method specified: choose from 'NetMHCcons', 'NetMHCpan', 'NetMHC', 'PickPocket'\n";
	exit;
}

# Working directory
if (!defined $rdir) {
	$rdir = $ENV{'NETMHCcons'};
}
elsif (defined $rdir and $rdir eq "") {
	print "ERROR: insufficient arguments for option -rdir\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $rdir) {
	print "ERROR: directory $rdir doesn't exist\n";
	exit;
}

# NetMHC method
if (!defined $NetMHC) {	
	$NetMHC = $NETMHC_env;
}
elsif (defined $NetMHC and $NetMHC eq "") {
	print "ERROR: insufficient arguments for option -mhc\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $NetMHC) {
	print "ERROR: program $NetMHC doesn't exist\n";
	exit;
}

# Specifying NetMHCpan method
if (!defined $NetMHCpan) {
	$NetMHCpan = $NETMHCpan_env;
}
elsif (defined $NetMHCpan and $NetMHCpan eq "") {
	print "ERROR: insufficient arguments for option -pan\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $NetMHCpan) {
	print "ERROR: program $NetMHCpan doesn't exist\n";
	exit;
}

# Specifying PickPocket method
if (!defined $PickPocket) {
	$PickPocket = $PICKPOCKET_env;
}
elsif (defined $PickPocket and $PickPocket eq "") {
	print "ERROR: insufficient arguments for option -pick\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $PickPocket) {
	print "ERROR: program $PickPocket doesn't exist\n";
	exit;
}

# Number of peptides
if(!defined $n_pep ) {
	$n_pep = 50;
}
unless (&isInt($n_pep)) {
	print "ERROR: number of peptides should be an integer\n";
	exit;
}

# Number of binders
if (!defined $n_bind) {
	$n_bind = 10;
}
unless (&isInt($n_bind)) {
	print "ERROR: number of binders should be an integer\n";
	exit;
}

# Affinity threshold for strong binders
if(!defined $thr_aff_S ) {
	$thr_aff_S = 50;
}
unless (&isAnumber($thr_aff_S)) {
	print "ERROR: threshold affinity for strong binders should be a number\n";
	exit;
}

# Affinity threshold for weak binders
if(!defined $thr_aff_W ) {
	$thr_aff_W = 500;
}
unless (&isAnumber($thr_aff_W)) {
	print "ERROR: threshold affinity for weak binders should be a number\n";
	exit;
}

# Rank threshold for weak binders
if (!defined $thr_rank_W) {
	$thr_rank_W = 2;
}
unless (&isAnumber($thr_rank_W)) {
	print "ERROR: threshold rank value for weak binders should be a number\n";
	exit;
}

# Rank threshold for strong binders
if (!defined $thr_rank_S) {
	$thr_rank_S = 0.5;
}
unless (&isAnumber($thr_rank_S)) {
	print "ERROR: threshold rank value for strong binders should be a number\n";
	exit;
}

# If filter for the output is defined
if (!defined $filter) {
	$filter = 0;
}
elsif (defined $filter and $filter eq "") {
	print "ERROR:insufficient arguments for option -filter\nCheck ussage of the program using -h option\n";
	exit;
}

# Rank threshold for filtering output
unless (&isAnumber($rank_f)) {
	print "ERROR: threshold rank value for filtering should be a number\n";
	exit;
}
# Affinity threshold for filtering output
unless (&isAnumber($aff_f)) {
	print "ERROR: threshold affinity value for filtering should be a number\n";
	exit;
}


# If filter is set to "Yes" or 1 (for the command line), give the default values for rank and affinity
my $filter_message;
if ($filter == 1) {
	if (!defined $rank_f) {
		$rank_f = 2;
	}
	if (!defined $aff_f) {
		$aff_f = 500;
	}
}
elsif ($filter == 0 and (defined $rank_f or defined $aff_f)) {
	if (!defined $w) {
		$filter_message = "# Filter was set to 0, input will not be filtered. Threshold value(s) will be ignored\n";
	}
}
	
# Distance to the nearest neighbour
if (!defined $dist) {
	$dist = 0.1;
}
unless (&isAnumber($dist)) {
	print "ERROR: distance should be a number\n";
	exit;
}

# Input file format
my $input_format = "";
if (!defined $inptype) {
	$inptype = 0;
	$input_format = "FASTA";
}
elsif ($inptype == 0) {
	$input_format = "FASTA";
}
elsif ($inptype == 1) {
	$input_format = "PEPTIDE";
}
unless (($inptype == 0 or $inptype == 1)) {
	print "ERROR: input type should be 1 for peptide format\n";
	exit;
}

# Peptide length
if (!defined $length) {
	$length = 9;
}
elsif (defined $length and $length eq "") {
	print "ERROR: insufficient arguments for option -length\nCheck ussage of the program using -h option\n";
	exit;
}
my @lengths = split (",", $length);
foreach my $length (@lengths) {
	unless (&isInt($length)) {
		print "ERROR: length of peptides should be an integer or a set of integers divided by commas\n";
		exit;
	}
	unless (($length >= 8 and $length <= 15)) {
		print "ERROR: peptide length should fall in a range [8-15]\n";
		exit;
	}
}

# File name for xls output
my $xlsfilename;
if (!defined $xlsfile) {
	$xlsfilename = "NetMHCcons_out.xls";
}
elsif (defined $xlsfile and $xlsfile eq "") {	
	print "ERROR: insufficient arguments for option -xlsfile\nCheck ussage of the program using -h option\n";
	exit;
}
else {
	$xlsfilename = $xlsfile;
}

# Temporary directory
if (!defined $tdir) {
	$tdir = "$TMPDIR/tmp_$sId"; 
}
elsif (defined $tdir and $tdir eq ""){
	print "ERROR: insufficient arguments for option -tdir\nCheck ussage of the program using -h option\n";
	exit;
}

# Training data file for NetMHC method
if (! defined $training_MHC) {
	$training_MHC = "$rdir/data/training.count";
}
elsif (defined $training_MHC and $training_MHC eq "") {
	print "ERROR: insufficient arguments for option -trMHC\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $training_MHC) {
	print "ERROR: file $training_MHC doesn't exist\n";
	exit;
}

# Training data file for NetMHCpan and PickPocket methods
if (! defined $training_Pan_Pick) {
	$training_Pan_Pick = "$rdir/data/training.count";
}
elsif (defined $training_Pan_Pick and $training_Pan_Pick eq "") {
	print "ERROR: insufficient arguments for option -trPanPick\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $training_Pan_Pick) {
	print "ERROR: file $training_Pan_Pick doesn't exist\n";
	exit;
}

# Saving the names of all known alleles
my %all_alleles = ();
open (IN, "<", "$rdir/data/alleles.dat") or die "Can not open the file with all possible allele names $!\n";
	while (defined (my $line = <IN>)) {
		chomp $line;
		$all_alleles{$line} = 1;
	}
close IN;
	
# Saving the alleles in question (one allele or a comma separated file with many alleles)
my @alleles = ();
if (! defined $alleles) {
	$alleles = ("HLA-A02:01");
}
elsif (defined $alleles and $alleles eq "") {
	print "ERROR: insufficient arguments for option -a\nCheck ussage of the program using -h option\n";
	exit;
}
@alleles = split (",", $alleles);
foreach my $query_allele (@alleles) {
	unless (&isAllele($query_allele)) {
		print "ERROR: allele name $query_allele is in a wrong format or doesn't exist in allele list\n";
		exit;
	}
}

# If a1 option is used
if (! defined $a1) {
	$a1 ="";
}
elsif (defined $a1 and $a1 eq "") {
	print "ERROR: insufficient arguments for option -a1\nCheck ussage of the program using -h option\n";
	exit;
}

# If a2 option is used
if (! defined $a2) {
	$a2 ="";
}
elsif (defined $a2 and $a2 eq "") {
	print "ERROR: insufficient arguments for option -a2\nCheck ussage of the program using -h option\n";
	exit;
}

# If a3 option is used
if (! defined $a3) {
	$a3 ="";
}
elsif (defined $a3 and $a3 eq "") {
	print "ERROR: insufficient arguments for option -a3\nCheck ussage of the program using -h option\n";
	exit;
}

# If a4 option is used
if (! defined $a4) {
	$a4 ="";
}
elsif (defined $a4 and $a4 eq "") {
	print "ERROR: insufficient arguments for option -a4\nCheck ussage of the program using -h option\n";
	exit;
}

# Input file	
if (! defined $file and ! defined $help) {
	print "ERROR: no input data\n";
	print "Usage: NetMHCcons [-h] [args] -f [fastafile/peptidefile]\n";
	exit;
}
elsif (defined $file and $file eq ""){
	print "ERROR: insufficient arguments for option -f\nCheck ussage of the program using -h option\n";
	exit;
}
elsif (defined $file) {
	unless (-e $file) {
		print "ERROR: input file $file doesn't exist\n";
		exit;
	}
}

# When the help option is defined
if (defined $help) {
	&usage();
	exit;
}

# Test if the training file has all required columns (Peptide, number of peptides available for the allele, number of binders available for the allele
open (IN, "<", $training_MHC) or die "Can not open the file with the training data $! \n";
while (defined (my $training_file_line = <IN>)) {
	chomp $training_file_line;
	my @test = split (" ", $training_file_line);
	unless (scalar @test == 3) {
		print "ERROR: file with the training data has a wrong formating\n";
		exit;
	}
}
close IN;

# Testing if the input file is in FASTA format when no type specified and in peptide format when peptide type is specified
my %expected_affinities = ();
my $flag_expected_aff = 0;
if ($inptype == 0) {
	open (IN, "<", $file) or die "Can not open input file $! \n";
	my $first = 0;
	while ($first == 0 and defined (my $input_line =<IN>)) {
		if ($input_line !~ m/^>/) {
			print "ERROR: Input file is not in FASTA format\n";
			exit;
		}
		$first = 1;
	}
	close IN;
}
elsif ($inptype == 1) {
	my @lengths_peptides = ();
	open (IN, "<", $file) or die "Can not open input file $! \n";
	while (defined (my $input_line =<IN>)) {
		chomp $input_line;
		if ($input_line !~/^\s/) {
			my @tmp = split (" ", $input_line);
			my $pep = $tmp[0];
			if (defined ($tmp[1]) and $tmp[1] ne "") {
				$expected_affinities{$pep} = $tmp[1];
				$flag_expected_aff = 1;
			}		
			push (@lengths_peptides, length($pep));
			$length = length($pep);
			if ($pep !~ m/^[A-Za-z]{8,15}$/) {
				print "ERROR: Wrong format of the input file - check if the file contains peptides [8-15]\n";
				exit;
			}
		}
	}
	close IN;
	for (my $i = 0; $i < $#lengths_peptides; $i++) {
		if ($lengths_peptides[$i] != $lengths_peptides[$i+1]) {
			print "ERROR: Wrong format of the input file - Peptide lenght must be equal for all peptides\n";
			exit;
		}
	}
}

# Creating temporary directory
mkdir $tdir;
unless (-e $tdir) {
	print "ERROR: directory $tdir doesn't exist\n";
	exit;
}

# If a1 or a2 or a3 or a4 allele options are specified
if ($method eq 'NetMHCcons' and $a1 ne "") {
	@alleles = split (",", $a1);
}
elsif ($method eq 'NetMHCpan' and $a2 ne "") {
	@alleles = split (",", $a2);
}
elsif ($method eq 'NetMHC' and $a3 ne "") {
	@alleles = split (",", $a3);
}
elsif ($method eq 'PickPocket' and $a4 ne "") {
	@alleles = split (",", $a4);
}
foreach my $query_allele (@alleles) {
	unless (&isAllele($query_allele)) {
		print "ERROR: allele name $query_allele is in a wrong format or doesn't exist in allele list\n";
		exit;
	}
}

# Define a variable to save hlaseq option to be parsed to NetMHCpan and PickPocket programs
my $hlaseq_option = "";
if (defined $hlaseq and $hlaseq ne "") {
	$hlaseq_option = "-hlaseq $hlaseq";
	@alleles = ("USER_DEF");
}
elsif (defined $hlaseq and $hlaseq eq "") {
	print "ERROR: insufficient arguments for option -hlaseq\nCheck ussage of the program using -h option\n";
	exit;
}

########################
# Saving required data #
########################

# Hash for number of data points available for each allele in the trianing set
my %peptides = ();
# Hash for numbers of binders available for each allele in the training set
my %binders = (); 

# Open the file with training information for NetMHC method and save data from it
open (IN, "<", $training_MHC) or die "Can not open the file with training data for reading $! \n";
while (defined (my $line = <IN>)) {
	chomp $line;
	my @data = split (" ", $line);
	my $allele = $data[0];
	$peptides{$allele} = $data[1];
	$binders{$allele} = $data[2];
}
close IN;

# Open the file with training pseudosequences and saving pseudosequences and allele names
my %pseudo2allele = ();
open (IN, "<", "$rdir/data/training.pseudo") or die "Can not open the file with training data for reading $! \n";
while (defined (my $line = <IN>)) {
	chomp $line;
	my @data = split (" ", $line);
	my $pseudo = $data[1];
	$pseudo2allele{$pseudo} = $data[0];	
}
close IN;

# Make a temporary file with the list of alleles from the training set to be used to find distances
system("cat $training_Pan_Pick > $tdir/training_alleles.dat");

#######################################
# Printing the results -initial lines #
#######################################

my $pos = sprintf("%6s", "pos");
my $allele_print = sprintf("%12s", "Allele");
my $peptide = sprintf("%15s", "peptide");
my $identity = sprintf("%16s", "Identity");
my $score = sprintf("%13s", "1-log50k(aff)");
my $affinity = sprintf ("%12s", "Affinity(nM)");
my $rank_print = sprintf("%6s", "\%Rank");
my $level = " BindingLevel";

my $new_thr_aff_S = sprintf("%.3f", $thr_aff_S);
my $new_thr_aff_W = sprintf("%.3f", $thr_aff_W);
my $new_aff_f;
if (defined $aff_f) {
	$new_aff_f = sprintf("%.3f", $aff_f);
}
my $expected_affinity = "";
if ($flag_expected_aff == 1) {
	$expected_affinity = sprintf("%11s", "Exp_Binding");
}

# Define hashes to be used for saving the output into the xls file
my %pos =();
my %pep = ();
my %log = ();
my %nm = ();
my %rank = ();
my %bl = ();
my %prot_id = ();

########################################################################################################################################
################################### HLA sequence option is used ########################################################################
########################################################################################################################################

# Starting calculations 
# Obtaining the results when the hla sequence is pasted or uploaded by the user (only one sequence at a time is possible)      
if ($hlaseq_option ne "") {	
	my $pseudosequence = "";
	
	system ("$rdir/bin/mhcfsa2psseq.$PLATFORM -p $rdir/data/all_varcontacts.nlist -r $rdir/data/B0702.fsa -m $rdir/data/BLOSUM50  $hlaseq | grep -v '#' > $tdir/tmp_pseudo_$sId");
	open (IN, "<", "$tdir/tmp_pseudo_$sId") or die "Can not open the file with created pseudosequence $tdir/tmp_pseudo_$sId $!\n";
	while (defined (my $line = <IN>)) {
		chomp $line;
		$pseudosequence = $line;
	}
	close IN;

	# check if pseudosequence contains gaps
	if ($pseudosequence =~ m/-/) {
		print "\nError. MHC input sequence incomplete. Pseudo sequence $pseudosequence contains gaps\n";
		exit;
	}
	#check if pseudosequence is composed of letters only
	if ($pseudosequence !~ m/^\w+$/) {
		print "\nError. The error occurred obtaining pseudosequence from MHC sequence\n";
		exit;
	}	
	# Getting the allele from the pseudosequence
	my $allele = "";
	if  (exists $pseudo2allele{$pseudosequence}) {
		$allele = $pseudo2allele{$pseudosequence};
	}
		
	# Calculate distance to the nearest neighbour
	my ($nn_dist, $neighbour);	
		
	system("echo $pseudosequence | $rdir/bin/pseudofind -m $rdir/data/BLOSUM50 -r $tdir/training_alleles.dat -p $rdir/data/MHC_pseudo.dat -s -S | grep -v '#' > $tdir/pseudodist_$sId.out");
	open(IN, "<", "$tdir/pseudodist_$sId.out") or die "Can not open the file with pseudodistances$!\n";
	while(defined(my $line =<IN>)) {
		chomp $line;
		if ($line =~ m/^\S+\s+(\S+)\s+(\d+\.\d+)/) {
			$neighbour = $1;
			$nn_dist = $2;
		}
	}
	close IN;	
		
	if ($filter == 1) {
		if (!defined $aff_f) {
			$aff_f= 500;
		}
	}	
       
	# Here the query allele is empty because "hlaseq" option is specified, so hlaseq option is given
	my $query_allele = "";
	my $RESULT="";
	# Define the variable to save which method (or combination of methods in Consensus) was finally used
	my $final_method = "";
	
	# Defining variables to find the number of strong and weak binders
	my $count_strong = 0;
	my $count_weak = 0;

####################################################	
### From here it depends which method was used #####
####################################################
	
	###########################
	# If NetMHC was specified #
	###########################
	if ($method eq "NetMHC") {		
			
		# Setting the query allele
		$query_allele = &change_allele_format($allele);
		
		# Getting the results from NetMHC method		
		my @MHC_Results = ();
		if (! defined $query_allele or ! exists $peptides{$query_allele} ){
			print "\nERROR: NetMHC can not produce output for the uploaded sequence\n";
			exit;
		}
		foreach my $length (@lengths){		
			my @MHC_Results_tmp = &runMHC ($NetMHC, $inptype, $length, $query_allele, $file);
			push (@MHC_Results, @MHC_Results_tmp);			
		}
		for(my $i = 0; $i <= $#MHC_Results; $i++) {			
			my @tmp_MHC = split (" ", $MHC_Results[$i]);
			my $m = 1;
			# Take results from the result lines
			if ($MHC_Results[$i] =~ m/\s+(WB|SB)\s+/) {
				$m = 0;
			}			
			# Changing the format of output
			my $pos = sprintf ("%6s", $tmp_MHC[0]);
			my $allele = sprintf("%12s", "USER_DEF");
			my $peptide = sprintf("%15s", $tmp_MHC[1]);
			my $identity = sprintf("%16s", $tmp_MHC[5-$m]);
			my $log_score_MHC = $tmp_MHC[2];
			my $log_score = sprintf("%13.3f", $log_score_MHC);
			my $affinity = sprintf("%12.2f", $tmp_MHC[3]);
			
			# Check if expected affinities for peptides were specified and prepare them for printing if they were
			my $expected = "";			
			if ($flag_expected_aff == 1) {
				$expected = $expected_affinities{$tmp_MHC[1]};
				$expected = sprintf("%12.3f", $expected);
			}		
		
			# Finding the level of binding
			my $level ="";					
			if ($affinity <= $thr_aff_S) {
				$level = "<=SB";				
			}
			elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
				$level = "<=WB";				
			}
			$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
		}
	}
	
	##############################
	# If NetMHCpan was specified #
	##############################
	if ($method eq "NetMHCpan") {
		# Getting results from NetMHCpan method
		my @Pan_Results = ();
		foreach my $length (@lengths) {
			my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele, $file);
			push (@Pan_Results, @Pan_Results_tmp);
		}
		foreach my $result_line (@Pan_Results) {
			if ($result_line =~ m/^\s+\d+\s+/) {
				my @tmp_Pan = split (' ', $result_line);
				
				# Changing the format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				my $log_score = sprintf("%13s", $tmp_Pan[4]);
				my $affinity = sprintf("%12s", $tmp_Pan[5]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}		
												
				# Finding the level of binding and counting strong and weak binders
				my $level ="";
				my $count_strong = 0;
				my $count_weak = 0;			
				if ($affinity <= $thr_aff_S) {
					$level = "<=SB";					
				}
				elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
					$level = "<=WB";					
				}
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
			}
		}
	}
	##############################
	#If PickPocket was specified #
	##############################
	if ($method eq "PickPocket") {
		# Getting results from PickPocket method
		my @Pick_Results = ();
		foreach my $length (@lengths) {
			my @Pick_Results_tmp = &runPick ($PickPocket, $inptype, $length, $hlaseq_option, $query_allele, $file);
			push (@Pick_Results, @Pick_Results_tmp);
		}
		
		for(my $i = 0; $i <= $#Pick_Results; $i++) {
			my @tmp_Pick = ();
			# Take results from the result lines
			@tmp_Pick = split (' ', $Pick_Results[$i]);
			
			# Change the format of the output	
			my $pos = sprintf ("%6s", $tmp_Pick[0]);
			my $allele = sprintf("%12s", $tmp_Pick[1]);
			my $peptide = sprintf("%15s", $tmp_Pick[2]);
			my $identity = sprintf("%16s", $tmp_Pick[3]);
			my $log_score_Pick = $tmp_Pick[4];
			my $log_score = sprintf("%13.3f", $log_score_Pick);
			
			# Check if expected affinities for peptides were specified and prepare them for printing if they were
			my $expected = "";			
			if ($flag_expected_aff == 1) {
				$expected = $expected_affinities{$tmp_Pick[2]};
				$expected = sprintf("%12.3f", $expected);
			}			
			
			# Getting IC50 affinity values in nM
			my $affinity_Pick = exp( (1-$log_score_Pick)*log(50000));
			my $affinity = sprintf("%12.2f", $affinity_Pick);
			
			# Finding the level of binding
			my $level ="";				
			if ($affinity <= $thr_aff_S) {
				$level = "<=SB";				
			}
			elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
				$level = "<=WB";				
			}
			$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
		}
	}	
	
	#############################
	#If NetMHCcons was specified
	############################
	if ($method eq "NetMHCcons") {
			
		# If the distance to the nearest neighbour is equal to 0, it means the allele is part of the training set and we use NetMHC+NetMHCpan
		if ($nn_dist eq "0.000") {
			$final_method = "NetMHC+NetMHCpan";
			
			# Getting results for NetMHCpan and NetMHC
			my @Pan_Results = ();
			my @MHC_Results = ();
		
			foreach my $length (@lengths){
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele, $file);
				push (@Pan_Results, @Pan_Results_tmp);
			}			
	
			foreach my $length (@lengths){		
				# Change allele format to the old version
				my $query_allele_mhc = &change_allele_format($allele);
				my @MHC_Results_tmp = &runMHC ($NetMHC, $inptype, $length, $query_allele_mhc, $file);
				push (@MHC_Results, @MHC_Results_tmp);
			}
			
			# Checking if the two methods produced output for the same peptides before combining them
			my $message = &CompareOutputs (\@Pan_Results, \@MHC_Results, "NetMHCpan", "NetMHC");
			if ($message ne "") {
				print $message;
				exit;
			}

			for(my $i = 0; $i <= $#Pan_Results; $i++) {
				# Define arrays for saving the results for different methods
				my @tmp_Pan = ();
				my @tmp_MHC = ();
				# Take results from the result lines
				@tmp_Pan = split (' ', $Pan_Results[$i]);
				@tmp_MHC = split (' ', $MHC_Results[$i]);
			
				# Changing format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}
				
				# Getting log score of the consensus method - average of NetMHC and NetMHCpan methods
				my $log_score_MHC = $tmp_MHC[2];
				my $log_score_Pan = $tmp_Pan[4];
				my $log_score_Cons = ($log_score_MHC + $log_score_Pan)/2;
				my $log_score = sprintf("%13.3f", $log_score_Cons);				
			
				# Getting IC50 affinity values in nM
				my $affinity_Cons = exp( (1-$log_score_Cons)*log(50000));
				my $affinity = sprintf("%12.2f", $affinity_Cons);			
				
				# Finding the level of binding
				my $level ="";				
				if ($affinity <= $thr_aff_S) {
					$level = "<=SB";					
				}
				elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
					$level = "<=WB";					
				}
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
			}
		}		
	
		# If the distance to the nearest neighbour is less than "dist" (can be specified as an option on the command line), use NetMHCpan
		elsif ($nn_dist < $dist and $nn_dist > 0) {
			$final_method = "NetMHCpan";
			
			# Getting results for NetMHCpan method
			my @Pan_Results = ();
			foreach my $length (@lengths) {
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele, $file);
				push (@Pan_Results, @Pan_Results_tmp);
			}
			foreach my $result_line (@Pan_Results) {			
				if ($result_line =~ m/^\s+\d+\s+/) {					
					my @tmp_Pan = split (' ', $result_line);
					
					# Changing the format of the output					
					my $pos = sprintf ("%6s", $tmp_Pan[0]);
					my $allele = sprintf("%12s", $tmp_Pan[1]);
					my $peptide = sprintf("%15s", $tmp_Pan[2]);
					my $identity = sprintf("%16s", $tmp_Pan[3]);
					my $log_score = sprintf("%13s", $tmp_Pan[4]);
					my $affinity = sprintf("%12s", $tmp_Pan[5]);
					
					# Check if expected affinities for peptides were specified and prepare them for printing if they were
					my $expected = "";			
					if ($flag_expected_aff == 1) {
						$expected = $expected_affinities{$tmp_Pan[2]};
						$expected = sprintf("%12.3f", $expected);
					}
									
					# Finding the level of binding
					my $level ="";				
					if ($affinity <= $thr_aff_S) {
						$level = "<=SB";						
					}
					elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
						$level = "<=WB";						
					}
					$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
				}
			}
		}
		
		# If the distance to the nearest neighbour is more than "dist" (can be specified as an option on the command line) - use NetMHCpan + PickPocket
		else {
			$final_method = "NetMHCpan+PickPocket";
			
			# Getting results for NetMHCpan and PickPocket methods
			my @Pan_Results = ();
			my @Pick_Results = ();
			foreach my $length (@lengths) {
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele, $file);
				my @Pick_Results_tmp = &runPick ($PickPocket, $inptype, $length, $hlaseq_option, $query_allele, $file);
				push (@Pan_Results, @Pan_Results_tmp);
				push (@Pick_Results, @Pick_Results_tmp);
			}
		
			# Checking if the two methods produced output for the same peptides before combining them
			my $message = &CompareOutputs (\@Pan_Results, \@Pick_Results, "NetMHCpan", "PickPocket");
			if ($message ne "") {
				print $message;
				exit;
			}

			for(my $i = 0; $i <= $#Pan_Results; $i++) {
				# Define arrays for saving the results for different methods
				my @tmp_Pan = ();
				my @tmp_Pick = ();
				# Take results from the result lines
				@tmp_Pan = split (' ', $Pan_Results[$i]);
				@tmp_Pick = split (' ', $Pick_Results[$i]);
			
				# Changing format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}				
			
				# Getting log score of the consensus method - average of PickPocket and NetMHCpan methods
				my $log_score_Pick = $tmp_Pick[4];
				my $log_score_Pan = $tmp_Pan[4];
				my $log_score_Cons = ($log_score_Pick + $log_score_Pan)/2;
				my $log_score = sprintf("%13.3f", $log_score_Cons);
			
				# Getting IC50 affinity values in nM
				my $affinity_Cons = exp( (1-$log_score_Cons)*log(50000));
				my $affinity = sprintf("%12.2f", $affinity_Cons);			
				
				# Finding the level of binding
				my $level ="";			
				if ($affinity <= $thr_aff_S) {
					$level = "<=SB";					
				}
				elsif ($affinity <= $thr_aff_W and $affinity > $thr_aff_S) {
					$level = "<=WB";					
				}			
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
			}	
		}
	}
	
	##### Modifying the results if the filtering or the sorting was used ########
	my $RESULT_MOD = "";
	my @result_lines = split ("\n", $RESULT);
	foreach my $result_line (@result_lines) {
		my @scores = split (" ", $result_line);
		my $pos = sprintf ("%6s", $scores[0]);
		my $allele = sprintf("%12s", $scores[1]);
		my $peptide = sprintf("%15s", $scores[2]);
		my $identity = sprintf("%16s", $scores[3]);
		my $log_score = sprintf("%13.3f", $scores[4]);
		my $affinity = sprintf("%12.2f", $scores[5]);
		my $level = "";
		my $expected = "";
		
		if (defined $scores[7]) {
			$expected = sprintf("%12.3f", $scores[6]);
			$level = $scores[7];
		}
		elsif (defined $scores[6] and !defined $scores[7] and $flag_expected_aff != 1) {
			$level = $scores[6];
		}
		elsif (defined $scores[6] and !defined $scores[7] and $flag_expected_aff == 1) {
			$expected = sprintf("%12.3f", $scores[6]);
		}
			
		# Finding the number of strong and weak binders
		if ($level eq "<=SB") {
			$count_strong++;			
		}
		elsif ($level eq "<=WB") {
			$count_weak++;			
		}	
	
		#If the filter for filtering output was defined
		if ($filter == 1) {
			if ($affinity <= $aff_f) {
				$RESULT_MOD .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
			}
			if (($level eq "<=SB") and ($affinity > $aff_f)) {
				$count_strong--;
			}
			if (($level eq "<=WB") and ($affinity > $aff_f)) {
				$count_weak--;
			}
		}
		else {
			$RESULT_MOD .= "$pos $allele $peptide $identity $log_score $affinity $expected $level\n";
		}
	}
				
	# Printing result lines into the file
	open (OUT, ">", "$tdir/results.out") or die "Can not open the file for writing $!\n";
	print OUT $RESULT_MOD;
	close OUT;

	# If the sort option was specified, sort the results based on affinity
	my $FINAL_RESULT = "";
	if (defined $sort)  {
		system("cat $tdir/results.out | sort -nrk5 > $tdir/final_results.out");
		open (IN, "<", "$tdir/final_results.out") or die "Can not open the file final_results.out for reading $!\n";
		while (defined (my $line = <IN>)) {
			$FINAL_RESULT .= $line;
		}
		close IN;
	}
	else {
		$FINAL_RESULT = $RESULT_MOD;
	}	
	
	print "# Method: $method\n\n" , 
      	"# Input is in $input_format format\n\n" ,
      	"# Peptide length $length\n\n" ,
      	"# Threshold for Strong binding peptides (IC50)\t$new_thr_aff_S nM\n",
      	"# Threshold for Weak binding peptides (IC50)\t$new_thr_aff_W nM\n";
      	if ($filter == 1) {
		print "\n# Threshold for filtering output (IC50)\t\t$new_aff_f nM\n";
	}
      	print "\n# Thresholds for \%Rank are ignored\n\n";
      	if (defined $filter_message) {
		print "\n$filter_message\n";
	}
	
	### Printing the results for MHC sequence
	print "\n# MHC sequence from FASTA file\n";
	print "\n# Distance to the nearest neighbour ( $neighbour ) in the training set: $nn_dist\n";
	
	if ($method eq "NetMHCcons") {
		print "\n# NetMHCcons = $final_method\n\n";
	}	
	print "-------------------------------------------------------------------------------\n",
      	      "$pos $allele_print $peptide $identity $score $affinity $expected_affinity $level\n",
              "-------------------------------------------------------------------------------\n";	
	print $FINAL_RESULT;
	print "-------------------------------------------------------------------------------\n",
	      "Number of strong binders: $count_strong Number of weak binders: $count_weak\n",
	      "-------------------------------------------------------------------------------\n";
	
	# Saving output into hashes of arrays to be used for the .xls output
	if (defined $xls) {
		
		my @result_lines = split ("\n", $RESULT);
	
		foreach my $result_line (@result_lines) {
			chomp $result_line;
			my @tmp = split (" ", $result_line);
			push (@{$pos{"USER_DEF"}}, $tmp[0]);
			push (@{$pep{"USER_DEF"}}, $tmp[2]);
			push (@{$prot_id{"USER_DEF"}}, $tmp[3]);
	  		push (@{$log{"USER_DEF"}}, $tmp[4]);
			push (@{$nm{"USER_DEF"}}, $tmp[5]);
		
			if (defined $tmp[6] and $flag_expected_aff !=1) {
				push (@{$bl{"USER_DEF"}}, 1);
			}
			elsif (defined $tmp[7]) {
				 push (@{$bl{"USER_DEF"}}, 1);
			}
			else {
				push (@{$bl{"USER_DEF"}}, 0);
			}
		}
	}
	
} ## if ($hlaseq_option ne "")

########################################################################################################################################
################################### Finish HLA sequence option #########################################################################
########################################################################################################################################

# Obtaining the results for the alleles chosen from the list (multiple alleles separated by commas are possible)
else {     

my $print_count = 0;
	
# Go through each allele specified by the user
foreach my $query_allele (@alleles) {
	# Calculate distance to the nearest neighbour
	my ($nn_dist, $neighbour);	
		
	system("echo $query_allele | $rdir/bin/pseudofind -m $rdir/data/BLOSUM50 -r $tdir/training_alleles.dat -p $rdir/data/MHC_pseudo.dat -s -S | grep -v '#' > $tdir/pseudodist_$sId.out");
	open(IN, "<", "$tdir/pseudodist_$sId.out") or die "Can not open the file with pseudodistances$!\n";
	while(defined(my $line =<IN>)) {
		chomp $line;
		if ($line =~ m/^\S+\s+(\S+)\s+(\d+\.\d+)/) {
			$neighbour = $1;
			$nn_dist = $2;
		}
	}
	close IN;	

	# Define variable for saving result lines
	my $RESULT ="";	
	# Define the variable to save which method (or combination of methods in Consensus) was finally used
	my $final_method = "";	
	# Defining variables to find the number of strong and weak binders
	my $count_strong = 0;
	my $count_weak = 0;

	
###########################
# If NetMHC was specified #
###########################

	if ($method eq "NetMHC") {
		# change the format of alleles for NEtMHC method		
		my $query_allele_n2o = &change_allele_format($query_allele);

		my @MHC_Results = ();
		if (! exists $peptides{$query_allele_n2o}){
			print "\nERROR: NetMHC can not produce output for allele $query_allele\n";
			exit;
		}
		foreach my $length (@lengths){		
			my @MHC_Results_tmp = &runMHC ($NetMHC, $inptype, $length, $query_allele_n2o, $file);
			push (@MHC_Results, @MHC_Results_tmp);
			
		}
		for(my $i = 0; $i <= $#MHC_Results; $i++) {			
			my @tmp_MHC = split (" ", $MHC_Results[$i]);
			my $m = 1;
			# Take results from the result lines
			if ($MHC_Results[$i] =~ m/\s+(WB|SB)\s+/) {
				$m = 0;
			}
			
			# Change format of the output
			my $pos = sprintf ("%6s", $tmp_MHC[0]);
			# my $allele = sprintf("%12s", $tmp_MHC[6-$m]);
			my $allele = sprintf("%12s", $query_allele);
			my $peptide = sprintf("%15s", $tmp_MHC[1]);
			my $identity = sprintf("%16s", $tmp_MHC[5-$m]);
			my $log_score_MHC = $tmp_MHC[2];
			my $log_score = sprintf("%13.3f", $log_score_MHC);
			my $affinity = sprintf("%12.2f", $tmp_MHC[3]);
			
			# Check if expected affinities for peptides were specified and prepare them for printing if they were
			my $expected = "";			
			if ($flag_expected_aff == 1) {
				$expected = $expected_affinities{$tmp_MHC[1]};
				$expected = sprintf("%12.3f", $expected);
			}
			
			# Finding the rank
			my $rank = "";
			my @RANKS = ();
			my @SCORES = ();
			open (IN, "<", "$rdir/data/thresholds/NetMHC/$query_allele_n2o.thr") or die "Can not open the file $rdir/data/thresholds/NetMHC/$query_allele_n2o.thr $!\n";
			while (defined (my $line =<IN>)){
				chomp $line;
				my @tmp = split (" ", $line);
				push (@RANKS, $tmp[1]);
				push (@SCORES, $tmp[3]);
			}
			close IN;
			my $flag = 0;
			for (my $i = 0; $i <= $#RANKS; $i++) {
				if ($log_score >= $SCORES[$i] and $flag == 0) {
					$flag = 1;
					$rank = $RANKS[$i];
				}
				if ($i == $#RANKS and $log_score < $SCORES[$i]) {
					$rank = $RANKS[$#RANKS];
				}
			}
			$rank = sprintf("%6.2f", $rank);
			
			## Finding the level of binding
			my $level ="";			
			if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
				$level = "<=SB";
			}
			elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
				$level = "<=WB";
			}
			$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
		}
	}
	
##############################
# If NetMHCpan was specified #
##############################

	if ($method eq "NetMHCpan") {
		# Getting results for NetMHCpan

		# Translate query names if nessesary
		my $query_allele_o2n = &change_allele_format_2($query_allele);

		my @Pan_Results = ();
		foreach my $length (@lengths) {
			my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
			push (@Pan_Results, @Pan_Results_tmp);
		}
		foreach my $result_line (@Pan_Results) {
			if ($result_line =~ m/^\s+\d+\s+/) {
				my @tmp_Pan = split (' ', $result_line);
				
				# Changing format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				#my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $allele = sprintf("%12s", $query_allele );
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				my $log_score = sprintf("%13s", $tmp_Pan[4]);
				my $affinity = sprintf("%12s", $tmp_Pan[5]);
				my $rank = sprintf("%8s", $tmp_Pan[6]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}				
			
				# Finding the level of binding and counting strong and weak binders
				my $level ="";
				if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
					$level = "<=SB";
				}
				elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
					$level = "<=WB";				
				}
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";			
			}
			
		}
	}
	
###############################
# If PickPocket was specified #
###############################

	if ($method eq "PickPocket") {	

		# Translate query names if nessesary
        	my $query_allele_o2n = &change_allele_format_2($query_allele);
	
		# Getting results for PickPocket
		my @Pick_Results = ();
		foreach my $length (@lengths) {
			my @Pick_Results_tmp = &runPick ($PickPocket, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
			push (@Pick_Results, @Pick_Results_tmp);
		}
		
		for(my $i = 0; $i <= $#Pick_Results; $i++) {
			my @tmp_Pick = ();
			# Take results from the result lines
			@tmp_Pick = split (' ', $Pick_Results[$i]);
			
			# Changing format of the output	
			my $pos = sprintf ("%6s", $tmp_Pick[0]);
			#my $allele = sprintf("%12s", $tmp_Pick[1]);
			my $allele = sprintf("%12s", $query_allele);
			my $peptide = sprintf("%15s", $tmp_Pick[2]);
			my $identity = sprintf("%16s", $tmp_Pick[3]);
			my $log_score_Pick = $tmp_Pick[4];
			my $log_score = sprintf("%13.3f", $log_score_Pick);
			
			# Check if expected affinities for peptides were specified and prepare them for printing if they were
			my $expected = "";			
			if ($flag_expected_aff == 1) {
				$expected = $expected_affinities{$tmp_Pick[2]};
				$expected = sprintf("%12.3f", $expected);
			}
				
			# Getting IC50 affinity values in nM
			my $affinity_Pick = exp( (1-$log_score_Pick)*log(50000));
			my $affinity = sprintf("%12.2f", $affinity_Pick);
						
			# Finding the rank
			my $rank = "";
			my @RANKS = ();
			my @SCORES = ();
			open (IN, "<", "$rdir/data/thresholds/PickPocket/$query_allele_o2n.thr") or die "Can not open the file $rdir/data/thresholds/PickPocket/$query_allele_o2n.thr $! \n";
			while (defined (my $line =<IN>)) {
				chomp $line;
				my @tmp = split (" ", $line);
				push (@RANKS, $tmp[1]);
				push (@SCORES, $tmp[3]);
			}
			close IN;
			my $flag = 0;
			for (my $i = 0; $i <= $#RANKS; $i++) {
				if ($log_score >= $SCORES[$i] and $flag == 0) {
					$flag = 1;
					$rank = $RANKS[$i];
				}
				if ($i == $#RANKS and $log_score < $SCORES[$i]) {
					$rank = $RANKS[$#RANKS];
				}
			}
			$rank = sprintf("%6.2f", $rank);
			
			# Finding the level of binding
			my $level ="";			
			if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
				$level = "<=SB";				
			}
			elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
				$level = "<=WB";				
			}
			$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
		}		
	}	
	
#####################################################################################
#If NetMHCcons was specified - from here everything is related to NetMHCcons method #
#####################################################################################

	if ($method eq "NetMHCcons") {	
	
#####################
# Condition 1: Allele is part of the training set
#####################
	
	# If allele is part of the training data and has less than "np" (can be specified as an option on the command line) data points and less than "nb" (can be specified as an option on the command line) binders - run NetMHCpan

	# Translation allele name from old to new
	my $query_allele_o2n = change_allele_format_2($query_allele);
	# Translation allele name from new to old 
	my $query_allele_n2o = change_allele_format($query_allele);

	if (exists $peptides{$query_allele_n2o}) {		
		my $no_peps = $peptides{$query_allele_n2o}; 
		my $no_bind = $binders{$query_allele_n2o};	
		if ($no_peps < $n_pep and $no_bind < $n_bind) {
			$final_method = "NetMHCpan";
			my @Pan_Results =();
			foreach my $length (@lengths) {
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
				push (@Pan_Results, @Pan_Results_tmp);
			}
			foreach my $result_line (@Pan_Results) {
				if ($result_line =~ m/^\s+\d+\s+/) {
					my @tmp_Pan = split (' ', $result_line);
					
					# Changing format of the output
					my $pos = sprintf ("%6s", $tmp_Pan[0]);
					#my $allele = sprintf("%12s", $tmp_Pan[1]);
					my $allele = sprintf("%12s", $query_allele);
					my $peptide = sprintf("%15s", $tmp_Pan[2]);
					my $identity = sprintf("%16s", $tmp_Pan[3]);
					my $log_score = sprintf("%13s", $tmp_Pan[4]);
					my $affinity = sprintf("%12s", $tmp_Pan[5]);
					my $rank = sprintf("%8s", $tmp_Pan[6]);
					
					# Check if expected affinities for peptides were specified and prepare them for printing if they were
					my $expected = "";			
					if ($flag_expected_aff == 1) {
						$expected = $expected_affinities{$tmp_Pan[2]};
						$expected = sprintf("%12.3f", $expected);
					}					
			
					# Finding the level of binding
					my $level ="";			
					if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
						$level = "<=SB";						
					}
					elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
						$level = "<=WB";						
					}
					$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
				}			
			
			}
		}
		# If allele is part of the training set and has more than "np" (can be specified as an option on the command line) data points and more than "nb" (can be specified as an option on the command line) binders - run NetMHCpan + NetMHC
		else {
			$final_method = "NetMHC+NetMHCpan";			
			# Getting results for NetMHC and NetMHCpan
			my @Pan_Results = ();
			my @MHC_Results = ();
		
			foreach my $length (@lengths){
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
				push (@Pan_Results, @Pan_Results_tmp);
			}		

			foreach my $length (@lengths){		
				#change allele format to the old version
				#my $query_allele_mhc = &change_allele_format($query_allele);
				#my @MHC_Results_tmp = &runMHC ($NetMHC, $inptype, $length, $query_allele_mhc, $file);
				my @MHC_Results_tmp = &runMHC ($NetMHC, $inptype, $length, $query_allele_n2o, $file);

				push (@MHC_Results, @MHC_Results_tmp);
			}
			
			# Checking if the two methods produced output for the same peptides before combining them
			my $message = &CompareOutputs (\@Pan_Results, \@MHC_Results, "NetMHCpan", "NetMHC");
			if ($message ne "") {
				print $message;
				exit;
			}
			my $flag = 0;
			for(my $i = 0; $i <= $#Pan_Results; $i++) {
				#define arrays for saving the results fro different methods
				my @tmp_Pan = ();
				my @tmp_MHC = ();
				# Take results from the result lines
				@tmp_Pan = split (' ', $Pan_Results[$i]);
				@tmp_MHC = split (' ', $MHC_Results[$i]);
			
				# Changing format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				#my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $allele = sprintf("%12s", $query_allele);
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}				
			
				# Getting log score of the consensus method - average of NetMHC and NetMHCpan methods
				my $log_score_MHC = $tmp_MHC[2];
				my $log_score_Pan = $tmp_Pan[4];
				my $log_score_Cons = ($log_score_MHC + $log_score_Pan)/2;
				my $log_score = sprintf("%13.3f", $log_score_Cons);
			
				# Getting IC50 affinity values in nM
				my $affinity_Cons = exp( (1-$log_score_Cons)*log(50000));
				my $affinity = sprintf("%12.2f", $affinity_Cons);			
				
				# Finding the rank
				my $rank = "";
				my @RANKS = ();
				my @SCORES = ();
				my $query_allele1;
			
				#open (IN, "<", "$rdir/data/thresholds/MHC_Pan/$query_allele.thr") or die "Can not open the file $rdir/data/thresholds/MHC_Pan/$query_allele.thr $! \n";
				open (IN, "<", "$rdir/data/thresholds/MHC_Pan/$query_allele_n2o.thr") or die "Can not open the file $rdir/data/thresholds/MHC_Pan/$query_allele_n2o.thr $! \n";
				while (defined (my $line =<IN>)) {
					chomp $line;
					my @tmp = split (" ", $line);
					push (@RANKS, $tmp[1]);
					push (@SCORES, $tmp[3]);
				}
				close IN;
				my $flag = 0;
				for (my $i = 0; $i <= $#RANKS; $i++) {
					if ($log_score >= $SCORES[$i] and $flag == 0) {
						$flag = 1;
						$rank = $RANKS[$i];
					}
					if ($i == $#RANKS and $log_score < $SCORES[$i]) {
						$rank = $RANKS[$#RANKS];
					}
				}
				$rank = sprintf("%6.2f", $rank);
				
				# Finding the level of binding
				my $level ="";				
				if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
					$level = "<=SB";					
				}
				elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
					$level = "<=WB";					
				}
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
			}
		}
	}

#####################
# Condition 2: Allele is not part of the training set
#####################

	# If allele is not part of the training set, calculate distance to the nearest neighbour
	else {	
		# If the distance to the nearest neighbour is less than "dist" (can be specified as an option on the command line), use NetMHCpan
		if ($nn_dist < $dist) {
			$final_method = "NetMHCpan";
			
			# Getting results for NetMHCpan method
			my @Pan_Results = ();
			foreach my $length (@lengths) {
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
				push (@Pan_Results, @Pan_Results_tmp);
			}
			foreach my $result_line (@Pan_Results) {
				if ($result_line =~ m/^\s+\d+\s+/) {
					# Obtaining result lines
					my @tmp_Pan = split (' ', $result_line);
					
					# Changing format of the output					
					my $pos = sprintf ("%6s", $tmp_Pan[0]);
					#my $allele = sprintf("%12s", $tmp_Pan[1]);
					my $allele = sprintf("%12s", $query_allele);
					my $peptide = sprintf("%15s", $tmp_Pan[2]);
					my $identity = sprintf("%16s", $tmp_Pan[3]);
					my $log_score = sprintf("%13s", $tmp_Pan[4]);
					my $affinity = sprintf("%12s", $tmp_Pan[5]);
					my $rank = sprintf("%8s", $tmp_Pan[6]);
					
					# Check if expected affinities for peptides were specified and prepare them for printing if they were
					my $expected = "";			
					if ($flag_expected_aff == 1) {
						$expected = $expected_affinities{$tmp_Pan[2]};
						$expected = sprintf("%12.3f", $expected);
					}				
					
					# Finding the level of binding
					my $level ="";				
					if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
						$level = "<=SB";					
					}
					elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
						$level = "<=WB";						
					}
					$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
				}			
			}
		}
		# If the distance to the nearest neighbour is more tha "dist" (can be specified as an option on the command line) - use NetMHCpan + PickPocket
		else {
			$final_method = "NetMHCpan+PickPocket";
			
			# Getting results for NetMHCpan and PickPocket method
			my @Pan_Results = ();
			my @Pick_Results = ();
			foreach my $length (@lengths) {
				my @Pan_Results_tmp = &runPan ($NetMHCpan, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
				my @Pick_Results_tmp = &runPick ($PickPocket, $inptype, $length, $hlaseq_option, $query_allele_o2n, $file);
				push (@Pan_Results, @Pan_Results_tmp);
				push (@Pick_Results, @Pick_Results_tmp);
			}
		
			# Checking if the two methods produced output for the same peptides before combining them
			my $message = &CompareOutputs (\@Pan_Results, \@Pick_Results, "NetMHCpan", "PickPocket");
			if ($message ne "") {
				print $message;
				exit;
			}

			for(my $i = 0; $i <= $#Pan_Results; $i++) {
				# Define arrays for saving the results fro different methods
				my @tmp_Pan = ();
				my @tmp_Pick = ();
				# Take results from the result lines
				@tmp_Pan = split (' ', $Pan_Results[$i]);
				@tmp_Pick = split (' ', $Pick_Results[$i]);
			
				# Change format of the output
				my $pos = sprintf ("%6s", $tmp_Pan[0]);
				#my $allele = sprintf("%12s", $tmp_Pan[1]);
				my $allele = sprintf("%12s", $query_allele);
				my $peptide = sprintf("%15s", $tmp_Pan[2]);
				my $identity = sprintf("%16s", $tmp_Pan[3]);
				
				# Check if expected affinities for peptides were specified and prepare them for printing if they were
				my $expected = "";			
				if ($flag_expected_aff == 1) {
					$expected = $expected_affinities{$tmp_Pan[2]};
					$expected = sprintf("%12.3f", $expected);
				}				
			
				# Getting log score of the consensus method - average of PickPocket and NetMHCpan methods
				my $log_score_Pick = $tmp_Pick[4];
				my $log_score_Pan = $tmp_Pan[4];
				my $log_score_Cons = ($log_score_Pick + $log_score_Pan)/2;
				my $log_score = sprintf("%13.3f", $log_score_Cons);
			
				# Getting IC50 affinity values in nM
				my $affinity_Cons = exp( (1-$log_score_Cons)*log(50000));
				my $affinity = sprintf("%12.2f", $affinity_Cons);			
				
				# Finding the rank
				my $rank = "";
				my @RANKS = ();
				my @SCORES = ();
				open (IN, "<", "$rdir/data/thresholds/Pan_Pick/$query_allele_o2n.thr") or die "Can not open the file $rdir/data/thresholds/Pan_Pick/$query_allele_o2n.thr $! \n";
				while (defined (my $line =<IN>)) {
					chomp $line;
					my @tmp = split (" ", $line);
					push (@RANKS, $tmp[1]);
					push (@SCORES, $tmp[3]);
				}
				close IN;
				my $flag = 0;
				for (my $i = 0; $i <= $#RANKS; $i++) {
					if ($log_score >= $SCORES[$i] and $flag == 0) {
						$flag = 1;
						$rank = $RANKS[$i];
					}
					if ($i == $#RANKS and $log_score < $SCORES[$i]) {
						$rank = $RANKS[$#RANKS];
					}
				}
				$rank = sprintf("%6.2f", $rank);
				
				# Finding the level of binding
				my $level ="";			
				if (($affinity <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
					$level = "<=SB";					
				}
				elsif (($affinity <= $thr_aff_W and $affinity > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
					$level = "<=WB";					
				}
				$RESULT .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
			}	
		}
	}
	} # if ($method eq "NetMHCcons")
		
	##### Modifying the results if the filtering or the sorting was used ########
	my $RESULT_MOD = "";
	my @result_lines = split ("\n", $RESULT);
	foreach my $result_line (@result_lines) {		
		my @scores = split (" ", $result_line);
		my $pos = sprintf ("%6s", $scores[0]);
		my $allele = sprintf("%12s", $scores[1]);
		my $peptide = sprintf("%15s", $scores[2]);
		my $identity = sprintf("%16s", $scores[3]);
		my $log_score = sprintf("%13.3f", $scores[4]);
		my $affinity = sprintf("%12.2f", $scores[5]);
		my $rank = sprintf("%6.2f", $scores[6]);
		my $level = "";
		my $expected = "";
		
		if (defined $scores[8]) {
			$expected = sprintf("%12.3f", $scores[7]);
			$level = $scores[8];
		}
		elsif (defined $scores[7] and !defined $scores[8] and $flag_expected_aff != 1) {
			$level = $scores[7];
		}
		elsif (defined $scores[7] and !defined $scores[8] and $flag_expected_aff == 1) {
			$expected = sprintf("%12.3f", $scores[7]);
		}
		
		# Finding the number of strong and weak binders
		if ($level eq "<=SB") {
			$count_strong++;		
		}
		elsif ($level eq "<=WB") {
			$count_weak++;			
		}
		
		# If the filter for filtering output was defined
		if ($filter == 1) {
			if ($rank <= $rank_f or $affinity <= $aff_f) {
				$RESULT_MOD .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
			}
			if (($level eq "<=SB") and ($rank > $rank_f and $affinity > $aff_f)) {
				$count_strong--;
			}
			if (($level eq "<=WB") and ($rank > $rank_f and $affinity > $aff_f)) {
				$count_weak--;
			}
		}
		else {
			$RESULT_MOD .= "$pos $allele $peptide $identity $log_score $affinity $rank $expected $level\n";
		}
	}

	# Printing result lines into the file
	open (OUT, ">", "$tdir/results.out") or die "Can not open the file for writing $!\n";
	print OUT $RESULT_MOD;
	close OUT;

	# If the sort option was specified, sort the results based on affinity
	my $FINAL_RESULT = "";
	if (defined $sort)  {
		system("cat $tdir/results.out | sort -nrk5 > $tdir/final_results.out");
		open (IN, "<", "$tdir/final_results.out") or die "Can not open the file final_results.out for reading $!\n";
		while (defined (my $line = <IN>)) {
			$FINAL_RESULT .= $line;
		}
		close IN;
	}
	else {
		$FINAL_RESULT = $RESULT_MOD;
	}	
	
	if ($print_count == 0) {
		print "# Method: $method\n\n" , 
      		"# Input is in $input_format format\n\n" ,
      		"# Peptide length $length\n\n" ,
      		"# Threshold for Strong binding peptides (IC50)\t$new_thr_aff_S nM\n",
      		"# Threshold for Weak binding peptides (IC50)\t$new_thr_aff_W nM\n\n",
      		"# Threshold for Strong binding peptides (\%Rank)\t$thr_rank_S%\n",
      		"# Threshold for Weak binding peptides (\%Rank)\t$thr_rank_W%\n";
		if ($filter == 1) {
			print "\n# Threshold for filtering output (\%Rank)\t$rank_f%\n",
			"# Threshold for filtering output (IC50)\t\t$new_aff_f nM\n";
		}
		if (defined $filter_message) {
			print "\n\n$filter_message\n";
		}
		$print_count++;		
	}
	
	### Printing the results for the allele
	print "\n# Allele: $query_allele\n";
	print "\n# Distance to the nearest neighbour ( $neighbour ) in the training set: $nn_dist\n";
	
	if ($method eq "NetMHCcons") {
		print "\n# NetMHCcons = $final_method\n\n";
	}
	print "---------------------------------------------------------------------------------------\n",
      	      "$pos $allele_print $peptide $identity $score $affinity $rank_print $expected_affinity $level\n",
              "---------------------------------------------------------------------------------------\n";	
	print $FINAL_RESULT;
	print "---------------------------------------------------------------------------------------\n",
	      "Number of strong binders: $count_strong Number of weak binders: $count_weak\n",
	      "---------------------------------------------------------------------------------------\n";
	
	# Using the result lines to save output into hashes of arrays in order to use them for .xls output later
	if (defined $xls) {
		
		my @lines = split ("\n", $RESULT);
		
		foreach my $result_line (@lines) {
			chomp $result_line;
			my @tmp = split (" ", $result_line);
		
			push (@{$pos{$query_allele}}, $tmp[0]);
			push (@{$pep{$query_allele}}, $tmp[2]);
			push (@{$prot_id{$query_allele}}, $tmp[3]);
	  		push (@{$log{$query_allele}}, $tmp[4]);
			push (@{$nm{$query_allele}}, $tmp[5]);
			push (@{$rank{$query_allele}}, $tmp[6]);
			if (defined $tmp[7] and $flag_expected_aff !=1) {
				push (@{$bl{$query_allele}}, 1);
			}
			elsif (defined $tmp[8]) {
				 push (@{$bl{$query_allele}}, 1);
			}
			else {
				push (@{$bl{$query_allele}}, 0);
			}
			
		}
		
	}

}  #closing foreach $query_allele

} # else from (if $hlaseq_option = "TRUE")

# Acessing and printing  the results into the tab separated file for .xls output if the option "xls" was defined
if (defined $xls) {
	foreach my $allele1 (@alleles) {
		foreach my $allele2 (@alleles) {
			if ($#{$pos{$allele1}} != $#{$pos{$allele2}}) {
				print "ERROR occured when creating Excel file: $allele1 and $allele2 resulted into different number of peptides\n";
				exit;
			}
		}
	}		
	my $first = "FALSE";
	my $file_name;
	if (defined $w) {
		$file_name = "/usr/opt/www/pub/CBS/services/NetMHCcons-1.1/tmp/$$"."_NetMHCcons.xls";
	}
	else {
		$file_name = $xlsfilename;
	}
	open (OUT, ">", $file_name) or die "Can not open the file$!\n";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			print OUT "\t\t\t\t$query_allele";
			$first = "TRUE";
		}
		else {
			print OUT "\t\t\t$query_allele";
		}
	}
	$first = "FALSE";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			if (!exists $rank{$alleles[0]}) {
				print OUT "\nPos\tPeptide\tID\t1-log50k\tnM";
			}
			else {
				print OUT "\nPos\tPeptide\tID\t1-log50k\tnM\tRank";
			}			
			$first = "TRUE";
		}
		else {
			print OUT "\t1-log50k\tnM\tRank";
		}
	}
	print OUT "\tAve\tNB\n";

	for (my $i = 0; $i <= $#{$pos{$alleles[0]}}; $i++) {
		my $nb = $bl{$alleles[0]}[$i];
		my $log_sum = $log{$alleles[0]}[$i];
		if (!exists $rank{$alleles[0]}) {
			print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$nm{$alleles[0]}[$i]\t";
		}
		else {
			print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$nm{$alleles[0]}[$i]\t$rank{$alleles[0]}[$i]\t";
		}
		for (my $n = 1; $n <= $#alleles; $n++){
			$nb += $bl{$alleles[$n]}[$i];
			$log_sum += $log{$alleles[$n]}[$i];
			print OUT "$log{$alleles[$n]}[$i]\t$nm{$alleles[$n]}[$i]\t$rank{$alleles[$n]}[$i]\t";
		}
		my $avg = $log_sum/scalar(@alleles);
		$avg = sprintf("%.4f", $avg);
		print OUT "$avg\t$nb\n";
	}
	close OUT;
	my $short_name = "/services/NetMHCcons-1.1/tmp/$$"."_NetMHCcons.xls";
	if (defined $w) {
		print "Link to output xls file <a href='$short_name'>NetMHCcons_out.xls</a>\n";
	}
	
}
#################### Finished creating .xls file #####################

# Deleting temporary directory if the dirty mode was not chosen
if (!defined $dirty) {
	system ("rm -r $tdir");
}

#####################
# Subroutines
#####################

# Function to run NetMHCpan method. Parameters: $Pan - full path to the NetMHCpan program, $allele - allele for which we need predictions, $file - file with peptides
sub runPan {
	my ($Pan, $inptype, $length, $hlaseq, $allele, $file) = @_;
	my $tmpfile = "$tdir/NetMHCpan_$sId.out";
	my @output =();
	#system ("$Pan -inptype $inptype -l $length $hlaseq -a $allele -f $file  | grep -v '#' > $tmpfile");
	#system ("$Pan -inptype $inptype -l $length $hlaseq -a $allele -f $file -ic50 | grep -v '#' > $tmpfile");
	#sleep(10);
	
	open(IN, '-|', "$Pan -inptype $inptype -l $length $hlaseq -a $allele -f $file -ic50 | grep -v '#'" );
	
	#open (IN, "<", $tmpfile) or die "Can't open the file with NetMHCpan output $tmpfile $!\n";
	while (defined (my $line = <IN>)) {
		chomp $line;
		if ($line =~ m/^\s*\d+\s+\w+/) {
			push (@output, $line);
			
		}
		
	}
	
	close IN;
	unlink $tmpfile;
	return @output;
}
	
# Function to run PickPocket method. Parameters: $Pick - full path to the PickPocket program, $allele - allele for which we need predictions, $file - file with peptides
sub runPick {
	my ($Pick, $inptype, $length, $hlaseq, $allele, $file) = @_;
	my $tmpfile = "$tdir/PickPocket_$sId.out";
	my @output;
	#system ("$Pick -inptype $inptype -l $length $hlaseq -a $allele -f $file | grep -v '#' > $tmpfile");
	#sleep(10);
	#open (IN, "<", $tmpfile) or die "Can't open the file with PickPocket output $tmpfile $!\n";
	open(IN, '-|', "$Pick -inptype $inptype -l $length $hlaseq -a $allele -f $file | grep -v '#'");
	while (defined (my $line = <IN>)) {
		chomp $line;
		if ($line =~ m/^\s*\d+\s+\w+/) {
			push (@output, $line);
			
		}
	}
	close IN;
	unlink $tmpfile;
	return @output;
}		

# Function to run NetMHC method. Parameters: $MHC - full path to the NetMHC program, $allele - allele for which we need predictions, $file - file with peptides
sub runMHC {
	my ($MHC, $inptype, $length, $allele, $file) = @_;
	my $tmpfile = "$tdir/NetMHC_$sId.out";
	my @output;
	
	if ($inptype == 0) {
		$inptype = "";
	}
	else {
		$inptype = "-p";
	}
	#system ("$MHC $inptype -a $allele -l $length $file > $tmpfile");
	#sleep(10);
	#open (IN, "<", $tmpfile) or die "Can't open the file with NetMHC output $tmpfile $!\n";
	
	
	open(IN, '-|', "$MHC $inptype -a $allele -l $length $file");
	while (defined (my $line = <IN>)) {
		chomp $line;
		if ($line =~ m/^\s*\d+\s+\w+\s+\d+\.*\d*/) {
			push (@output, $line);
			
		}
	}
	close IN;
	#if (scalar @output == 0) {
		#die "\nnetMHC-3.3 produced an error\n";
	#}
	unlink $tmpfile;
	return @output;
}		
		
# program usage	
sub usage {
	print "\nUsage: ./NetMHCcons [-h] [args] [fastafile/peptidefile]\n";
	print "Command line options:\n\n";
	printf ("%-16s\t%-30s\t%-12s\n",  "PARAMETER", "DEFAULT VALUE", "DESCRIPTION");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-rdir dirname]", "$NETMHCcons", "Home directory for NetMHCcons");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-tdir dirname]", "$TMPDIR/tmp_\$\$", "Temporary directory");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-mhc filename]", "$NetMHC", "Full path to NetMHC program");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-pan filename]", "$NetMHCpan", "Full path to NetMHCpan program");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-pick filename]", "$PickPocket", "Full path to PickPocket program");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-trMHC filename]", "$rdir/data/training.count", "File with the training data for NetMHC method");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-trPanPick filename]", "$rdir/data/training.count", "File with the training data for NetMHCpan and PickPocket methods");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-method name]", "NetMHCcons", "Specify method [NetMHCcons,NetMHCpan,NetMHC,PickPocket]");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-np int]", "50", "Threshold for the number of peptides in the training set at which different methods are chosen");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-nb int]", "10", "Threshold for the number of binders in the training set at which different methods are chosen");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-affS float]", "50.000", "Threshold for strong binders (IC50)");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-affW float]", "500.000", "Threshold for weak binders (IC50)");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-rankS float]", "0.5", "Threshold for strong binders (\%Rank)");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-rankW float]", "2", "Threshold for weak binders (\%Rank)");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-filter int]", "0", "Filter output [1]");	
	printf ("%-16s\t%-30s\t%-12s\n",  "[-affF float]", "500", "Threshold for filtering output (IC50). Used only when filter is 1");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-rankF float]", "2", "Threshold for filtering output (\%Rank). Used only when filter is 1");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-dist float]", "0.1", "Threshold for the distance to the nearest neighbour at which different methods are chosen");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-inptype int]", "0", "Input type [0] FASTA [1] Peptide");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-xls]", "0", "Save output into xls file");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-xlsfile filename]", "NetMHCcons_out.xls", "File name for xls output");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-length int_array]", "9", "Peptide length from a range [8-15]. Several lengths can be chosen, giving them separated by commas, ex.:8,9,10");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-hlaseq filename]", " ", "File with full length HLA sequence");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-s]", "0", "Sort output on descending affinity");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-a name]", "HLA-A02:01", "HLA allele");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-f filename]", " ", "File with the input data");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-dirty]", "0", "Dirty mode, leave tmp dir+files");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-w]", "0", "w option for webface");
	printf ("%-16s\t%-30s\t%-12s\n",  "[-h]", "0", "Print this message and exit");
}	
			
## Check format of a number
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

## Check for integer
sub isInt {
    my $test=shift;
    
    if ($test =~ m/^\d+$/ && $test>=0) {
	return 1; }
    else {
	return 0; }
}		

# Check allele format

sub isAllele {
	my $test=shift;
	
	if (exists $all_alleles{$test}) {
		return 1;
	}
	else {
		return 0;
	}
}		
		
sub change_allele_format {
	my $allele = shift;
	#changed allele to return
	my $changed_allele; 
	
	# Getting converted names for the old allele names
	my %names;

	open (IN, "<", "$rdir/data/allele_translation") or die "Can not open the file $rdir/data/allele_translation $!\n";
	while (defined (my $line = <IN>)) {
		chomp $line;
		my @tmp = split (" ", $line);
		$names{$tmp[0]} = $tmp[1];
	}
	close IN;

	if (defined $allele and exists $names{$allele}) {
		$changed_allele = $names{$allele};
	}
	else {
		$changed_allele = $allele if defined $allele;
	}
	return $changed_allele;
}

sub change_allele_format_2 {
        my $allele = shift;
        #changed allele to return
        my $changed_allele;

        # getting converted names for the Mamu alleles
        my %names;

        open (IN, "<", "$rdir/data/allele_translation") or die "Can not open the file $rdir/data/allele_translation $!\n";
        while (defined (my $line = <IN>)) {
                chomp $line;
                my @tmp = split (" ", $line);
                $names{$tmp[1]} = $tmp[0];
        }
        close IN;

        if (defined $allele and exists $names{$allele}) {
                $changed_allele = $names{$allele};
        }
        else {
                $changed_allele = $allele if defined $allele;
        }
        return $changed_allele;
}


sub CompareOutputs {
	my ($ref1, $ref2, $method1, $method2) = @_;
	my $i1;
	my $i2;
	if ($method1 eq "NetMHC") {
		$i1 = 1;
	}
	elsif($method1 eq "NetMHCpan" or $method1 eq "PickPocket") {
		$i1 = 2;
	}
	if ($method2 eq "NetMHC") {
		$i2 = 1;
	}
	elsif($method2 eq "NetMHCpan" or $method2 eq "PickPocket") {
		$i2 = 2;
	}
	my @peplist1 = ();
	my @peplist2 = ();
	
	foreach my $line (@{$ref1}) {
		my @tmp = split (" ", $line);
		push @peplist1, $tmp[$i1];
	}
	foreach my $line (@{$ref2}) {
		my @tmp = split (" ", $line);
		push @peplist2, $tmp[$i2];
	}

	my $err_message = "";
	if (scalar(@peplist1) != scalar(@peplist2)) {
		$err_message = "The two methods $method1 and $method2 produced different outputs, number of peptides not the same\n";
	}
	else {
		for (my $a = 0; $a <= $#peplist1; $a++){
			if ($peplist1[$a] ne $peplist2[$a]){
				$err_message = "The two methods $method1 and $method2 produced different outputs, some peptides are not the same\n";
			}
		}
	}

	return $err_message;
}		
				
		
