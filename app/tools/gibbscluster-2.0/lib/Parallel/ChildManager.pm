=head1 NAME

Parallel::ChildManager - A parallel processing child manager

=head1 SYNOPSIS


Example 1: Job mode

  use Parallel::ChildManager;

  my $cm = new ChildManager($max_process, $time_limit);

  foreach my $data (@all_data) {
    # Starts a parallel job (child)
    $cm->start("/foo/bar/job $data"); }
  $cm->wait_all_children;
  # or
  $cm->start("job1", "job2", "job3");
  $cm->wait_all_children;
  # or
  $cm->start(@all_jobs_in_an_array);
  $cm->wait_all_children;


Example 2: Empty call, forking mode

  use Parallel::ChildManager;

  my $cm = new ChildManager($max_process);

  foreach my $data (@all_data) {
    # Forks and returns the pid for the child:
    my $pid = $cm->start; 
    if ($pid) {
       # This is the parent
       # Do some parent code
       }
    else {
       # This is the child
       # Do child code with $data
       exit; # Must end with exit or you start executing main parent code
       }
  }
  $cm->wait_all_children;

Example 3: Children returns results directly back to the parent

  use Parallel::ChildManager;

  my $cm = new ChildManager($max_process);
  $cm->set_start_type('IPC');

  foreach my $data (@all_data) {
    # Forks and returns the pid for the child:
    my $pid = $cm->start; 
    if ($pid) {
       # This is the parent
       # Do some parent code
       }
    else {
       # This is the child
       # Do child code with $data
       $cm->return_result(\@results_array); # ends child execution
       }
  }
  $cm->wait_all_children;
  @all_results_references = $cm->get_all;
  
Example 4: Children returns results directly back to the parent in sequential order

  use Parallel::ChildManager;

  my $cm = new ChildManager($max_process, $timelimit);
  $cm->set_start_type('IPC');

  open(IN, '<', $file) or die "$!";
  my @pool;
  while (defined (my $line = <IN>)) {
     push(@pool, $line);
     next if scalar @pool < $poolsize;
     my $pid = $cm->start;
    if ($pid) {
       # This is the parent
       @pool = ();
       while (ref (my $res = $cm->get_next('NOBLOCK'))) {
		&store_results($res); }
    } else {
       # This is the child
       # Do child code with @pool data
       $cm->return_result(\@results_array); # ends child execution
    }
  }
  close IN;
  if (@pool and 0 == $cm->start) {
     # Do child code with last @pool data
     $cm->return_result(\@results_array); # ends child execution
   }
   while (ref (my $res = $cm->get_next('BLOCK'))) {
      &store_results($res); }

=head1 DESCRIPTION

This module is intended for running processes/jobs in parallel in an easy-to-use manner. It is developed with a twofold purpose in mind.
Primary use is "heavy-weight" job processing (can be any program) on a parallel machine in a batch-like manner. An example could be training many different neural networks in parallel and selecting the best of them. Basically any CPU intensive jobs.
Secondary use is for more "light-weight" operations (usually written in Perl) which can benefit from being done in parallel. An example could be a downloader, which retrieves many files in parallel. Basically any operations where the bottleneck is NOT the CPU or disk access (at least this goes for a single CPU machine).

In both uses, pains are taken to ensure that no process/job will hang if caught in an endless loop, as the user (if she is wise) can set a time limit for any job. Likewise it is possible (and advisable) to specify the maximum concurrent processes. The package will automatically schedule jobs, so maximum processes are running at any time.

Killing the parent will kill any children, grandchildren, etc. Killing a child will also kill that childs children.

The package was developed on SGI Origin 2000 with Perl 5.6, but works with 5.0. It makes use of the SGI automatic process scheduler to run jobs on free CPU's, that is; there is no code in the package which is SGI only. It all runs on a basic unix. But there is one parameter ($process_limit) which can be set to certain values and ensure a certain behaviour which only is valid under SGI/IRIX. This is checked for.

NOTE: You cannot have more than one instance of Parallel::ChildManager. That instance can be reused in child processes, though, in any which way you want it. 
NOTE: The package makes use of SIGTERM, so you should not use that signal. There is also some use of SIGINT and SIGQUIT, but these can be overridden.


=head1 METHODS

=over 4

=item new [$process_limit [, $time_limit [, $recursion_depth]]]

Instantiate a new Parallel::ChildManager object. There can only be one instance, but that can be used by the children too (it will be cleaned of children, but parameters are inherited).
$process_limit specifies the maximum number of children to fork off. The default is -1, which means unlimited number of concurrent children. Note: On multiprocessor SGI/IRIX this number can be negative (from -2 and down) and this means that at most that many processors will be used, still leaving CPU time to other jobs.
$time_limit specifies in minutes (can be fractional), how long the child is allowed to run, before its is terminated (killed). This is intended for limiting CPU usage of faulty children. The default is -1 and means no time limit.
$recursion_depth set to 0 as default and simply means that children can not have children of their own via the 'start' method (empty call, see 'start'). This is a simple means to ensure that process spawning does not run amok.

Example: $cm = new ChildManager(5,10);
Example: $cm = new ChildManager;


=item set_max_process $process_limit

This method sets the process limit. Can be called anytime. If lowering the limit eventual children are allowed to finished their buisness.

Example: $cm->set_max_process(5);


=item set_max_time $time_limit

This method sets the time limit. Can be called anytime. Will take effect at next death of a child.

Example: $cm->set_max_time(10);


=item set_child_priority

You can set the priority on the childs with this method. The number (priority) has to be meaningful to your system. Basically this invokes a Perl 'setpriority' call.

Example: $cm->set_child_priority(10);


=item set_recursion_depth $recursion_depth

This method sets the recursion depth. Can be called anytime. Will take effect at next spawning of a child. Makes only sense when using empty call of 'start' method using the default CHILD start type.

Example: $cm->set_recursion_depth(1);


=item set_start_type $start_type

This method influences heavily the way simple 'light-weight' forking behaves. The default value is set to CHILD which means that empty start calls just forks and returns the pid of the child to the parent. But if set to READ or WRITE, empty start calls returns a filehandle to the parent (still 0 to the child). The parent can either read from or write to the child through the filehandle with the childs STDOUT or STDIN connected to the other end of the pipe. If set to IPC the child will actually transfer back results to the parent through the return_result method. The results will be available to the parent through the get_next, get_ready or get_all methods.

Example: $cm->set_start_type('READ');


=item set_throttle $level

Sometimes the jobs you run are competing for a resource. This method will automatically cut down in a linear fashion on the number of max process if the running time of the longest running is greater than 10% (level 1), 20% (level 2) (and so forth), of the maximum allowed time limit. Default is 0 and means no throttle.

Example: $cm->set_throttle(4);


=item set_warning $level

The default is 0 which means no warnings, but can be set to 1 which is critical warnings (likely data loss) or 2 which also includes infomation on erronous but handled situations.

Example: $cm->set_warning(1);


=item start [$job] [@jobs]

This method spawns the child. There are really four different uses of the 'start' method. The first  is an empty call of the method, see example 2 in the synopsis and below;

Example: $pid = $cm->start;

This spawns a child, which is the (almost) exact copy of the parent, the only difference is the PID of the child is returned to the parent by the 'start' method, and 0 is returned to the child. That can be used to distinguish between parent and child and thus determine which code should be executed (example 2). If the child wants to have children, then you should set the recursion depth to 1, 2 for grandchildren etc... See also set_start_type.

Example: $cm->$set_start_type('READ'); $filehandle = $cm->start;

This later addition to the package spawns a child that the parent can read from or write to via the returned file handle. Undef is returned to the child, which uses STDOUT or STDIN to read/write from/to. Remember to close the file handle. Perhaps you might also consider unbuffered I/O.

Example: $cm->$set_start_type('IPC'); $pid = $cm->start;

IPC makes writeback of variables (results) possible from child to parent. This is a rather stong feature. See synopsis example 3.

Example: $cm->start("/foo/bar/job options");

This is the fourth form of the method call, see example 1 in the synopsis. This spawns a child which does a perl 'system' call with the parameter given. If more parameters are given (or an array), more children are spawned (according to the limitations set in the 'new' method call. All jobs are handled automagically with regard to timeout and manual killing. The method returns the childs PID (or an array of PIDs). The child never returns from this form of the 'start' method.
If child spawning does not succeed because of no free entries in the process table, 'start' waits a bit and tries again. If it fails for any other reason, it dies.


=item wait_next_child BLOCK | NOBLOCK

You can call this method to wait for the next child to terminate. This is can be a blocking wait (BLOCK), a nonblocking wait (NOBLOCK).  Returns the Process ID or -1 if no children or 0 for no children to be reaped in case of nonblocking waits. The method handles timeout on children (ie. you will not wait forever, if you have set a time limit and the child is in a never ending loop).

Example: $pid = $cm->wait_next_child('BLOCK');


=item wait_all_children

You can call this method to wait for all currently running children and children in queue. This is a blocking wait and usually used when using the package to run a lot of jobs in parallel.

Example: $cm->wait_all_children;


=item return_result

This method transfers the child's results back to the parent. Can only be used when the start type is set to IPC. Calling the method ends child execution. The parent utilises one (and only one) of the methods get_all, get_ready or get_next to get the data. The results has to packed as a reference to a single data structure like an array or a hash.

Example: $cm->return_result(\@result);


=item get_all

This method returns an array with references to your childrens results (from using return_result) in order of execution. If there is no data for a child due to some error (segmentation fault for instance) then undef is returned for that particular position in the array. This method must only be called after all children has finished, f.ex. after wait_all_children. You might want to consider your amount of data. get_all does not mix with get_ready or get_next.

Example: my @res = $cm->get_all;


=item get_ready

Returns a reference to the results of the next child to finish. The order is somewhat random. This is a blocking method unless you specify the parameter NOBLOCK. The data are removed after usage, so use this method if you are worried about memory. Returns undef if no results due to child crash. If using NOBLOCK it will return the number of active children if no result is currently ready. This means that it returns a reference if results are available or a number, where 0 will indicate no more results. get_ready does not mix with get_all or get_next.

Example: my $reference = $cm->get_ready;


=item get_next

Returns a reference to the results of the next child in order of execution. This is a blocking method unless you specify the parameter NOBLOCK. The data are removed after usage, so use this method if you are worried about memory. Returns undef if no results due to child crash, and 0 if no more results. If using NOBLOCK it will return the number of active children if the next result is currently unavailable. get_next does not mix with get_ready or get_all.

Example: my $reference = $cm->get_next(NOBLOCK);


=item get_status $level

Returns a line (perhaps multiline) describing the status of the run. Level 0 is one line status. Level 1 adds a list of individual problems. Level 2 adds a helping analysis.

Example: print STDERR $cm->get_status(0);


=item kill_late_children

Mainly used for internal administration of children. Calling this method terminates children, which have reached the time limit. The method does not wait for the children.

Example: $cm->kill_late_children;


=item kill_all

Mainly used for internal administration of children. Using this method terminates the current caller and all of his children.

Example: $cm->kill_all;


=item examine_load

Used for internal administration. This method figures out a conservative number (integer) of currently idle processors.

Example: $cm->examine_load;


=item reap

Used for internal administration of children when using READ and WRITE in set_start_type. Reaps all terminated children, ALSO those which disappeared due to instabilities.

Example: $cm->reap;


=head1 COPYRIGHT

Copyright (c) 2002-2009 Peter Wad Sackett

All right reserved. This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

All thanks to dLux (Szabó, Balázs) <dlux@kapu.hu> and his Parallel:ForkManager package, from which I stole initial ideas and methods shamelessly.


=head1 BUGS

It seems that children disappear, i.e. are harvested by some other mechanism when using the READ/WRITE start calls. There is currently a workaround in place.


=head1 AUTHOR

Peter Wad Sackett <pws@cbs.dtu.dk>

=cut

package ChildManager;
use strict;
use POSIX ":sys_wait_h";
use Socket;
use Data::Dumper;

$Parallel::ChildManager::VERSION = '1.4';

my $packagehash;


sub BEGIN {
   setpgrp(0,0); # New process group, useful for kills.
   $SIG{'INT'} = \&kill_all;
   $SIG{'QUIT'} = \&kill_all;
   $SIG{'TERM'} = \&kill_all;
}


sub new {
   my $self = shift @_;
   my ($maxprocess, $maxtime, $maxdepth);
   if (defined $packagehash) {
      warn "Only one ChildManager object is allowed.\n";
      return undef; }
   if (@_) {
      $maxprocess = int(shift @_);
      $maxprocess = -1 if $maxprocess == 0 ||
                          ($maxprocess < 0 && $ENV{'HOSTTYPE'} ne 'iris4d'); }
   else {
      $maxprocess = -1; }
   if (@_) {
      $maxtime = int(60 * shift @_);
      $maxtime = -1 if $maxtime < 1; }
   else {
      $maxtime = -1; }
   if (@_) {
      $maxdepth = int(shift @_);
      $maxdepth = 0 if $maxdepth < 1; }
   else {
      $maxdepth = 0; }
   $packagehash = {
      'max_proc' => $maxprocess,
      'max_time' => $maxtime,
      'rec_depth'=> $maxdepth,
      'orig_max_proc' => $maxprocess,
      'throttle' => 0,
      'warning' => 0,
      'workercount' => 0,
      'start_type' => 'CHILD',
      'robin_count' => 0,
      'ready_cpus' => 0,
      'last_check' => 0,
      'cur_proc' => 0,
      'dead_proc' => [],
      'status' => {},
      'proc_tab' => {}, };
   return bless($packagehash, ref($self) || $self);;
}


sub set_max_process {
   my ($self, $maxprocess) = @_;
   if ($maxprocess > 1 || $maxprocess == -1 ||
      ($maxprocess < -1 && $ENV{'HOSTTYPE'} eq 'iris4d')) {
      $self->{'max_proc'} = $maxprocess;
      $self->{'orig_max_proc'} = $maxprocess; }
   else {
      warn "The maximum process limit have to be a\n";
      warn "positive integer or -1 (for no limit).\n"; }
}


sub set_recursion_depth {
   my ($self, $recdepth) = @_;
   unless ($recdepth =~ m/^\d+$/) {
      warn "The maximum recursion depth have to be a non-negative\n";
      warn "integer (0 is no recursion, ie. children having children).\n"; }
   else {
      $self->{'rec_depth'} = $recdepth; }
}


sub set_child_priority {
   my ($self, $priority) = @_;
   if ($priority =~ m/^-?\d+$/) {
      $self->{'child_priority'} = $priority; }
   else {
      warn "Priority has to be a number meaningful to the system.\n"; }
}


sub set_max_time {
   my ($self, $maxtime) = @_;
   if (60 * $maxtime < 1 and $maxtime != -1) {
      warn "The maximum time limit is in minutes and have to be\n";
      warn "positive (fractionals allowed) or -1 (for no limit).\n"; }
   else {
      $self->{'max_time'} = int(60 * $maxtime); }
}


sub set_throttle {
   my ($self, $level) = @_;
   if ($level != m/^\d$/) {
      $self->{'throttle'} = $level; }
   else {
      warn "The throttle level is a number between 0 and 9, both inclusive\n"; }
}


sub set_warning {
   my ($self, $level) = @_;
   if ($level != m/^\^[012]$/) {
      $self->{'warning'} = $level; }
   else {
      warn "The warning is a number 0 for no warnings, 1 for critical warnings, 2 includes handled errors\n"; }
}


sub set_start_type {
   my ($self, $type) = @_;
   $type = uc($type);
   unless (' CHILD READ WRITE IPC ' =~ m/ $type /) {
      warn "Only parameters CHILD, READ, WRITE or IPC are permitted\n"; }
   else {
      $self->{'start_type'} = $type; }
   if ($type eq 'IPC') {
      socket($self->{'socket'}, PF_INET, SOCK_STREAM, getprotobyname('tcp')) or
         die("Socket call failed, $!");
      setsockopt($self->{'socket'}, SOL_SOCKET, SO_REUSEADDR, pack("l", 1)) or
         die("Setting socket options failed, $!");
      my $tries = 0;
      while ($tries < 10) {
         $self->{'socketport'} = int(2000+rand(3000));
         bind($self->{'socket'}, sockaddr_in($self->{'socketport'}, INADDR_LOOPBACK)) and last;
	 $tries++; }
      $tries == 10 and die("Binding socket failed, $!");
      listen($self->{'socket'},SOMAXCONN) or
         die("Listen socket call failed, $!");
      $self->{data} = {};
      $self->{order} = [];
#      $self->{proc_order} = {};  # debug
      $Data::Dumper::Indent = 0;
      $Data::Dumper::Purity = 1;
   }
}


sub start {
   my $self = shift @_;
   my ($pid, $perlfork, $dofork, @new_procs, @pid_list);
   # Get parameters
    if (@_) {
       $perlfork = 0;
       while (@_) {
          push(@new_procs, shift @_); } }
    else {
       $perlfork = 1; 
       if ($self->{'start_type'} eq 'READ' or $self->{'start_type'} eq 'WRITE') {
          $self->{'robin_count'}++;
	  if ($self->{'robin_count'} >= $self->{'max_proc'}>>2) {
	     $self->{'robin_count'} = 0;
	     $self->reap; } } }
   if ($perlfork && $self->{'rec_depth'} < 0) {
      warn "A child/process is trying to fork/create more children.\n";
      warn "Use 'set_recursion_depth' method to allow this.\n";
      return -1; }
   # IPC special - harvesting ready children is more important than starting new jobs
#   if ($self->{'start_type'} eq 'IPC') {
#      1 while $self->wait_next_child('NOBLOCK'); }
   # Ready to fork ??
   while (1) {
      if ($self->{'cur_proc'} < $self->{'max_proc'}) {
         $dofork = 1; }
      elsif ($self->{'cur_proc'} == 0) {
         $dofork = 1; }
      elsif ($self->{'max_proc'} == -1) {
         $self->wait_next_child('NOBLOCK');
         $dofork = 1; }
      elsif ($self->{'cur_proc'} < -$self->{'max_proc'}) {
         $self->examine_load;
         if ($self->{'ready_cpus'} > 0) {
            $dofork = 1;
            $self->{'ready_cpus'}--; }
         else {
            $dofork = -1; } }
      else {
         $dofork = 0; }
      if ($dofork == 1) {
         while (1) { # fork loop
	    my $filehandle;
	    unless ($perlfork) {
               $pid = fork; }
	    elsif ($self->{'start_type'} eq 'READ') {
	       $pid = open($filehandle, "-|"); }
	    elsif ($self->{'start_type'} eq 'WRITE') {
	       $pid = open($filehandle, "|-"); }
	    else {
	       $pid = fork;  }
            if ($pid) {  # parent
	       if (exists $self->{status}->{$pid}) {
	          $self->{status}->{$pid} .= 'P';
	          warn "Reuse of process id, child discarded, no loss incurred\n" if $self->{warning} == 2;
		  waitpid($pid, 0);
		  next; }
	       $self->{status}->{$pid} = ++$self->{workercount} . ' ';
               $self->{'cur_proc'}++;
               $self->{'proc_tab'}{$pid} = time;
	       unless ($perlfork) {
                  last; }
	       elsif ($self->{'start_type'} eq 'READ' or $self->{'start_type'} eq 'WRITE') {
	          return $filehandle; }
	       else {
	          push(@{$self->{'order'}}, $pid) if $self->{'start_type'} eq 'IPC';
#		  $self->{'proc_order'}->{$pid} = ++$self->{'workercount'} if $self->{'start_type'} eq 'IPC'; # debug
#		  warn "Forking worker: $self->{'workercount'}\tpid->$pid\t\tNext out: $self->{'proc_order'}->{$self->{'order'}->[0]}\n" if $self->{'start_type'} eq 'IPC';  # debug
	          return $pid;  } }
            if (defined $pid) {  # child
	       exit if exists $self->{status}->{$$};
               setpgrp(0,0); # New process group, useful for kills.
               setpriority(0, 0, $self->{'child_priority'})
                  if exists $self->{'child_priority'};
               if ($perlfork) {
	          $self->{workercount}++;
                  $self->{'rec_depth'}--;
                  $self->{'cur_proc'} = 0;
                  undef $self->{'proc_tab'};
                  $self->{'proc_tab'} = {};
		  if ($self->{'start_type'} eq 'IPC') {
		     close $self->{'socket'};
		     $self->{data} = {};
     		     $self->{order} = []; }
                  return 0 if $self->{'start_type'} eq 'CHILD' or $self->{'start_type'} eq 'IPC';
		  return undef; }
               system(shift @new_procs);
               exit; }
            if ($! =~ m/No more process/) { # recoverable error
               sleep 10;
               next; }
            die "Can't fork: $!\n"; }
         shift @new_procs;
         push(@pid_list, $pid);
         unless (@new_procs) {
            return @pid_list if $#pid_list > 0;
            return shift @pid_list; } }
      elsif ($dofork == 0) {  # wait for next child to die
         $self->wait_next_child('BLOCK'); }
      else {
          unless ($self->wait_next_child('NOBLOCK')) {
	     sleep 60;
	     $self->wait_next_child('NOBLOCK'); }
      }
   }
}


sub wait_next_child {
   my ($self, $blocktype) = @_;
   my ($pid, $timeleft);
   die "Undefined blocking type parameter: $blocktype"
      unless defined $blocktype and ' BLOCK NOBLOCK ' =~ / $blocktype /;
   for (;;) {
      $timeleft = $self->kill_late_children;
      # Compute throttle
      if ($self->{throttle} > 0) {
         my $part = 9 - int(10*$timeleft/$self->{max_time});
	 if ($self->{throttle} > $part) {
	    $self->{max_procs} = $self->{orig_max_procs}; }
	 else {
	    $part = int($self->{throttle}*$self->{max_time}/10);
	    $part = 1-($part-$timeleft)/($self->{max_time}-$part);
	    $part = int($part*$self->{orig_max_procs});
	    $part = 2 if abs($part) < 2; }
      }
      # pick up untimely dead processes
      if (@{$self->{'dead_proc'}}) {
         for (my $i = $#{$self->{'dead_proc'}}; $i >= 0; $i--) {
            $pid = waitpid(${$self->{'dead_proc'}}[$i], WNOHANG);
	    next if $pid <= 0;
	    $self->{status}->{$pid} .= 'D';
	    splice(@{$self->{'dead_proc'}}, $i, 1);
	    $self->{'data'}->{$pid} = undef if $self->{'start_type'} eq 'IPC';
	    $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
            $self->{'cur_proc'}--;
            delete $self->{'proc_tab'}->{$pid};
#	    warn "Reaping killed child: $self->{proc_order}->{$pid}\tpid->$pid\n" if $self->{start_type} eq 'IPC'; # debug
	    return $pid; }
      }
      # Now normal processes
      if ($self->{'start_type'} eq 'IPC') {
         $timeleft = 1 if $blocktype eq 'NOBLOCK';
	 my $client;
	 my $address;

	 eval {  # Implementing timeout on accept call to facilitate NOBLOCK
	    local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
	    alarm $timeleft;
	    $address = accept($client,$self->{'socket'});
	    alarm 0;
	 };

#         my $rvec = '';
#         vec($rvec,fileno($self->{'socket'}),1) = 1;
#         unless (select($rvec,undef,undef,$timeleft)) {
	 unless (defined $address) {
	    return 0 if $blocktype eq 'NOBLOCK';
	    next; }
#         $address = accept($client,$self->{'socket'}); # new line
	 
	 
	 my $id = <$client>;
         unless ($id =~ m/^<Child (\d+)>/) {
	    close $client;
	    warn "Unknown and unexpected connection from irrelevant source\n" if $self->{warning} == 2;
	    next;
	 }
	 my $return_pid = $1;
#	 warn "Getting data from $self->{proc_order}->{$return_pid}\tpid->$return_pid\t\tData in queue: ",
#	      scalar keys %{$self->{data}}, "\n"; # debug
	 substr($id, 1, 0, '/');
	 my $data = '';
	 my $line;
	 while (defined ($line = <$client>)) {
	    last if $id eq $line;
	    $data .= $line; }
	 close $client;
	 if (defined $line) {
	    my $VAR1;
	    eval $data;
	    $self->{'data'}->{$return_pid} = $VAR1;
	    $self->{status}->{$return_pid} .= 'R';
	 } else {
	    $self->{'data'}->{$return_pid} = undef;
	    warn "Unexpected end of child communication, losing data from one job\n" if $self->{warning};
	 }
     	 $pid = waitpid($return_pid, 0); # Blocking mode, but that should be OK
	 $self->{status}->{$pid} .= 'E';
         $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
         $self->{'cur_proc'}--;
         delete $self->{'proc_tab'}->{$pid};
     	 return $pid;
      } else {
         $pid = waitpid(-1, WNOHANG);
         if ($pid > 0) {
	    $self->{status}->{$pid} .= 'E';
            $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
            $self->{'cur_proc'}--;
            delete $self->{'proc_tab'}->{$pid};
            return ($pid);
         } elsif ($blocktype eq 'NOBLOCK') {
            return 0; 
         } elsif ($timeleft) {
	    eval {
	       local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
	       alarm($timeleft);
	       $pid = wait;
	       alarm(0); };
            if (defined $pid and $pid > 0) {
	       $self->{status}->{$pid} .= 'E';
               $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
               $self->{'cur_proc'}--;
               delete $self->{'proc_tab'}->{$pid};
               return ($pid); }
         } elsif ($self->{'max_time'} == -1) {
	    $pid = wait;
	    $self->{status}->{$pid} .= 'E';
            $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
            $self->{'cur_proc'}--;
            delete $self->{'proc_tab'}->{$pid};
            return ($pid);
         }
      }
   }
}


sub wait_all_children {
   my ($self) = @_;
   while ($self->{'cur_proc'}) {
      $self->wait_next_child('BLOCK'); }
}


sub reap {
   my $self = shift @_;
   while ($self->wait_next_child('NOBLOCK')) {}
   my @procs = keys %{$self->{'proc_tab'}};
   foreach my $pid (@procs) {
      unless (kill 0, $pid) {
         $self->{'ready_cpus'}++ if $self->{'last_check'} + 30 > time;
         $self->{'cur_proc'}--;
         delete $self->{'proc_tab'}->{$pid}; } }
}
      

sub kill_late_children {
   my ($self) = @_;
   return 0 if $self->{'max_time'} == -1 or $self->{'cur_proc'} == 0;
   my ($time, $pid, $timeused, $result);
   $time = time;
   $result = -1;
   foreach $pid (keys %{$self->{'proc_tab'}}) {
      $timeused = $time - $self->{'proc_tab'}->{$pid};
      if ($timeused >= $self->{'max_time'}) {
         if (kill('TERM', $pid)) {
	    warn "Killing a process ($pid) exceeding the time limit, likely data loss\n" if $self->{warning};
	    $self->{status}->{$pid} .= 'K';
	 } else {
	    warn "A process ($pid) died without returning data\n" if $self->{warning};
	 }
	 push(@{$self->{'dead_proc'}}, $pid);
         $self->{'proc_tab'}{$pid} = $time;
      } else {
         $result = $timeused if $result < $timeused; } }
   return ($self->{'max_time'} + 1 - $result) if $result != -1;
   return 0;
}


sub kill_all {
   my $pid;
   foreach $pid (keys %{$packagehash->{'proc_tab'}}) {
      kill 'TERM', -$pid; }
   $SIG{'TERM'} = 'IGNORE';
   kill 'TERM', -$$;
   exit;
}


sub return_result {
   my ($self, $reference) = @_;
   my $paddr = sockaddr_in($self->{'socketport'}, INADDR_LOOPBACK);
   my $socket;
#   warn "Child return: ", $self->{workercount}, "\tpid->$$\n"; # debug
   socket($socket, PF_INET, SOCK_STREAM, getprotobyname('tcp')) or
      die "Child: Socket call (initialization) failed: $!";
   my $tries = 0;
   while (1) {
      connect($socket, $paddr) and last;
      sleep 1;
      $tries++;
      if ($tries >= 100) {
         warn "Child: Connect call to parent failed $tries times, dying - data loss: $!\n" if $self->{warning};
	 exit;
      }
      ($tries >= 10 and $self->{warning} == 2) or
         warn "Child: Connect call to parent failed $tries times (still trying): $!\n";
      }
   print $socket "<Child $$>\n", Dumper($reference) , "\n</Child $$>\n" or
      ($self->{warning} == 2 and warn "Child: Streaming data to parent failed: $!\n");
   close ($socket) or
      ($self->{warning} == 2 and warn "Child: Close call on socket failed: $!\n");
   exit;
} 


sub get_all {
   my($self) = @_;
   # Enforce waiting for all children
   while ($self->{'cur_proc'}) {
      $self->wait_next_child('BLOCK'); }
   my @data;
   for (my $i = 0; $i <= $#{$self->{order}}; $i++) {
      $data[$i] = exists $self->{data}->{$self->{order}->[$i]} ? $self->{data}->{$self->{order}->[$i]} : undef;
   }
   return @data;
}


sub get_ready {
   my($self, $block) = @_;
   $block = 'BLOCK' unless defined $block;
   die "Undefined blocking type parameter: $block"
      unless ' BLOCK NOBLOCK ' =~ / $block /;
   my @key;
   while ((@key = keys %{$self->{data}}) == 0) {
      my $childno = scalar keys %{$self->{proc_tab}};
      return $childno if $block eq 'NOBLOCK';
      if ($childno) {
         $self->wait_next_child('BLOCK');
      } else {
         return undef;
      }
   }
   my $k = shift @key;
   my $reference = $self->{data}->{$k};
   delete $self->{data}->{$k};
   return $reference;
}


sub get_next {
   my($self, $block) = @_;
   $block = 'BLOCK' unless defined $block;
   die "Undefined blocking type parameter: $block"
      unless ' BLOCK NOBLOCK ' =~ / $block /;
   my $pidkey = ${$self->{order}}[0];
   return 0 unless defined $pidkey;
   my $turn = 0;
   for (;;) {
      if (exists $self->{data}->{$pidkey}) {
         my $reference = $self->{data}->{$pidkey};
         delete $self->{data}->{$pidkey};
	 shift @{$self->{order}};
	 return $reference; }
      elsif ($block eq 'NOBLOCK') {
         return scalar keys %{$self->{proc_tab}}; }
      if (exists $self->{proc_tab}->{$pidkey}) {
         $self->wait_next_child('BLOCK');
      } elsif ($turn == 0) {
         return undef if 0 == scalar keys %{$self->{proc_tab}};
	 $self->wait_next_child('BLOCK');
	 $turn++;
      } else {
         return undef;
      }
   }
}


sub get_status {
   my ($self, $level) = @_;
   $level = 0 unless defined $level;
   my $line = "Total started jobs: $self->{workercount}\tCurrent jobs: $self->{cur_proc}\t";
   my ($endjobs, $killedjobs, $diedjobs, $reusepid, $returndata) = (0,0,0,0,0);
   # Status fields; E = ended normally, D = died unexpected, K = killed by parent
   # R = returned data, P = process id reuse
   my @list;
   my $href = $self->{status};
   while (my($pid, $val) = each %{$href}) {
      $reusepid += $val =~ tr/P/P/;
      $endjobs += $val =~ tr/E/E/;
      $returndata += $val =~ tr/R/R/;
      if ($val =~ tr/D/D/) {
         $diedjobs++;
	 my $out = substr($val, 0, index($val, ' ')) . "\tpid: $pid\t";
	 if ($val =~ tr/K/K/) {
            $killedjobs++;
	    $out .= "Killed due to timeout\n";
	 } else {
	    $out .= "Died - not related to child manager\n";
	 }
	 push(@list, $out);
      }
   }
   $line .= "Ended jobs: $endjobs\tDied: $diedjobs (killed $killedjobs)\t";
   $line .= "Returned data: $returndata\t" if $self->{start_type} eq 'IPC';
   $line .= "Process id reused: $reusepid\n";
   return $line if $level == 0;
   if (@list) {
      @list = sort {$a<=>$b} @list;
      $line .= 'Job ' . join('Job ', @list);
   }
   return $line if $level == 1;
   $line .= "All seems fine\n" if ($self->{workercount} == $self->{cur_proc} + $endjobs);
   $line .= "Some jobs did not return data allthough they ended normally - technically unexpected\n"
      if $self->{start_type} eq 'IPC' and $endjobs != $returndata;
   if ($diedjobs) {
      if ($killedjobs) {
         $line .= $diedjobs == $killedjobs ? "All bad job death is due to timeout\n" :
                                             "Some jobs die from timeout, other from own internal problems\n";
      } else {
         $line .= "All job death is due to own internal problems\n";
      }
   }
   return $line;
}


sub examine_load {
   my($self) = @_;
   my($load, $cpu, @top);
   $load = time;
   return if $load < ($self->{'last_check'} + 30);
   $self->{'last_check'} = $load;
   if ($ENV{'HOSTTYPE'} ne 'iris4d') {
      warn "Probably not an SGI/IRIX system, setting ready processors to 0\n";
      $self->{'ready_cpus'} = 0;
      return; }
   @top = `top -b -u`;
   if ($top[0] =~ /load averages:\s+(\d+\.\d+)\s+/) {
      $load = $1; }
   else {
      warn "Can't recognize top output, setting ready processors to 0\n";
      $self->{'ready_cpus'} = 0;
      return; }
   if ($top[2] =~ /^(\d+)\s+CPUs:/) {
      $cpu = $1; }
   else {
      warn "Can't recognize top output, setting ready processors to 0\n";
      $self->{'ready_cpus'} = 0;
      return; }
   $load = int($cpu - $load);
   $load = 0 if $load < 0;
   $self->{'ready_cpus'} = $load;
}

1;
