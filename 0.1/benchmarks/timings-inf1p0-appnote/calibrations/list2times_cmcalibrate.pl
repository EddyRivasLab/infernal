if(scalar(@ARGV) != 1) { printf("usage: perl list2times_cmcalibrate.pl <list file>\n"); exit(1); }

($list_file) = @ARGV;

open(LIST, $list_file);
while($line = <LIST>) { 
    chomp $line;
    if($line =~ m/\w/) { 
	($key, $file, $nprocs) = split(/\s+/, $line);
	if(($key eq "") || ($file eq "") || ($nprocs eq "")) { die("ERROR list line: $line is invalid.\n"); }
	if(!(-e ($file))) { die("ERROR, search file $file does not exist."); }
	$time = `tail -1 $file`;
	if($time =~ /\# CPU time:\s*(\S+)\s+(\S+)\s+(\S+)\s+Elapsed:\s+(\S+)/) { 
	    $elapsed = $4;
	    ($hours, $minutes, $seconds) = split(":", $elapsed); 
	    $total_hours = $nprocs * (($hours) + ($minutes/60.) + ($seconds / 3600.));
	}
	else { die("ERROR, couldn't read time line $time\n"); }
	printf("%-25s %8.2f hours\n", $key, ($total_hours));
    }
}
	   
