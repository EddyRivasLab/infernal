if(scalar(@ARGV) != 2) { printf("usage: perl list2times.pl <list file> <size in Mb of search (2*dbsize)>\n"); exit(1); }

($list_file, $searchsize) = @ARGV;

open(LIST, $list_file);
while($line = <LIST>) { 
    chomp $line;
    if($line =~ m/\w/) { 
	($key, $file, $nprocs) = split(/\s+/, $line);
	if(($key eq "") || ($file eq "")) { die("ERROR list line: $line is invalid.\n"); }
	if(($nprocs eq "") || ($file eq "")) { die("ERROR list line: $line is invalid.\n"); }
	if(!(-e ($file))) { die("ERROR, search file $file does not exist."); }
	$time = `tail -1 $file`;
	if($time =~ /\# CPU time:\s*(\S+)\s+(\S+)\s+(\S+)\s+Elapsed:\s+(\S+)/) { 
	    $elapsed = $4;
	    ($hours, $minutes, $seconds) = split(":", $elapsed); 
	    $total_min = (60. * $hours) + $minutes + ($seconds / 60.);
	}
	else { die("ERROR, couldn't read time line $time\n"); }
	$min_per_search_Mb = ($total_min / $searchsize) * $nprocs;
	printf("%-25s %8.2f min/Mb\n", $key, ($min_per_search_Mb));
    }
}
	   
