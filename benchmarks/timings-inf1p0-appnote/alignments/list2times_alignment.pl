if(scalar(@ARGV) != 1) { printf("usage: perl list2times_alignments.pl <list file>\n"); exit(1); }

($list_file) = @ARGV;

open(LIST, $list_file);
while($line = <LIST>) { 
    chomp $line;
    if($line =~ m/\w/) { 
	($key, $file, $nseq) = split(/\s+/, $line);
	if(($key eq "") || ($file eq "") || ($nseq eq "")) { die("ERROR list line: $line is invalid.\n"); }
	if(!(-e ($file))) { die("ERROR, search file $file does not exist."); }
	$time = `tail -1 $file`;
	if($time =~ /\# CPU time:\s*(\S+)\s+(\S+)\s+(\S+)\s+Elapsed:\s+(\S+)/) { 
	    $elapsed = $4;
	    ($hours, $minutes, $seconds) = split(":", $elapsed); 
	    $sec_per_seq = (($hours * 3600.) + ($minutes * 60.) + $seconds) / $nseq;
	}
	else { die("ERROR, couldn't read time line $time\n"); }
	printf("%-25s %8.2f sec/seq\n", $key, ($sec_per_seq));
    }
}
	   
