#!/usr/bin/perl 
$total_seconds = 0;
while($line = <>) { 
    chomp $line;
    if($line =~ /CPU time\:\s+\S+\s+\S+\s+(\d+)\:(\d+)\:(\S+)\s+/) {
	($hours, $minutes, $seconds) = ($1, $2, $3);
	$total_seconds += ($hours * 3600) + ($minutes * 60) + $seconds;
    }
}
printf("%10.2f hours == %10.2f minutes == %10.2f seconds\n", 
       $total_seconds / 3600.,
       $total_seconds / 60.,
       $total_seconds);
