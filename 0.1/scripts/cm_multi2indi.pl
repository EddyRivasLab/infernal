$cmi = 0;
while($line = <>) {
    if($line =~ /^INFERNAL/) { 
	$inf_line = $line;
    }
    elsif($line =~ /^NAME\s+(\S+)/) { 
	chomp $line; 
	$new_file = $1 . ".cm";
	open(OUT, ">" . $new_file); 
	print OUT $inf_line;
	print OUT $line . "\n";
    }
    elsif($line =~ m/^\/\//) { 
	print OUT $line;
	close(OUT);
    }
    else {
	print OUT $line;
    }
}

