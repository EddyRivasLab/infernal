# Personal library of Perl functions.
#
# Thu Feb 22 09:18:39 1996 -- FASTA reading/writing added.



package SRE_perlstuff;

require "importenv.pl";

# Function: tempname
#
# Returns a unique temporary filename. 
#
# Should be robust. Uses the pid as part of the temp name
# to prevent other processes from clashing. A two-letter
# code is also added, so a given process can request
# up to 676 temp file names (26*26). An "sre" code is
# also added to distinguish these temp files from those
# made by other programs.
#
# Returns nothing if it fails to get a temp file name.
#
# Normally puts temp files to /tmp. This directory can
# be overridden by an environment variable TMPDIR
#
sub main'tempname {
    local ($dir, $name);
    if ($TMPDIR) { $dir = $TMPDIR; } else {$dir = "/tmp";}

    foreach $suffix ("aa".."zz") {
	$name = "$dir/sre$suffix$$";
        if (! (-e $name)) { 
            open(TMP,">$name") || die; # Touch it to reserve it.
	    close(TMP);
            return "$name"; 
        }
    }                           
}


# Function: alistat(file)
#
# Calls 'alistat'.
# Returns the average, high/low % identities, and length
# for an the alignment.
# Dies if 'alistat' fails.
#
sub main'alistat {
    local($file) = @_;
    local($output, $nseq, $nres, $minlen, $maxlen, $avlen);
    local($alilen, $id, $high, $low, $outlier);
    $output = `alistat $file`;
    die ("alistat failed; died") unless ($? == 0);
    if ($output =~ /Number of sequences:\s+(\d+)/)  { $nseq   = $1; }
    if ($output =~ /Total \# of residues:\s+(\d+)/) { $nres   = $1; }
    if ($output =~ /Smallest:\s+(\d+)/)             { $minlen = $1; }
    if ($output =~ /Largest:\s+(\d+)/)              { $maxlen = $1; }
    if ($output =~ /Average length:\s+(\S+)/)       { $avlen  = $1; }
    if ($output =~ /Alignment length:\s+(\d+)/)     { $alilen = $1; }
    if ($output =~ /Average identity:\s+(\d+)%/)    { $id     = $1; }
    if ($output =~ /Most related pair:\s+(\d+)%/)   { $high   = $1; }
    if ($output =~ /Most unrelated pair:\s+(\d+)%/) { $low    = $1; }
    if ($output =~ /Most distant seq:\s+(\d+)/)     { $outlier = $1; }
    ($nseq, $nres, $minlen, $maxlen, $avlen, $alilen, $id, $high, $low, $outlier);
}


# Function: seqstat(file)
#
# Calls 'seqstat'.
# Returns the number of sequences in the file,
# and their maximum and minimum length, and their avg. len.
# Dies if 'seqstat' fails.
#
sub main'seqstat {
    local($file) = @_;
    local($output, $nseq, $fromlen, $tolen, $avlen);
    $output = `seqstat $file`;
    die ("seqstat failed; died") unless ($? == 0);
    if ($output =~ /Number of sequences:\s+(\d+)/)      {$nseq    = $1; }
    if ($output =~ /Smallest:\s+(\d+)/)                 {$fromlen = $1; }
    if ($output =~ /Largest:\s+(\d+)/)                  {$tolen   = $1; }
    if ($output =~ /Average length:\s+(\S+)/)           { $avlen = $1;}
    ($nseq, $fromlen, $tolen, $avlen);
}

# Function: gdfcount(GDF file)
#
# Counts nseqs, ndomains, and nresidues covered in a 
# GDF format file.
# Usage: ($nseqs, $ndomains, $nresidues) = &gdfcount("foo.gdf");
#
sub main'gdfcount {
    local($file) = @_;
    local($nseqs, $ndomains, $nresidues, %hmmhit);

    $nseqs = $ndomains = $nresidues = 0;
    open(TMP,"$file") || die("file open failed for $file");
    while (<TMP>) {
        if (/^\s*\S+\s+(\d+)\s+(\d+)\s+(\S+)/) {
            $ndomains++;
            $nresidues += $2 - $1 + 1;
            if (! $hmmhit{$3}) { $nseqs++; }
            $hmmhit{$3} = 1;
        }
    }
    close(TMP);
    ($nseqs, $ndomains, $nresidues);
}

##########################################
# FASTA reading/writing functions
########################################## 

# Function: open_fasta($filename)
# 
# Opens a FASTA file. Keeps the open file handle 
# as a package global variable. (I.e. don't open
# more than one FASTA file at a time.)
#
sub main'open_fasta {
    local($fname) = @_;
    open(FASTA, $fname) || die("Failed to open FASTA file $fname\n");
    $SavedLine = "";
    1;
}
sub main'close_fasta {
    close(FASTA);
    1;
}

# Function: read_fasta(*name, *desc, *seq)
#
# Sequential read of FASTA file.
# Returns 0 when there are no more sequences.
#
sub main'read_fasta {
    local(*name, *desc, *seq) = @_;

    while (($SavedLine =~ /^>/) || ($SavedLine = <FASTA>))
    {
	if (($name, $desc) = ($SavedLine =~ /^>\s*(\S+)\s*(.*)$/))
	{
	    $seq   = "";
	    while ($SavedLine = <FASTA>)
	    {
		if ($SavedLine =~ /^>/) { last; }
		$SavedLine =~ s/[ \n\t]//g;	# strip whitespace
		$seq = $seq.$SavedLine;
	    }
	    return 1;
	}
    }
    0;
}

# Function: write_fasta(*FILEHANDLE, $name, $desc, $seq)
#
sub main'write_fasta {
    local(*FASTA, $name, $desc, $seq) = @_;
    local($pos, $line, $length);

    $length = length($seq);

    print FASTA ">$name $desc\n";
    for ($pos = 0; $pos < $length; $pos += 60)
    {
	$line = substr($seq,$pos,60);
	print FASTA $line, "\n";
    }
    1;
}


1;				# This "1" must be here or require fails.
