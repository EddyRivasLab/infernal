#! /usr/bin/perl 

# The top level script that runs a rmark benchmark.
#
# Usage:
#   ./rmark-master.pl <execdir> <resultdir> <benchmark prefix> <benchmark script>
# 
# <execdir>: Directory for finding executables to be benchmarked. This
#    is passed on to the <benchmark script> without modification.
#
# <modeldir>: Directory for finding model files for the benchmark. This
#    is passed on to the <benchmark script> without modification.
#    Some programs don't require model files (e.g. BLAST), but they're
#    benchmark scripts will still take a <modeldir> argument.
#
# <resultdir>: A directory for holding all <ncpu> temporary files
#   created by the benchmark. This name should be short and unique;
#   it will also be used to construct job names on the cluster, as
#   <resultdir.$i>.
#
# <optsfile>: File with options for the search program. 
#   For example, a non-filtered cmsearch 1.0 run would have a single
#   line in this file:
#   '--fil-no-hmm --fil-no-qdb"
#
# <benchmark prefix>:  The script will look for files <prefix>.tbl,
#    <prefix>.msa, <prefix>.fa, and <prefix>.pos, defining a RMARK
#     benchmark set. If -P option is used, a positive-only benchmark is
#     being conducted, <prefix>.pfa and <prefix>.ppos will be used
#     instead of <prefix>.fa and <prefix>.pos.
#
# <benchmark script>: This script is executed for each process, 
#    on an appropriately constructed subset of the benchmark queries.
#
#    It must take the following arguments:
#    <execdir> <modeldir> <resultdir> <optsfile> <tblfile> <msafile> <posfile> <fafile> <outfile>
#
# Command-line options:
# -P     : run a positive-only benchmark, only positive sequences will be searched
# -N <n> : specify <N> processes be used for all families, default is to do one process per family
# -M <n> : pass -M <n> onto search module, telling it to run MPI with <n> <= 8 processors.
#          only valid if --mpi exists in the $optsfile
#
# Examples of rmark benchmark:
#   ./rmark-master.pl ../src/ cmsearch-results cmsearch-df.opts rmark3 ./rmark-cmsearch
#   ./rmark-master.pl -P -N 1 ../src/ cmsearch-po-results cmsearch-df.opts rmark3 ./rmark-cmsearch

use Getopt::Std;
getopts('PN:M:');
$do_posonly = 0;
$ncpu_set  = 0;
$ncpus = 0;
$mpi_nprocs = 8; #default, we'll only use it if --mpi exists in <optsfile>
if (defined $opt_P) { $do_posonly = 1; }
if (defined $opt_N) { $ncpu_set = 1; $ncpus = $opt_N; }
if (defined $opt_M) { 
    $mpi_nprocs = $opt_M; 
    if($mpi_nprocs < 2 || $mpi_nprocs > 8) { die "ERROR, with -M <n>, <n> must be between 2 and 8"; }
}

$usage =  "Usage: perl rmark-master.pl\n\t<executable dir>\n\t<model dir>\n\tresult dir, for output, must not exist>\n\t";
$usage .= "<options file, 1 line>\n\t<benchmark prefix>\n\t<benchmark script>\n\n";
$options_usage  = "Options:\n\t";
$options_usage .= "-P     : run a positive-only benchmark, only positive sequences will be searched\n\t";
$options_usage .= "-N <n> : use <n> processors, default is to use one per model\n\t";
$options_usage .= "-M <n> : pass -M <n> onto search module, telling it to run MPI with <n> <= 8 processors\n\t";
$options_usage .= "         only valid if --mpi exists in <optsfile>\n\n";

if(@ARGV != 6)
{
    print $usage;
    print $options_usage;
    exit();
}

($execdir, $modeldir, $resultdir, $optsfile, $benchmark_pfx, $rmark_script) = @ARGV;

$tbl      = "$benchmark_pfx.tbl";
$msafile  = "$benchmark_pfx.msa";
$fafile   = "$benchmark_pfx.fa";
$pfafile  = "$benchmark_pfx.pfa";
$posfile  = "$benchmark_pfx.pos";
$pposfile = "$benchmark_pfx.ppos";

if ($do_posonly) { 
    $fafile  = $pfafile;
    $posfile = $pposfile;
}
if (-e $resultdir) { die("$resultdir exists");}
system("mkdir $resultdir");

# Suck in the master table
open(BENCHMARK_TBL, $tbl) || die;
$n = 0;
while (<BENCHMARK_TBL>) 
{
    ($msaname[$n], $pid, $L, $nseq) = split;
    $alen{$msaname[$n]} = $L;
    $n++;
}
close BENCHMARK_TBL;

if(! $ncpu_set) { $ncpu = $n; }

# Sort it by alen - this helps load balance.
sub by_alen { $alen{$b} <=> $alen{$a} }
@sorted_msaname = sort by_alen @msaname;

# Create <ncpu> subtables.
for ($i = 0; $i < $n; $i++)
{
    $subtbl[$i % $ncpu] .= $sorted_msaname[$i];
    $subtbl[$i % $ncpu] .= "\n";
}

# Output the <ncpu> subtables
for ($i = 0; $i < $ncpu; $i++)
{
    open(SUBTBL, ">$resultdir/tbl.$i") || die ("Failed to create $resultdir/tbl.$i");
    print SUBTBL $subtbl[$i];
    close SUBTBL;
}

# Determine if we're using MPI or not
$do_mpi = 0;
open(OPTS, $optsfile) || die "couldn't open options file $optsfile"; 
$searchopts = <OPTS>;
if($searchopts =~ m/\-\-mpi/) { $do_mpi = 1; }
close(OPTS);
chomp $searchopts;

# Submit all the individual rmark jobs
for ($i = 0; $i < $ncpu; $i++)
{
    if($do_mpi) { # turn exclusivity on, so we get all processors on our node, to run MPI with
	system("qsub -V -cwd -b y -N $resultdir.$i -j y -o $resultdir/tbl$i.sge -l excl=true '$rmark_script -M $mpi_nprocs $execdir $modeldir $resultdir $optsfile $resultdir/tbl.$i $msafile $posfile $fafile $resultdir/tbl$i.out'");
    }
    else { 
	system("qsub -V -cwd -b y -N $resultdir.$i -j y -o $resultdir/tbl$i.sge '$rmark_script $execdir $modeldir $resultdir $optsfile $resultdir/tbl.$i $msafile $posfile $fafile $resultdir/tbl$i.out'");
    }
}

