#! /usr/bin/perl 
#
# Given a list file with 3-tuples for each series to plot
# create a pdf of a ROC curve.
# 
# Format of list file:
# <series_name> <root> <color>
#
# <root>/<root>.xy, <root>/<root>.mer and <root>/<root>.time must exist.
#
# Usage:  perl rmark-Rroc.pl <listfile> <name for pdf>
# Example usage as pipe into R:
# > perl rmark-Rroc.pl rmark-2.list rmark-2.pdf 1 rmark-ROC
my $usage = "Usage: perl rmark-rocR.pl <listfile> <pdfname> <1/0 yes/no draw error-bars> <plot title>\n";
$usage .= "\nFormat of list file:\n\t<series_name> <root> <color>\n\n";
$usage .= "\nExample:\n\tinf1p02-df r2-i1p02-df red\n\n";
$usage .= "<root>/<root>.xy, <root>/<root>.mer and <root>/<root>.time must exist\n\n";

if(scalar(@ARGV) != 4) { printf("$usage"); exit(1); }
($listfile, $pdf, $do_errorbars, $main) = @ARGV;
$n = 0;

@R = ();

$xlabel = "errors per query";
$ylabel = "fractional coverage of positives";

push(@R, "ylimit<-c(0,1)\n");
push(@R, "pdf(\"$pdf\", height=8.5, width=11)\n");

open(LIST, $listfile) || die "ERROR, could not open $listfile"; 
@nametimemerA = ();
while(<LIST>) { 
    chomp $_;
    if(m/\S+/) {
	@fields  = split(' ', $_, 3);
	if(scalar(@fields) != 3) { die "ERROR, format of $listfile invalid"; }
	($name, $root, $color) = @fields;
	$n++;
	@xA = ();
	@yA = ();
	@dy1A = ();
	@dy2A = ();
	$xyfile = $root . "/" . $root . ".xy";
	$merfile = $root . "/" . $root . ".mer";
	$timefile = $root . "/" . $root . ".time";
	if(! -e $xyfile) { die "ERROR, $xyfile does not exist"; }
	if(! -e $merfile) { die "ERROR, $merfile does not exist"; }
	if(! -e $timefile) { die "ERROR, $timefile does not exist"; }

	# process time file
	open(TIME, $timefile) || die "ERROR, could not open time file $timefile";
	while($line = <TIME>) { 
	    chomp $line;
	    if($line =~ s/^total\s+//) { 
		$line =~ s/\s+.+$//;
		$time = $line;
	    }
	}
	close(TIME);

	# process mer file
	open(MER, $merfile) || die "ERROR, could not open mer file $merfile";
	while($mer = <MER>) { 
	    if($mer =~ s/^\s*\*summary\*\s+//) { 
		$mer =~ s/\s+$//; 
		push(@nametimemerA,  "\"" . $name  . " " . $time . "h MER: " . $mer . "\"");
	    }
	}
	close(MER);
       
	push(@colorA,     "\"" . $color . "\"");
	open(XY,   $root . ".xy")   || die "ERROR, could not open xy file $root.xy";
	while(<XY>) { 
	    if(/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) { 
		($x, $y, $dy1, $dy2) = ($1, $2, $3, $4);
		push(@xA, $x);
		push(@yA, $y);
		push(@dy1A, $y + $dy1);
		push(@dy2A, $y - $dy2);
	    }
	}
	push(@R, "x"   . $n . "<-c" . return_vec_line(\@xA)   . "\n");
	push(@R, "y"   . $n . "<-c" . return_vec_line(\@yA)   . "\n");
	push(@R, "dy1" . $n . "<-c" . return_vec_line(\@dy1A) . "\n");
	push(@R, "dy2" . $n . "<-c" . return_vec_line(\@dy2A) . "\n");
	if($n == 1) { 
	    push(@R, "plot(x$n, y$n,   type=\"l\", log=\"x\", ylim=ylimit, col=\"$color\", main=\"$main\", xlab=\"$xlabel\", ylab=\"$ylabel\")\n");
	}
	else { 
	    push(@R, "points(x$n, y$n, type=\"l\", col=\"$color\")\n");
	}
	if($do_errorbars) {
	    push(@R, "points(x$n, dy1$n, type=\"l\", lty=2, lwd=0.4, col=\"$color\")\n");
	    push(@R, "points(x$n, dy2$n, type=\"l\", lty=2, lwd=0.4, col=\"$color\")\n");
	}
	push(@R, "axis(2, labels=FALSE, at=c(0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0))\n");
	push(@R, "axis(4, labels=FALSE, at=c(0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0))\n");
	close(XY);
    }
}
close(LIST);

push(@R, "legend(0.02, 0.25, c" . return_vec_line(\@nametimemerA) . ", lty=1, cex=0.8, col=c" . return_vec_line(\@colorA) . ", text.col=c" . return_vec_line(\@colorA) . ")\n");
push(@R, "dev.off()\n");

$Rinput = join("", @R);
open(OUT, ">rmark-Rroc.tmp") || die "ERROR, couldn't open rmark-Rroc.tmp for writing";
print OUT $Rinput;
close(OUT);

$status = system("cat rmark-Rroc.tmp | R --vanilla");
if ($status != 0) { die "FAILED: cat rmark-Rroc.tmp | R --vanilla"; }
#unlink "rmark-Rroc.tmp";

sub return_vec_line
{
    ($arr_ref) = $_[0];
    $return_line = "(";
    for($i = 0; $i < (scalar(@{$arr_ref})-1); $i++) { $return_line .= $arr_ref->[$i] . ", "; }
    $return_line .= $arr_ref->[(scalar(@{$arr_ref}-1))] . ")";
    return $return_line;
}
