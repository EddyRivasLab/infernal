#! /usr/bin/perl

package i1;

sub ParseTbl {
    my ($tblfile)    = @_;
    my (@fields);

    $ntbl     = 0;
    @tname    = ();
    @tacc     = ();
    @qname    = ();
    @qacc     = ();
    @model    = ();
    @mfrom    = ();
    @mto      = ();
    @sfrom    = ();
    @sto      = ();
    @strand   = ();
    @trunc    = ();
    @pass     = ();
    @hitgc    = ();
    @hitbias  = ();
    @hitsc    = ();
    @hitE     = ();
    @tdesc    = ();
    
    if (! open(TBLFILE, $tblfile)) { print "FAIL: couldn't open table file"; exit 1 ; }
    while (<TBLFILE>)
    {
	if (/^\#/) { next; }
	chop;
	@fields = split(' ', $_, 17);

	$tname[$ntbl]     = $fields[0];
	$tacc[$ntbl]      = $fields[1];
	$qname[$ntbl]     = $fields[2];
	$qacc[$ntbl]      = $fields[3];
	$model[$ntbl]     = $fields[4];
	$mfrom[$ntbl]     = $fields[5];
	$mto[$ntbl]       = $fields[6];
	$sfrom[$ntbl]     = $fields[7];
	$sto[$ntbl]       = $fields[8];
	$strand[$ntbl]    = $fields[9];
	$trunc[$ntbl]     = $fields[10];
	$pass[$ntbl]      = $fields[11];
	$hitgc[$ntbl]     = $fields[12];
	$hitbias[$ntbl]   = $fields[13];
	$hitsc[$ntbl]     = $fields[14];
	$hitE[$ntbl]      = $fields[15];	
	$tdesc[$ntbl]     = $fields[16];
	$ntbl++;
    }
    close TBLFILE;
    1;
}

1;
