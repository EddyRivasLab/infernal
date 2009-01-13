while(<>)
{
    if(/total:search_runtime\(secs_elapsed\): (\S+)/)
    {
	$secs{$fam} += $1;
	$total_secs += $1;
    }
    elsif(/(\S+):search_runtime\(secs_elapsed\): (\S+)/)
    {
	$fam = $1;
    }
}
foreach $fam (sort keys(%secs))
{
    $hrs{$fam} = $secs{$fam} / 3600;
    $total_hrs += $hrs{$fam};
    printf("$fam:search_runtime(secs_elapsed): $secs{$fam}\n");
    printf("$fam:search_runtime(hrs_elapsed): %.4f\n", $hrs{$fam});
}
printf("total:search_runtime(secs_elapsed): $total_secs\n");
printf("total:search_runtime(hrs_elapsed): %.4f\n", $total_hrs);
