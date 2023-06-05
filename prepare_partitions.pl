#!/usr/bin/perl
use strict;
use warnings;

# make sure to source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
# prior to running this script

opendir my $dir, "./MPAS_O_Shallow_Water/Convergence2/CoastalKelvinWave/" or die "Cannot open directory: $!";

my @graphfiles = grep { /^graph_\d+x\d+\.info$/ } readdir($dir);
closedir($dir);

# mkdir my $pdir, "./graphpartitions/";

foreach my $graphfile (@graphfiles) {
	$graphfile =~ /^graph_(\d+)x\d+\.info$/;
	# mkdir my $resolutiondir, "${pdir}/${1}x${1}/";
	opendir my $resolutiondir, "graphpartitions/${1}x${1}/";# or die "problem creating directory!";
	for my $p (0..12) {
		my $npieces = 2**$p;
		# print $npieces . $graphfile . $1 . "\n";
		system("/global/common/software/e3sm/anaconda_envs/base/envs/e3sm_unified_1.8.0_nompi/bin/gpmetis ${npieces} ../../${graphfile}")
	}
	closedir $resolutiondir;
}



