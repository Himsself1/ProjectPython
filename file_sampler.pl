#! /usr/bin/perl
use warnings;

open IN , "<ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf";
open OUT, ">chr22_3000_lines.vcf";
print OUT "##casjkfcja","\n","##casjkfcja","\n","##casjkfcja","\n","##casjkfcja","\n";
$no=0;
while (<IN>) {
	if (/^##/) {
		next;
	} 
	else {
	
		if ($no<3001) {
			print OUT $_;
			$no++
		}
		else {
			exit;
		}
	}
}