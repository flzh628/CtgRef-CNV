#!/usr/bin/perl -w
use strict;
use Math::Round;

unless(@ARGV==2){
	print STDOUT "\n\tAuthor: Fengli Zhao\tAGIS, CAAS  flzh628\@gmail.com\n";
	print STDOUT "\tVersion: 2.0\t\tDate: 2017/10/28  2017/3/9 \n";
	print STDOUT "\tThis script is designed to correct the depth of Wins from the result of 'GC-Dpth-slidingWin.pl' based on GC content.\n";
	print STDOUT "\tUsage: perl $0 GC-Dpth-slidingWin-result GC-correct-result\n\n";
	exit;
}

open IN,"$ARGV[0]" or die $!;
my (%Dpth,%Gc,$mn,$tot_dep,$dpav);
while(<IN>){
	chomp;
	next if $.==1;
	my ($chr,$id,$gc,$depth)=split;
	$Dpth{$chr}{$id}=$depth;
	$Gc{$chr}{$id}=$gc;
	$mn++;
	$tot_dep+=$depth;
}
close IN;
$dpav=sprintf "%.2f",$tot_dep/$mn;

my (%Dp_gc,%Dp_chr,%Pos_chr,%depthgc,%numgc,%AverDp);
for my $chr(sort keys %Gc){
	for my $pos(sort keys %{$Gc{$chr}}){
		my $gc_n=round($Gc{$chr}{$pos});
		$depthgc{$gc_n}+=$Dpth{$chr}{$pos};
		$numgc{$gc_n}++;
		$Dp_chr{$chr}+=$Dpth{$chr}{$pos};
		$Pos_chr{$chr}++; 
	}
	$AverDp{$chr}=sprintf "%.2f",$Dp_chr{$chr}/$Pos_chr{$chr};
}

for my $gc(sort keys %depthgc){
	$Dp_gc{$gc}=sprintf "%.2f",$depthgc{$gc}/$numgc{$gc};
}

open OUT1,">$ARGV[1].supp" or die $!;#20171028
print OUT1 "Genome-AverDpth\t$dpav\n";
for my $chr(sort keys %AverDp){
	print OUT1 "$chr-AverDpth\t$AverDp{$chr}\n";
}
print OUT1 "\nGC\tGC_averdp\n";
for my $igc(sort {$a <=> $b} keys %Dp_gc){
	print OUT1 "$igc\t$Dp_gc{$igc}\n";
}
close OUT1;
open OUT,">$ARGV[1]" or die $!;
print OUT "Chr\tWin\tGC\tDepth\tDp_Cor\n";
for my $chr(sort keys %Dpth){
	for my $pos(sort {$a <=> $b} keys %{$Dpth{$chr}}){
		print OUT "$chr\t$pos\t$Gc{$chr}{$pos}\t$Dpth{$chr}{$pos}\t";
		my $gc_int=round($Gc{$chr}{$pos});
		$Dp_gc{$gc_int}=0.001 if $Dp_gc{$gc_int}==0;
		my $depth=($dpav*$Dpth{$chr}{$pos})/$Dp_gc{$gc_int};
		my $depthn=sprintf "%.2f",$depth;
		print OUT "$depthn\n";
	}
}
close OUT;
