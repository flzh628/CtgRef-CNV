#!/usr/bin/perl -w
#Author: Fengli Zhao	AGIS, CAAS  flzh628@gmail.com	2017/4/15
use strict;
use Statistics::Basic qw(:all);

my $ref=shift or die "\n\tUsage: $0 Ref Dpth-file.gz Win-width Res-Dir slidWinDp-res \n\n";
my $dpth=shift or die $!;
my $win_width=shift or die $!;
my $dir=shift or die $!;
my $output=shift or die $!;

if (!-e "$ref.fai"){
        my $cmd="samtools faidx $ref";
        system($cmd);
}
if (!-e $dir){
	mkdir $dir,0755;
}

open DP,"zcat $dpth |" or die $!;
my %Dp;
while(<DP>){
	chomp;
	my ($chr,$pos,$dp)=split;
	my $win=sprintf "%d",$pos/$win_width;
	my $dpnew=sprintf "%.4f",$dp/$win_width;
	$Dp{$chr}{$win}+=$dpnew;
}
close DP;

open OUT1,">$dir/$output" or die $!;
print OUT1 "CHR\tWindow\tDepth\n";
for my $chr(sort keys %Dp){
	my $chrlen=Chr_len($ref,$chr);
	my $w2=sprintf "%d",$chrlen/$win_width;
	for my $win(0 .. $w2){
		print OUT1 "$chr\t$win\t";
		if (defined $Dp{$chr}{$win}){
			my $dpi=sprintf "%.4f",$Dp{$chr}{$win};
			print OUT1 "$dpi";
		}else{
			print OUT1 "0";
		}
		print OUT1 "\n";
	}
}
close OUT1;

open IN,"$dir/$output" or die $!;
my (%DpArr,$num,$totdp);
while(<IN>){
        chomp;
        next if $.==1;
        my ($chr,$dp)=(split)[0,2];
	push @{$DpArr{$chr}},$dp;
        $num++;
        $totdp+=$dp;
}
close IN;

open OUT2,">$dir/$output-ChrDp" or die $!;
my $averdp=sprintf "%.2f",$totdp/$num;
print OUT2 "Chr\tAvDp\tMedian\tSd\tVar\tGenomic: $averdp\n";
for my $chr (sort keys %DpArr){
        my $median=median(@{$DpArr{$chr}});
        my $mean=mean(@{$DpArr{$chr}});
        my $var=variance(@{$DpArr{$chr}});
        my $sd=stddev(@{$DpArr{$chr}});
        print OUT2 "$chr\t$mean\t$median\t$sd\t$var\n";
}
close OUT2;
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++Get the length of a Chr++++++++++++++++++++
sub Chr_len{
	my ($ref,$chr)=@_;
	my $seq=`samtools faidx $ref $chr`;
	my $seq_new=(split/\n/,$seq,2)[1];
	$seq_new=~ s/\n//g;
	my $len=length($seq_new);
	return $len;
}
