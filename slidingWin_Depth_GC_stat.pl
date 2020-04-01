#!/usr/bin/perl -w
#Author: Fengli Zhao   2017/3/8  2017/3/10  20190905

use strict;

my $ref=shift or die "\n\tUsage: $0 Ref-Genome Bam Win-width Step Chr-Array (default: all the Chrs) >Res-slidWin-GC-Dpth \n\n";
my $bam=shift or die $!;
my $win_width=shift or die $!;
my $step=shift or die $!;
#my @chrarray=@ARGV or die $!;
my @chrarray;

if (@ARGV){
	@chrarray=@ARGV;
}else{
	open REF,"$ref" or die $!;
	$/=">";
	<REF>;
	while(<REF>){
		chomp;
		my $chrid=(split/\n/,$_,2)[0];
		push @chrarray,$chrid;
	}
	close REF;
	$/="\n";
}

if (!-e "$ref.fai"){
	my $cmd="samtools faidx $ref";
	system($cmd);
}
if (!-e "$bam.bai"){
	my $cmd="samtools index $bam";
	system($cmd);
}

print "Chr\tWindow\tGC_%\tDpth_aver\n";
foreach my $chr(@chrarray){
	my $chrlen=Chr_len($ref,$chr);
	for (my $i=1;$i<$chrlen-$win_width+$step;$i+=$step){ # 20190905: For the first win, $i=1 (Wins:1-300;301-600) 
		my ($gc,$atcg_len)=GC_calc($ref,$chr,$i,$win_width); # 20190905: The last win, get and use the actual length
		my $dpth_aver=Dpth_count($bam,$chr,$i,$atcg_len); # 20190905: The last win, use the actual length
		print "$chr\t$i\t$gc\t$dpth_aver\n";
	}
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Get the length of a Chr
sub Chr_len{
	my ($ref,$chr)=@_;
	my $seq=`samtools faidx $ref $chr`;
	my $seq_new=(split/\n/,$seq,2)[1];
	$seq_new=~ s/\n//g;
	my $len=length($seq_new);
	return $len;
}	

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calculate the GC content of a Win
sub GC_calc{
	my ($ref,$chr,$pos,$win_width)=@_;
	my $end=$pos+$win_width-1;
	my $seq=`samtools faidx $ref $chr:$pos-$end`;
	my $seq_new=(split/\n/,$seq,2)[1];
	my $num_g=($seq_new=~ s/G/G/g);
	my $num_c=($seq_new=~ s/C/C/g);
	my ($seq_len,$gc)=(0,0);
	$seq_new=~s/\n//g;
	$seq_len=length($seq_new);
	$gc=100*($num_g + $num_c)/$seq_len; # 20190905
	my $gc_new=sprintf("%.2f",$gc);
	return ($gc_new,$seq_len);
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Count the average depth of a Win (MQ >= 30)
sub Dpth_count{
	my ($bam,$chr,$pos,$seqseq_len)=@_;
	my $end=$pos+$seqseq_len-1;
	my $minMQ=30;
	my $cmd="samtools depth $bam -Q $minMQ -r $chr:$pos-$end";
	open DP,"$cmd |" || die $!;
	my $dp=0;
	while(<DP>){
		my @arr=split;
		$dp+=$arr[2];
	}
	close DP;
	my $dpth_aver=sprintf "%.2f",$dp/$seqseq_len; # 20190905
	return $dpth_aver;
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
