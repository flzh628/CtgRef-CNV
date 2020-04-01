#!/usr/bin/perl -w

=for text
#########################################################################################################################################

=head1	
Script_Name:

=for text
	Name: slidWinDp2CNVs.pl
	Note: The core script of CtgRef-CNV pipeline to Call && Filter CNVs.

=head1	Description:

=for text
	1. Get the Wins of CNVs in each Chr ...
	2. Connect adjacent Wins into Intervals: DEL && DUP ...
	3. Extend candidate Intervals of DEL and DUP ...
	4. Begin to call and filter CNVs ...

=head1	Version:

=for text
	Author: Fengli Zhao	AGIS, CAAS	flzh628@gmail.com / zhaofengli@caas.cn
	Version: 13.cn_idx + Add the repeat seq of Nip to filter DUP
	Note: 2018/11/15: Add parameters of the bin_size and the min length of CNVs
	      2018/11/16: Add the parameter of ChrArr
	Date: 2018/11/16 2018/11/15 2018/4/13 2018/4/11 2018/4/10 2018/4/9 2018/3/30 2018/3/29 2017/12/25 2017/11/30 2017/11/20 2017/10/30 
	      2017/4/22 2017/4/21 2017/4/20

=head1	Usage:

slidWinDp2CNVs.pl [options]

	-in	slidWinDp-Res
	-bin	size of slid-Win (bp)
	-minlen	the min length of CNVs (bp)
	-dpcv	Aver-Depth and Coverage results of each Chr.
	-chrArr	chr name, splitted by commas 
	-bam	Bam(to Reference)
	-dir	CtgCNV-Dir
	-all	preCNV && Call && Filter the CNVs
	-preCNV	Wins && Extand
	-call	Call the CNVs
	-flt	Filter the CNVs
	-algn	nucmer-algn-flt
	-repeat	Position file of repeat seq in Nip
	-fltdir	CtgCNV-flt-Dir

=head1 Example:

=for text
	slidWinDp2CNVs.pl -in slidWinDp -bin 250 -minlen 1000 -dpcv DpthCvrg -dir CNV-Dir -bam Bam -all -algn nucmer-flt \ 
			      -fltdir CNV-flt-Dir -r Os7-Rep-Pos -chrArr Chr1,Chr2,Chr3 

=head1

=for text
#########################################################################################################################################

=cut

use strict;
use POSIX qw(strftime);
use Getopt::Long;

my ($slidwindp,$bin,$minlen,$dpcv,$chrArr,$refbam,$dir,$all,$preCNV,$call,$flt,$algnflt,$repeat,$fltdir,$help);
GetOptions(
	"in=s"=>\$slidwindp,
	"bin=s"=>\$bin,
	"minlen=s"=>\$minlen,
	"dpcv=s"=>\$dpcv,
	"chrArr=s"=>\$chrArr,
	"bam=s"=>\$refbam,
	"dir=s"=>\$dir,
	"all!"=>\$all,
	"preCNV!"=>\$preCNV,
	"call!"=>\$call,
	"flt!"=>\$flt,
	"algn=s"=>\$algnflt,
	"repeat=s"=>\$repeat,
	"fltdir=s"=>\$fltdir,
	"help|?|h"=>\$help,
);
die `pod2text $0` if ((defined $all) && (!defined $slidwindp || !defined $bin || !defined $minlen || !defined $dpcv || !defined $refbam || !defined $repeat));
die `pod2text $0` if ((defined $all) && (!defined $dir || !defined $fltdir || !defined $dpcv));
die `pod2text $0` if ((defined $preCNV) && (!defined $slidwindp || !defined $dir));
die `pod2text $0` if ((defined $call) && (!defined $slidwindp || !defined $bin || !defined $minlen || !defined $dpcv || !defined $refbam || !defined $dir || defined $preCNV || !defined $repeat));
die `pod2text $0` if ((defined $flt) && (!defined $fltdir || !defined $dir || defined $preCNV));
die `pod2text $0` if (defined $help || (!defined $all && !defined $preCNV && !defined $call && !defined $flt));

my $name=$1 if $slidwindp=~/.*(X\d+)-nip-slidWinDp\d+.*/;
if (defined $bin){
	die "\n\tError: Please make sure the ratio of minlen/bin is a nonzero integer !\n\n" if ($minlen % $bin);
}

############################################ Preparatory Work ############################################

my (%DpWin,%AverDpth,%AvDpFlt,$newalgnflt,$pid,$vsz,$rss,$mem);
if (!-e $dir){
        mkdir $dir,0755;#CtgCNV-Res-Dir
}

my (@ChrArr,$newchr);
if (defined $chrArr){
	@ChrArr=split/,/,$chrArr;
	$newchr=$chrArr;
	$newchr=~s/,/-/g;
}

if(($all) || ($preCNV) || ($call)){
	my $t1=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	if (defined $chrArr){
		print STDERR "=" x 40,"[Begin to process $name-$newchr data]","=" x 40,"\n";
	}else{
		print STDERR "=" x 40,"[Begin to process $name data]","=" x 40,"\n";
	}
	print STDERR "[$t1] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tReading the AverDp of each Chr and Dpth of each Win, building the Res and Flt DIR ...\n";
	
	open DP,"$slidwindp" or die $!;
	while(<DP>){
		chomp;
		my @arrwindp=split;
		$DpWin{$arrwindp[0]}{$arrwindp[1]}=$arrwindp[2];
	}
	close DP;
	
	open AA,"$slidwindp-ChrDp" or die $!;
	while(<AA>){
		chomp;
		my ($chrid,$averDp)=split;
		$AverDpth{$chrid}=$averDp;
	}
	close AA;
}

if(($all) || ($call)){
	open BB,"tail -n +2 $dpcv | " or die $!; # ChrNum	AverDp	ChrCoverage	
	while(<BB>){				 # 10	48.31	0.9354
        	chomp;
        	my ($inf,$averDp)=(split)[0,1];
        	my $chrid="Chr$inf";
        	$AvDpFlt{$chrid}=$averDp;
	}
	close BB;
	
	if (! -e "$refbam.bai"){
        	my $cmd1="samtools index $refbam &";
        	system($cmd1);
	}
}

if(($all) || ($flt)){
	if (!-e $fltdir){
		mkdir $fltdir,0755;#CtgCNV-Res-flt-Dir
	}
	$newalgnflt=$1 if $algnflt=~/.*(X\d+-Algn-Flt).*/;
	if (defined $chrArr){
		system("grep '\\b$ChrArr[0]\\b' $algnflt | awk '{print \$4\"\t\"\$5\"\t\"\$6}' >$dir/$newalgnflt-Ref-$newchr");
		if ($#ChrArr>=1){
			for my $i (1 .. $#ChrArr){
				system("grep '\\b$ChrArr[$i]\\b' $algnflt | awk '{print \$4\"\t\"\$5\"\t\"\$6}' >>$dir/$newalgnflt-Ref-$newchr");
			}
		}
		system("msort -k m1,n2,n3 $dir/$newalgnflt-Ref-$newchr >$dir/$newalgnflt-Ref-$newchr.msort");
	}else{
		system("awk '{print \$4\"\t\"\$5\"\t\"\$6}' $algnflt >$dir/$newalgnflt-Ref");
		system("msort -k m1,n2,n3 $dir/$newalgnflt-Ref >$dir/$newalgnflt-Ref.msort");
	}
}

############################################# Detect Wins of CNVs #############################################

if (($all) || ($preCNV)){
	my $t2=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	print STDERR "[$t2] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tBegin to detect candidate Intervals of CNVs ...\n";

	&Del_Dup_Win;#Detect candidate Intervals
	&Del_Dup_Merge;#Extend candidate Intervals
}

################################################## Call CNVs ##################################################

if (($all) || ($call)){
	my $t3=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	print STDERR "[$t3] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tBegin to call CNVs ...\n";

	my $get_lines;
	if (defined $chrArr){
		open WIN,"$dir/$name-DelDup-Merge-$newchr" or die $!;
		open RES,">$dir/$name-DelDup-Result-$newchr" or die $!;
		$get_lines=`wc -l $dir/$name-DelDup-Merge-$newchr`;
	}else{
		open WIN,"$dir/$name-DelDup-Merge" or die $!;
		open RES,">$dir/$name-DelDup-Result" or die $!;
		$get_lines=`wc -l $dir/$name-DelDup-Merge`;
	}
	my $total_lines=$1 if $get_lines=~/(\d+)\s\S+\n/;
	my ($now_lines,$perc_id)=(0,0); my @perc;
	for my $m(1 .. 9){
		my $linenum=int($total_lines*$m*0.1);
		push @perc,$linenum;
	}
	push @perc,$total_lines;
	print RES "CNV\tCHR\tBegin\tEnd\tLength\tCN_index\tDp_bam\tCvrg_%\tDp10_%\t1.5Dp%\tAlgn_Cvrg\tAlgn_max_vcrg\tRep_Ovlp_Rate\n";
	print STDERR "\t";
	my $min_num=$minlen/$bin;
	while(<WIN>){
		chomp;
		my ($typechr,$win1,$win2)=split;
		my ($type,$chr)=($1,$2) if $typechr=~/(\w+)-(\w+)/;
		$now_lines++;
		next if ($win2-$win1+1 < $min_num); # length
		my $aver_dp=Aver_Dp_win($chr,$win1,$win2);
		my $beg=$bin*$win1;
		my $end=$bin*($win2+1);
		my $len=$end-$beg;
		my $cn_index=sprintf "%.2f",$aver_dp/$AverDpth{$chr};##20180117；
		my ($depth,$cvrg,$dp10rate,$dupcvrg)=AverDp_Cvrg_bam($refbam,$chr,$beg,$end,$AvDpFlt{$chr});
		my ($algncvrg,$algncvrg_max)=Get_Algn_Coverage($chr,$beg,$end);
		my $rep_rate=Get_Rep_Ovlp_rate($chr,$beg,$end);
		print RES "$type\t$chr\t$beg\t$end\t$len\t$cn_index\t$depth\t$cvrg\t$dp10rate\t$dupcvrg\t$algncvrg\t$algncvrg_max\t$rep_rate\n";
		if ($perc_id<=$#perc){
			if ($now_lines>=$perc[$perc_id]){
				my $now_time=strftime "%H:%M:%S",localtime;
				my $now_percent=($perc_id+1)*10;
				print STDERR "\t[$now_time]: $now_percent%";
				$perc_id++;
			}
		}
	}
	print STDERR "\n";
	close WIN;
	close RES;

	my $t4=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	if (defined $chrArr){
		print STDERR "[$t4] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tOK! The CNV calling of $name-$newchr is done! \n";
	}else{
		print STDERR "[$t4] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tOK! The CNV calling of $name is done! \n";
	}
}

####################################################### Filter CNVs ####################################################

if (($all) || ($flt)){
	my $t5=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	print STDERR "[$t5] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tBegin to filter CNVs ...\n";

	if (defined $chrArr){
		open RES,"$dir/$name-DelDup-Result-$newchr" or die $!;
		open FLT,">$fltdir/$name-DelDup-Res-Flt-$newchr" or die $!;
	}else{
		open RES,"$dir/$name-DelDup-Result" or die $!;
		open FLT,">$fltdir/$name-DelDup-Res-Flt" or die $!;
	}
	while(<RES>){
		chomp;
		if ($.==1){
			print FLT "$_\n";
		}else{
			my ($type,$chr,$beg,$end,$len,$cn_index,$depth,$cvrg,$dp10rate,$dupcvrg,$algncvrg,$algncvrg_max,$rep_rate)=split;
			if ($type eq 'DEL'){
				next if $cvrg>=50; next if $dp10rate>=40; next if $cn_index>=0.4; next if $algncvrg>=50; next if $algncvrg_max>=40; next if $rep_rate>=50;
				print FLT "$type\t$chr\t$beg\t$end\t$len\t$cn_index\t$depth\t$cvrg\t$dp10rate\t$dupcvrg\t$algncvrg\t$algncvrg_max\t$rep_rate\n";
			}else{
				next if $cn_index<1.7;  next if $algncvrg<90; next if $rep_rate>=50;
				next if ($algncvrg_max<30 && $len>=10000); next if ($algncvrg_max<80 && $len<10000); #20180418
				print FLT "$type\t$chr\t$beg\t$end\t$len\t$cn_index\t$depth\t$cvrg\t$dp10rate\t$dupcvrg\t$algncvrg\t$algncvrg_max\t$rep_rate\n";
			}
		}
	}
	close RES;
	close FLT;

	my $t6=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	print STDERR "[$t6] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tThe CNVs have been filtered ...\n";
}

my $t7=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t7] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tCongratulations! The analysis is finished ...\n";
if (defined $chrArr){
	print STDERR "=" x 35,"[The CNV calling if $name-$newchr is finished!]","=" x 35,"\n";
}else{
	print STDERR "=" x 35,"[The CNV calling if $name is finished!]","=" x 35,"\n";
}

############################################## END of main program ###############################################

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Del_Dup_Win{	
	my $t01=strftime "%Y-%m-%d %H:%M:%S",localtime;
        ($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
        print STDERR "\t[$t01] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tGet the Wins of CNVs in each Chr ...\n";
	
	if (defined $chrArr){
		system("grep '\\b$ChrArr[0]\\b' $slidwindp >$dir/$name-nip-slidDp-$newchr");
		if ($#ChrArr>=1){
			for my $i (1 .. $#ChrArr){
				system("grep '\\b$ChrArr[$i]\\b' $slidwindp >>$dir/$name-nip-slidDp-$newchr");
			}
		}
		open IN,"$dir/$name-nip-slidDp-$newchr" or die $!;
		open DOTS,">$dir/$name-DelDup-dot-$newchr" or die $!;
	}else{
		open IN,"grep -v 'CHR' $slidwindp |" or die $!;
		open DOTS,">$dir/$name-DelDup-dot" or die $!;
	}
	while(<IN>){
		chomp;
		my @e=split;
		if ($e[2] <0.5*$AverDpth{$e[0]}){#For Deletion
			print DOTS "DEL-$e[0]\t$e[1]\n";
		}elsif($e[2] >=1.5*$AverDpth{$e[0]}){#For Duplication
			print DOTS "DUP-$e[0]\t$e[1]\n";
		}
	}
	close IN;
	close DOTS;
	
	if (defined $chrArr){
		system("msort -k m1,n2 $dir/$name-DelDup-dot-$newchr >$dir/$name-DelDup-dot-$newchr.msort");
	}else{
		system("msort -k m1,n2 $dir/$name-DelDup-dot >$dir/$name-DelDup-dot.msort");
	}

	my $t02=strftime "%Y-%m-%d %H:%M:%S",localtime;
        ($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
        print STDERR "\t[$t02] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tConnect adjacent Wins into Intervals: DEL && DUP ...\n";
	if (defined $chrArr){
		open DOT,"$dir/$name-DelDup-dot-$newchr.msort" or die $!;
		open WIN,">$dir/$name-DelDup-Win-$newchr" or die $!;
	}else{
		open DOT,"$dir/$name-DelDup-dot.msort" or die $!;
		open WIN,">$dir/$name-DelDup-Win" or die $!;
	}
	my ($type,$dot);
	while(<DOT>){
		chomp;
		my @wins=split;
		if ($.==1){
			($type,$dot)=($wins[0],$wins[1]);
			print WIN "$type\t$dot";
			next;
		}else{
			if ($wins[0] eq $type){
				if ($wins[1] == $dot+1){
					$dot=$wins[1];
					if (eof(DOT)){
						print WIN "-$dot\n";
					}
					next;
				}else{
					print WIN "-$dot\t$wins[1]";
					$dot=$wins[1];
				}
			}else{
				print WIN "-$dot\n";
				($type,$dot)=($wins[0],$wins[1]);
				print WIN "$type\t$dot";
			}
		}
	}
	close DOT;
	close WIN;
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Del_Dup_Merge{
	my $t01=strftime "%Y-%m-%d %H:%M:%S",localtime;
	($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
	print STDERR "[$t01] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tExtending candidate Intervals of DEL and DUP ...\n";
	
	if (defined $chrArr){
		open WIN,"$dir/$name-DelDup-Win-$newchr" or die $!;
		open EXD,">$dir/$name-DelDup-Extend-$newchr" or die $!;
	}else{
		open WIN,"$dir/$name-DelDup-Win" or die $!;
		open EXD,">$dir/$name-DelDup-Extend" or die $!;
	}
	while(<WIN>){
		chomp;
		my ($type,$internal)=(split/\t/,$_,2);
		my @win=split/\t/,$internal;
		my ($typenew,$chr)=($1,$2) if $type=~/(\w+)-(\w+)/;
		for my $i(0 .. $#win){
			my ($left,$right);
			($left,$right)=($1,$2) if $win[$i]=~/^(\d+)-(\d+)$/;
			($left,$right)=($1,$1) if $win[$i]=~/^(\d+)$/;
			while(1){
				last if $left==0;
				my $leftnew=$left-1;
				my $diff=abs($DpWin{$chr}{$leftnew}-$DpWin{$chr}{$left});
				last if $diff>0.25*$AverDpth{$chr};
				my $newcnvdp=Aver_Dp_win($chr,$leftnew,$right);
				last if (($typenew eq 'DEL') && ($newcnvdp>0.5*$AverDpth{$chr}));
				last if (($typenew eq 'DUP') && ($newcnvdp<1.5*$AverDpth{$chr}));
				$left=$leftnew;
			}
			while(1){
				my $rightnew=$right+1;
				last if (! defined $DpWin{$chr}{$rightnew});
				my $diff=abs($DpWin{$chr}{$rightnew}-$DpWin{$chr}{$right});
				last if $diff>0.25*$AverDpth{$chr};
				my $newcnvdp=Aver_Dp_win($chr,$left,$rightnew);
				last if (($typenew eq 'DEL') && ($newcnvdp>0.5*$AverDpth{$chr}));
				last if (($typenew eq 'DUP') && ($newcnvdp<1.5*$AverDpth{$chr}));
				$right=$rightnew;
			}
			if ($left!=$right){
				print EXD "$type\t$left\t$right\n";
			}
		}
	}
	close WIN;
	close EXD;
	
	if (defined $chrArr){
		system("msort -k m1,n2,n3 $dir/$name-DelDup-Extend-$newchr >$dir/$name-DelDup-Extend-$newchr.msort");
		open EXD,"$dir/$name-DelDup-Extend-$newchr.msort" or die $!;
		open MRG,">$dir/$name-DelDup-Merge-$newchr" or die $!;
	}else{
		system("msort -k m1,n2,n3 $dir/$name-DelDup-Extend >$dir/$name-DelDup-Extend.msort");
		open EXD,"$dir/$name-DelDup-Extend.msort" or die $!;
		open MRG,">$dir/$name-DelDup-Merge" or die $!;
	}
	my ($type,$win1,$win2);##20171225: more faster
	while(<EXD>){
		chomp;
		my @extend=split;
		if ($.==1){
			($type,$win1,$win2)=($extend[0],$extend[1],$extend[2]);
			next;
		}else{
			if ($extend[0] eq $type){
				if ($extend[1]<=$win2){
					my @mn=sort {$a <=> $b} ($extend[2],$win2);
					$win2=$mn[1];
				}elsif($extend[1]==$win2+1){ #20180419  deBug
					$win2=$extend[2];
				}else{
					print MRG "$type\t$win1\t$win2\n";
					($win1,$win2)=($extend[1],$extend[2]);
				}
			}else{
				print MRG "$type\t$win1\t$win2\n";
				($type,$win1,$win2)=($extend[0],$extend[1],$extend[2]);
			}
		}
	}
	close EXD;
	close MRG;
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Aver_Dp_win{
	my ($chr,$w1,$w2)=@_;
	my ($dp,$averdp)=(0,0);
	for my $wins(sort {$a <=> $b} keys %{$DpWin{$chr}}){ # V9: 20180330
		next if $wins<$w1;
		last if $wins>$w2;
		$dp+=$DpWin{$chr}{$wins};
	}
	$averdp=sprintf "%.2f",$dp/($w2-$w1+1);
	return $averdp;
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub AverDp_Cvrg_bam{
        my ($bam,$chr,$beg,$end,$AverDp)=@_;
        my ($dp,$num1,$num2,$num3)=(0,0,0,0);
	my $mpQ=30;
        my $cmd="samtools depth $bam -Q $mpQ -r $chr:$beg-$end";
        open DP,"$cmd |" or die $!;
        while(<DP>){
                chomp;
                my @ee=split;
                $dp+=$ee[2];
		if ($ee[2]>=1){
			$num1++;
			if ($ee[2]>=10){
				$num2++;
				if ($ee[2]>=1.6*$AverDp){
					$num3++;
				}
			}
		}	
        }
        close DP;
        my $aver_dp=sprintf "%.2f",$dp/($end-$beg+1);
	my $coverage=sprintf "%.2f",100*$num1/($end-$beg+1);
	my $dp10_rate=sprintf "%.2f",100*$num2/($end-$beg+1);
	my $dupdp_rate=sprintf "%.2f",100*$num3/($end-$beg+1);
        return ($aver_dp,$coverage,$dp10_rate,$dupdp_rate);
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Get_Algn_Coverage{
	my ($chr,$beg,$end)=@_;
	my ($coverage,$maxcvrg)=(0,0);
	my @cvrg=(); my @cvrglen=('0'); 
	if (defined $chrArr){
		open AA,"grep '\\b$chr\\b' $dir/$newalgnflt-Ref-$newchr.msort |" or die $!;
	}else{
		open AA,"grep '\\b$chr\\b' $dir/$newalgnflt-Ref.msort |" or die $!;
	}
	my @match=<AA>;               
	foreach my $line(@match){     
		my ($left,$right)=(split/\t/,$line)[1,2];
		next if $right<$beg;
		last if $left>$end;
		my @digit=sort {$a <=> $b} ($beg,$end,$left,$right);
		my $overlap=$digit[2]-$digit[1]+1;
		push @cvrglen,$overlap;
		for my $numid($digit[1] .. $digit[2]){
			push @cvrg,$numid;
		}
	}
	close AA;
	my %count1=();
	@cvrg=grep { ++$count1{$_} < 2 } @cvrg;
	$coverage=sprintf "%.2f",100*($#cvrg+1)/($end-$beg+1);
	my @ovlp=sort {$b <=> $a} @cvrglen;
	$maxcvrg=sprintf "%.2f",100*$ovlp[0]/($end-$beg+1);
	return ($coverage,$maxcvrg);
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Get_Rep_Ovlp_rate{
	my ($chr,$beg,$end)=@_;
	my $rate=0;
	my @repcvrg=();
	open REP,"grep '\\b$chr\\b' $repeat |" or die $!;
	my @othmat=<REP>;
	foreach my $ovlpline(@othmat){
		my ($repleft,$repright)=(split/\t/,$ovlpline)[2,3];
		next if $repright<$beg;
		last if $repleft>$end;
		my @repdigit=sort {$a <=> $b} ($beg,$end,$repleft,$repright);
		for my $numid($repdigit[1] .. $repdigit[2]){
			push @repcvrg,$numid;
		}
	}
	close REP;
	my %count2=();
	@repcvrg=grep { ++$count2{$_} < 2 } @repcvrg;
	$rate=sprintf "%.2f",100*($#repcvrg+1)/($end-$beg+1);
	return $rate;
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Get_RSS_MEM{
        my ($pid,$vsz,$rss,$mem)=(0,0,0,0);
	if (defined $chrArr){
		$pid=`ps -ef|grep '$0' |grep '$name-nip-slidWinDp' |grep '\\b$ChrArr[0]\\b'|grep -v grep |awk '{print \$2}'`;
	}else{
        	$pid=`ps -ef|grep '$0' |grep '$name-nip-slidWinDp'|grep -v grep |awk '{print \$2}'`;
	}
        $pid=~s/\n//g;
        my $inf=`pidstat -r -p $pid`;
        my @new=split/\n/,$inf;
        my @newinf=split/\s+/,$new[3];
	if ($new[3] =~ /(A|P)M/){
		$mem=$newinf[7];
		$vsz=sprintf "%.2f",$newinf[5]/(1028*1028);
		$rss=sprintf "%.2f",$newinf[6]/(1028*1028);
	}else{
        	$mem=$newinf[6];
        	$vsz=sprintf "%.2f",$newinf[4]/(1028*1028);
        	$rss=sprintf "%.2f",$newinf[5]/(1028*1028);
	}
        return($pid,$vsz,$rss,$mem);
}
