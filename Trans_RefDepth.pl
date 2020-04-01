	#!/usr/bin/perl -w

use POSIX ();
use POSIX qw(strftime);

die "\n\tAuthor: Fengli Zhao  AGI,CAAS	flzh628\@gmail.com
	Version: 6.0	20180415  20180103  20171216  20171213 ...
	It is designed to filter the aln-res of nucmer, and use it and Ctg-slidWinDp 
	to convert the depth information from Ctg to Nip.
	Usage: $0 nucmer-coords slidWin-ctg-dp-cor NipDpgz-DIR Temp-DIR \n\n" unless @ARGV==4;

if(!-e $ARGV[2]){ #Res-DIR
	mkdir $ARGV[2],0755;
}
if(!-e $ARGV[3]){ #Temp-DIR
	mkdir $ARGV[3],0755;
}

my $id=$1 if $ARGV[0]=~/.*(X\d+)-nip.coords/;
#Filter and Uniform the Nucmer result
my $t1=strftime "%Y-%m-%d %H:%M:%S",localtime;
my ($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "=" x 56,"[Begin to process $id data]","=" x 56,"\n";
print STDERR "[$t1] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tUniform and Filter the nucmer result ....\n";
open IN,"$ARGV[0]" or die $!;
open TM,">$ARGV[3]/$id-Algn" or die $!;
while(<IN>){
        chomp;
        next unless /Chr/;
        s/_c//g;s/-RC//g;s/^\s+//g;
        s/\|//g;
        my ($ctg,$begctg,$endctg,$len_nip,$len_ctg,$idty,$chr,$beg,$end)=(split)[-1,2,3,4,5,6,7,0,1];
        next if ($idty<85); next if ($len_nip<100); next if ($len_ctg<100);
        my ($left,$right)=sort {$a <=> $b} ($begctg,$endctg);
        print TM "$ctg\t$left\t$right\t$chr\t$beg\t$end\n";
}
close IN;
close TM;

my $cmd1="msort -k m1,n2,rn3 $ARGV[3]/$id-Algn >$ARGV[3]/$id-Algn.msort";
system("$cmd1");

open MS,"$ARGV[3]/$id-Algn.msort" or die $!;
open OUT,">$ARGV[3]/$id-Algn-Flt" or die $!;
my ($ctg,$begc,$endc,$chr,$beg,$end);
while(<MS>){
        chomp;
        my @arr=split;
        if ($.==1){
                ($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
                next;
        }else{
                if ($arr[0] eq $ctg){
                        my $len1=$arr[2]-$arr[1]+1;
                        my $len2=$endc-$begc+1;
                        my ($short,$long)=sort {$a <=> $b} ($len1,$len2);
                        if ($arr[2] <= $endc){
                                if ($short/$long >=0.8){
                                        print OUT "$ctg\t$begc\t$endc\t$chr\t$beg\t$end\n";
                                        ($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
                                }else{
                                        next;
                                }
                        }elsif($arr[1]<$endc){
                                my @digit=sort {$a <=> $b} ($begc,$endc,$arr[1],$arr[2]);
                                my $rate=($digit[2]-$digit[1]+1)/$long;
                                if ($rate>=0.7){
                                        print OUT "$ctg\t$begc\t$endc\t$chr\t$beg\t$end\n";
                                        ($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
                                }else{
                                        if ($len1/$len2 <0.2){
                                                next;
                                        }elsif($len2/$len1 <0.2){
                                                ($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
                                        }else{
                                                print OUT "$ctg\t$begc\t$endc\t$chr\t$beg\t$end\n";
                                                ($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
					}
				}
			}else{
				print OUT "$ctg\t$begc\t$endc\t$chr\t$beg\t$end\n";
				($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
			}
		}else{
			print OUT "$ctg\t$begc\t$endc\t$chr\t$beg\t$end\n";
			($ctg,$begc,$endc,$chr,$beg,$end)=($arr[0],$arr[1],$arr[2],$arr[3],$arr[4],$arr[5]);
		}
	}
}
close MS;
close OUT;

my (%DotCtg,%DPctgB,%NumCtgB,%DPctgs,%NumCtgs,%DpNip);# %NumCtgB: Big Ctg win;20171216
my $t2=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t2] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tRead the algn-flt and store dots in the Ctg ....\n";
open AC,"$ARGV[3]/$id-Algn-Flt" or die $!;
while(<AC>){
	chomp;
	my ($ctg,$begctg,$endctg)=(split)[0,1,2];
	push @{$DotCtg{$ctg}},$begctg; push @{$DotCtg{$ctg}},$endctg;
	my $winctg="$begctg-$endctg";
	$NumCtgB{$ctg}{$winctg}++;
}
close AC;

my $t3=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t3] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tRead the Ctg-slidWinDp file, update its dots and store Dpth ....\n";
open DP,"$ARGV[1]" or die $!;
while(<DP>){
        chomp;
        next if $.==1;
        next if /Aver/; next if /=/; next if /^$/; next if /GC/; next if /^\d+/;#20171026
        my ($ctg,$win,$dpcor)=(split)[0,1,4];
        my $end=$win+299;
        push @{$DotCtg{$ctg}},$win; push @{$DotCtg{$ctg}},$end;
        $DPctgB{$ctg}{$win}=$dpcor;
}
close DP;

my $t4=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t4] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tAllot the Number and Dpth for small Wins of Ctg ....\n";
for my $ctg(keys %DotCtg){
        my %count=();
        my @new = grep {++$count{$_}<2;} @{$DotCtg{$ctg}};
        my @newarr= sort {$a <=> $b} @new;
        my $id=$#newarr-1;
        for my $i(0 .. $id){
                my $pos="$newarr[$i]-$newarr[$i+1]";
                $NumCtgs{$ctg}{$pos}=0;
                for my $win (keys %{$DPctgB{$ctg}}){
                        next if ($newarr[$i]>$win+299); next if $newarr[$i+1]<$win;#20171216
                        if (($newarr[$i] >= $win) && ($newarr[$i+1] <= $win+299)){
                                $DpCtgs{$ctg}{$pos}=$DPctgB{$ctg}{$win};
                        }
                }
		for my $winn(keys %{$NumCtgB{$ctg}}){
			my ($zuo,$you)=($1,$2) if $winn=~/(\d+)-(\d+)/;
			next if $zuo>$newarr[$i+1]; next if $you<$newarr[$i];
			if (($newarr[$i] >= $zuo) && ($newarr[$i+1] <= $you)){
				$NumCtgs{$ctg}{$pos}+=$NumCtgB{$ctg}{$winn};
			}
		}
        }
}

my $t5=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t5] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tRead the algn-flt again and convert depth from Ctg to Nip ....\n";
open BC,"$ARGV[3]/$id-Algn-Flt" or die $!;
while(<BC>){
        chomp;
        my ($chr,$begnip,$endnip,$ctg,$begctg,$endctg)=(split)[3,4,5,0,1,2];
	my $winnip="$begnip-$endnip";
	my $dpnew=0;
	for my $win(keys %{$DpCtgs{$ctg}}){#######20171218 $NumCtgs => $DpCtgs
		my ($beg,$end)=($1,$2) if $win=~/(\d+)-(\d+)/;
		next if $beg>$endctg; next if $end<$begctg;#20171216
		my $rate=($end-$beg+1)/($endnip-$begnip+1);###20171218 
		if (($beg >= $begctg) && ($end <= $endctg)){
			if (defined $NumCtgs{$ctg}{$win}){#######20171218 $DpCtgs => $NumCtgs
				$dpnew+=sprintf "%.2f",$rate*$DpCtgs{$ctg}{$win}/$NumCtgs{$ctg}{$win};
			}else{
				$dpnew+=0;
			}
		}
	}
	$DpNip{$chr}{$winnip}=$dpnew;
}
close BC;

my $t6=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t6] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tIntegrate and export the Dpth of uniq sites in Nip ....\n";
my %NewDp;
for my $chr(keys %DpNip){
	for my $win(keys %{$DpNip{$chr}}){
		my ($begn,$endn)=($1,$2) if $win=~/(\d+)-(\d+)/;
		for my $m($begn .. $endn){
			$NewDp{$chr}{$m}+=$DpNip{$chr}{$win};
		}
	}
}
open OUT,"|gzip -c - >$ARGV[2]/$id-Dpth.gz" or die $!;
for $chr(sort keys %NewDp){
	for my $dot (keys %{$NewDp{$chr}}){
		print OUT "$chr\t$dot\t$NewDp{$chr}{$dot}\n";
	}
}
close OUT;

my $t7=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t7] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tExport the results successfully !\n";
my $t8=strftime "%Y-%m-%d %H:%M:%S",localtime;
($pid,$vsz,$rss,$mem)=Get_RSS_MEM();
print STDERR "[$t8] PID:$pid:VSZ:$vsz:RSS:$rss:%MEM:$mem\tCongratulations! The step has finished! Please run the next! \n";
print STDERR "=" x 61,"[$id-is-Finished]","=" x 61,"\n";
POSIX::_exit(0);

#my $cmd="less $dpout |msort -k f1,n2 - |gzip -c -9 - >./$dpout-msort";
#system("$cmd");
#unlink "$dpout";
#rename "$dpout-msort","$dpout";

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub Get_RSS_MEM{
	my ($pid,$vsz,$rss,$mem)=(0,0,0,0);
	my $scriptname=$1 if $0=~/(.*).pl/;
	$pid=`ps -ef|grep '$scriptname' |grep '$id'|grep -v grep |awk '{print \$2}'`;
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
