#!/usr/bin/perl
use strict;
use Getopt::Long;
my $usage = <<USAGE;
Usage:
     $0 genome.fasta > telomere_info.txt

    大部分物种端粒序列的重复单元是TTAGGG/CCCTAA。本程序能在基因组中寻找端粒重复单元的串联重复序列，并给出位点信息。

    --split-length <int>    default: 100000
    --overlap-length <int>    default: 10000
    程序会将每条序列打断后进行重复单元搜索。这两个参数设置打断的序列长度和相邻两序列之间的重叠长度。

    --repeat-unit <string>    default: TTAGGG
    设置重复单元碱基序列，该重复单元的反向互补也将作为重复单元进行搜索。可以在端粒数据库（http://telomerase.asu.edu/sequences_telomere.html）中寻找目标端粒重复单元。
    vertebrate sp.      TTAGGG
    plants sp.          TTTAGGG
    Pezizomycotina      TTAGGG

    --min-repeat-num <int>    default: 4
    设置重复单元最小重复次数。
USAGE
if (@ARGV==0){die $usage}

my ($splitLength, $overlapLength, $repeatunit, $minRepeatNum);
GetOptions(
    "split-length:i" => \$splitLength,
        "overlap-length:i" => \$overlapLength,
            "repeat-unit:s" => \$repeatunit,
                "min-repeat-num:s" => \$minRepeatNum,
                 );
$splitLength ||= 100000;
$overlapLength ||= 10000;
$repeatunit ||= "TTAGGG";
$repeatunit = uc($repeatunit);
my $repeatunit_reverse = reverse $repeatunit;
$repeatunit_reverse =~ tr/ATCG/TAGC/;
$minRepeatNum ||= 4;
# 读取基因组序列
 open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
 my (%seq, $seq_id);
 while (<IN>) {
     chomp;
         if (m/^>(\S+)/) {
                 $seq_id = $1;
                     }
                         else {
                                 $_ = uc($_);
                                         $seq{$seq_id} .= $_;
                                             }
                                            }
                                             close IN;
# 将基因组序列打断
my (%seq_split, %seq_length);
foreach my $id (keys %seq) {
    my $seq = $seq{$id};
    my $length = length($seq);
    $seq_length{$id} = $length;
    my $pos = 0;
    while ($pos < $length) {
        $seq_split{$id}{$pos} = substr($seq, $pos, $splitLength + $overlapLength);
        $pos += $splitLength;
    }
}
print "SeqID\tSeqLength\tStart\tEnd\tLength\tType\n";
foreach my $id (sort keys %seq_split) {
    foreach my $pos (sort {$a <=> $b} keys %{$seq_split{$id}}) {
        my $seq = $seq_split{$id}{$pos};
        while ($seq =~ m/(($repeatunit){$minRepeatNum,})/g) {
            my $length = length($1);
            my $end = pos($seq);
            $end = $end + $pos;
            my $start = $end - $length + 1;
            print "$id\t$seq_length{$id}\t$start\t$end\t$length\t$repeatunit\n";
        }
        while ($seq =~ m/(($repeatunit_reverse){$minRepeatNum,})/g) {
            my $length = length($1);
            my $end = pos($seq);
            $end = $end + $pos;
            my $start = $end - $length + 1;
            print "$id\t$seq_length{$id}\t$start\t$end\t$length\t$repeatunit_reverse\n";
        }
    }
}
