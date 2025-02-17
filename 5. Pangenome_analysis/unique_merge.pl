$identity=$ARGV[3];
open IN, "<$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].2022.$identity.merge";
open OUT2, ">$ARGV[0].2022.$identity.log";
open IN2, "<$ARGV[1]" or die "$!";
$n=0;
$aa=1;
while (<IN2>){
        $line = $_;
        chomp($line);
        $merge[0][$aa]=$line;
        $rrrhash{$line}=$aa;
        $aa++;
}
close IN2;

print $aa."\n";

open IN3, "<$ARGV[2]" or die "$!";
$n=0;
$aaa=0;
while (<IN3>){
        $line = $_;
        chomp($line);
        ($gname,$glen) = split(/\s+/,$line);
        $hash_len{$gname}=$glen;
        $aaa++;
}
close IN3;

print $aaa."\n";


#for (1..$aa-1){
#       print "\t".$merge[0][$_];
#}
#print "\n";

$n=0;
while (<IN>){
        $line = $_;
        chomp($line);
                $line = $_;
        chomp($line);
        if ($line eq ""){
                next;
        }
#       $line=~s/^\s+//;
        @hang = split(/\s+/,$line);
        $lane_n=$#hang;
#       for (0..($lane_n)){
        for (0..3){
                $array[$n][$_]=$hang[$_];
        }
#       $array[$n][20]=$line;
        $n++;
}
close IN;

RRRR: for (0..$n-1){
        $hh=$_;
        if ($array[$hh][1] eq "Query:"){
                $gene=$array[$hh][2];
                $gene_sn=$array[$hh][3];
#               $gene_exon=$array[$hh][4];
#               $gene_len=$array[$hh][10];
                $gene_len=$hash_len{$gene};
                print OUT2 "$gene = $gene_len\n";
#               print OUT "\n".$gene."\t".$gene_exon."\t".$gene_len;
#               @len=();
#               $m=0;
                $m++;
                $merge[$m][0]=$gene;
#               %rice_hash=();
                next RRRR;
        }
        if (($array[$hh][0] eq $gene) && ($array[$hh][1] ne $array[$hh-1][1])){
                $match=$array[$hh][2];
        }
        if (($array[$hh][0] eq $gene) && ($array[$hh][1] eq $array[$hh+1][1])){
                $subject=$array[$hh][1];
                $len=$len+$array[$hh][3];
                next RRRR;
        }
        if (($array[$hh][0] eq $gene) && ($array[$hh][1] ne $array[$hh+1][1])){
                $subject=$array[$hh][1];
                $len=$len+$array[$hh][3];

                ($rice,$scaffold) = split(/_/,$subject);
#               $pair1=$gene."___".$subject;
#               $pair2=$subject."___".$gene;
#####           if (exists $rice_hash{$rice}){
#####                   $len=0;
#####                   next RRRR;
#####           }
#####           if (exists $hash{$pair1}){
#####                   $len=0;
#####                   next RRRR;
#####           }
#####           if (exists $hash{$pair2}){
#####                   $len=0;
#####                   next RRRR;
#####           }
                if (($len>$gene_len/2) && ($match >=$identity)){
                        if ($merge[$m][$rrrhash{$rice}] eq ""){
                                $merge[$m][$rrrhash{$rice}]=$subject;
                        }
                        $len=0;
                        next RRRR;
                }
        }
}
print $m."\n";

for (1..$m-1){
        $bb=$_;
        for (1..$aa-1){
                $merge[$bb][$aa+1]=$merge[$bb][$aa+1].$merge[$bb][$_];
                if ($merge[$bb][$_] ne ""){
                        $count++;
                }
        }
        $merge[$bb][$aa]=$count;
        $count=0;
}


RRR1: for (0..$m-1){
        $cc=$_;
        if (exists $temp{$merge[$cc][$aa+1]}){
                $countcount++;
                next RRR1;
        }
        $temp{$merge[$cc][$aa+1]}++;
        print OUT $merge[$cc][0];
        for (1..$aa){
                print OUT "\t".$merge[$cc][$_];
        }
        print OUT "\n";
}
print $countcount." merged\n";