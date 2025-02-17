open OUTPUT,">$ARGV[0].len" or die "$!";
#### open OUTPUT2,">$ARGV[0].error" or die "$!";

open INPUT,"<$ARGV[0]" or die "$!";
while (<INPUT>) {
        chomp;
        if (/^>(\S+)/) {
                if (exists $hash{$1}) {
####                    print OUTPUT2 "$1\n";
                }
        }
        else {
                $hash{$1}.=$_;
        }
}

#foreach (sort keys %hash) {
        #print OUTPUT "$_ =>\n$hash{$_}\n\n";
#       print OUTPUT "$_\t",length $hash{$_},"\n";
#}

close INPUT;

$n=0;
open INPUT,"<$ARGV[0]" or die "$!";
while (<INPUT>) {
        chomp;
        if (/^>(\S+)/) {
                print OUTPUT "$1\t",length $hash{$1},"\n";
                $len[$n]=length $hash{$1};
                $n++;
        }
}
close INPUT;

for (0..$n-1){
        $total_len=$len[$_]+$total_len;
}
print "Total length of the sequences in $ARGV[0] is: $total_len bp.\n";
print "Total $n sequences\n";

@len = reverse sort { $a <=> $b } @len;

#for (0..$n-1){
#       print "$len[$_] bp\n";
#}

for (0..$n-1){
        $total_len_n50=$len[$_]+$total_len_n50;
        if ($total_len_n50 >= $total_len/2){
                print "Max ctg: $len[0] bp\n";
                print "Minimal ctg: $len[$n-1] bp\n";
                $average=int(($total_len/$n)+0.5);
                print "Average ctg: $average bp\n";
                print "N50: $len[$_] bp\n";
                last;
        }
}
exit;
