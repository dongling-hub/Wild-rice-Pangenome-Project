open OUTPUT,">$ARGV[0].len" or die "$!";

open INPUT,"<$ARGV[0]" or die "$!";
while (<INPUT>) {
        chomp;
        if (/^>(\S+)/) {
                if (exists $hash{$1}) {
                }
        }
        else {
                $hash{$1}.=$_;
        }
}

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

exit;
