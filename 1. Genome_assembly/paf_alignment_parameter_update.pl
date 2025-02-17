#calculate identity
#
use strict;
open INPUT1,"<",$ARGV[0] or die "Can't open INPUT1:$!\n";
my @array1;
my $identity;
my $ratio;
my $id;
my $match;
my $align;
my $read_len;
my $alignlength;
my $lastid;
my %hash1;
open OUT,">",$ARGV[0].".aln" or die "Can't open OUT:$!\n";
print OUT "Read_ID\tRead_Length\tAlign_Length\tAlign_Identity\tAlign_Ratio\n";
NNNN:while(<INPUT1>){
        chomp;
        @array1=split/\t/,$_;
        if($hash1{$array1[0]}){
              #if(($array1[2]==$lastcost)&&($array1[3]==$lastcoed)){
             next NNNN;
        #       }else{
              #   $match += $array1[9];
        #   $alignlength=$array1[3]-$array1[2]+1;
              #   $align += $alignlength;
              # }
        }else{
                $hash1{$array1[0]}=1;
                $align=$array1[3]-$array1[2]+1;
                $match=$array1[9];
                $id=$array1[0];
                $read_len=$array1[1];
                if(($align != 0) && ($read_len != 0)){                           #?????
                  $identity=$match/$align;
      $ratio=$align/$read_len;
      print OUT $id."\t".$read_len."\t".$align."\t".$identity."\t".$ratio."\n"; 
    }

  }
  $lastid=$array1[0];
}
close INPUT1;
close OUT;
