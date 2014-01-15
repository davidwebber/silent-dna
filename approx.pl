#! /usr/bin/perl -w

# $ sudo cpanm Text::LevenshteinXS
# $ sudo cpanm Text::WagnerFischer

#use Text::LevenshteinXS qw(distance);
use Text::WagnerFischer qw(distance);

 #print distance("foo","four")."\n";
 # prints "2"

 @array = ("foo","baro","boaz");

 my @distances = distance("foo",@array);
 #print "@distances"."\n";
 # prints "3"

 my $var = min_distance("AAAGAAATCCAAGAACAAAAACTC","AGAAATCCAAGACAAAAA");
 print "var is $var \n";
 #                "AGAAATCCAAG ACAAAAA")."\n";

 #print distance("TAAAGAAATCCAAGACAAAAACTC","AGAAATCCAAGACAAAAA")."\n";

 #print distance("AAAGAAATCCAAGACAAAAACTC","AGAAATCCAAGACAAAAA")."\n";
 #print distance("AAGAAATCCAAGACAAAAACTC","AGAAATCCAAGACAAAAA")."\n";
 #print distance("AGAAATCCAAGACAAAAACTC","AGAAATCCAAGACAAAAA")."\n";
 #print distance("GAAATCCAAGACAAAAACTC","AGAAATCCAAGACAAAAA")."\n";

sub min_distance{
    my ($seq, $ref) = @_;
    my $n=0;
    my $m=0;
    my $length = length $seq;
    my $dist_min=100;
 
    my $last_dist;
    my $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
    #print $dist."\n";

    # I think this algorithm could be improved
    for ($n=1; $n<10; $n++){
        $last_dist = $dist;
        $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
        print "$n $m $dist"."\n";
        if ($dist <= $last_dist){
        	$last_dist=$dist;
        } else {
        	$n-=1;
        	$dist = $last_dist;
        	last; #break;
        }
    }

    for ($m=1; $m<10; $m++){
        $last_dist = $dist;
        $dist = distance( (substr $seq, $n, $length-$n-$m), $ref);
        print "$n $m $dist"."\n";
        if ($dist <= $last_dist){
        	$last_dist=$dist;
        } else {
        	$m-=1;
        	$dist = $last_dist;
        	last; #break;
        }
    }

    print "minimum distance is $dist \n";
    return $dist;

}
