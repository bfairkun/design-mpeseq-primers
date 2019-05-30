#!/usr/bin/env perl
BEGIN { $/ = "\n"; $\ = "\n"; }
LINE: while (defined($_ = <ARGV>)) {
    chomp $_;
    our(@F) = split(/\t/, $_, 0);
    use List::Util ('min', 'max');
    sub BEGIN {

    }
    {
        if ($F[5] eq '+') {
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $F[0], $F[1], min($F[2], $F[8]), $F[3], $F[4], $F[5];
        }
        elsif ($F[5] eq '-') {
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $F[0], max($F[1], $F[7]), $F[2], $F[3], $F[4], $F[5];
        }
    }
}
