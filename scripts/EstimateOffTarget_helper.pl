BEGIN { $/ = "\n"; $\ = "\n"; }
LINE: while (defined($_ = <ARGV>)) {
    chomp $_;
    our(@F) = split(/\t/, $_, 0);
    use List::Util ('min', 'max');
    sub BEGIN {
        $CurrentPrimer = 'AG';
    }
    $F[0] =~ /(.+):(\w+)$/;
    if ($2 ne $CurrentPrimer) {
        $CurrentPrimer = $2;
        $GCCount = $2 =~ tr/CG//;
        $PrimerTm = 81.25 + 0.41 * $GCCount / length($CurrentPrimer) * 100 - 675 / length($CurrentPrimer);
        printf "%s\tSelf\t%s\t0\t0\t%s\t%s\n", $_, $PrimerTm, $2, $1;
    }
    else {
        $Candidate = substr($2, $F[6]);
        $GCCount = $Candidate =~ tr/CG//;
        $CandidateTm = 81.25 + 0.41 * $GCCount / $F[3] * 100 - 675 / $F[3] - (100 - $F[2]);
        printf "%s\tNonSelf\t%s\t%s\t%s\t%s\t%s\n", $_, $CandidateTm, $PrimerTm - $CandidateTm, max(0, 1 - ($PrimerTm - $CandidateTm) / 20), $2, $1;
    }
}
