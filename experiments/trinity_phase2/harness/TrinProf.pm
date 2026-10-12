package TrinProf;
# Opt-in profiling: set TRINITY_PROF_LOG=/path to append one line per
# child command (wall, child user, child sys) plus a TOTAL line at exit.
use strict;
use warnings;
use Time::HiRes ();
our $LOG = $ENV{TRINITY_PROF_LOG};
my $T0 = Time::HiRes::time();

sub timed {
    my ($label, $code) = @_;
    return $code->() unless $LOG;
    my @t0 = times();
    my $w0 = Time::HiRes::time();
    my $ret = $code->();
    my @t1 = times();
    my $w1 = Time::HiRes::time();
    _log(join("\t", "CMD", sprintf("%.4f", $w1 - $w0),
              sprintf("%.3f", $t1[2] - $t0[2]), sprintf("%.3f", $t1[3] - $t0[3]),
              _label($label)));
    return $ret;
}

sub _label {
    my $c = shift;
    $c =~ s/\s+/ /g;
    return substr($c, 0, 400);
}

sub _log {
    open(my $fh, ">>", $LOG) or return;
    print $fh "$$\t$_[0]\n";
    close $fh;
}

END {
    if ($LOG) {
        my @t = times();
        _log(join("\t", "TOTAL", sprintf("%.4f", Time::HiRes::time() - $T0),
                  map { sprintf("%.3f", $_) } @t));
    }
}
1;
