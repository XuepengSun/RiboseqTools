package SXP::Util;
use strict;
use warnings;
use Term::ANSIColor;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(LOG);
 
sub LOG {
	$| = 0;
    my $print = shift;
    my $time = localtime;
    print STDERR color("green"),'[',$time,'] ',color("red"),$print,"\n",color("reset");
}
 
1;

