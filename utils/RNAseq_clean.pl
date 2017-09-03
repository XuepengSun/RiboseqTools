#!/usr/bin/perl -w
use File::Find;
use threads;
use Cwd;

unless(@ARGV == 3){print STDERR "\n\nUseage: perl $0 [rawReads dir] [output dir] [rRNA database index]\n\n";exit}

my ($source,$target,$rRNA) = @ARGV;
my (%data);

my $threads_trim = 20;
my $threads_map = 60;

my $tmpdir = $target.'/tmp';
unless (-d $tmpdir){mkdir $tmpdir}

my @files = glob("$source/*.gz");
$source=~s/\/$//;
$target=~s/\/$//;

while(1){
        last if scalar @files == 0;
        if(scalar threads->list(threads::running) < $threads_trim){
                my $file = shift @files;
                threads->new(\&trim,$file);
                foreach my $t(threads->list(threads::joinable)){$t->join()}
        }
}
while(scalar threads->list(threads::running) > 0){sleep(1)}
foreach my $t(threads->list(threads::joinable)){$t->join()}

my @trimF = glob("$tmpdir/*.gz");

foreach $f(@trimF){
	my ($i)=$f=~/$tmpdir\/(.*)_trimmed/;
	my $out = $target.'/'.$i.'.clean.fq.gz';
	$data{$i}++;
	system("bowtie2 -p $threads_map -x $rRNA -U $f --un-gz $out -S mapped.sam");	
}


open OUT,'>clean.summary.txt';
print OUT "file\traw\tadaptor_removed\trRNA_removed\n";

foreach $k(sort keys %data){
	my $raw = `zgrep -c \$ $source/$k*`;
	my $trim = `zgrep -c \$ $tmpdir/$k*`;
	my $clean = `zgrep -c \$ $target/$k*`;
	$raw=~s/^\s+|\s+$//g;
	$trim=~s/^\s+|\s+$//g;
	$clean=~s/^\s+|\s+$//g;
	print OUT $k,"\t",$raw,"\t",$trim,"\t",$clean,"\n";
}
close OUT;

sub trim{
	my $name = shift;
	$name=~/$source\/(.*)\..*\.gz/;
	system("trim_galore --illumina --no_report_file --suppress_warn -o $tmpdir $name");
}

