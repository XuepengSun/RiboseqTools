#!/usr/bin/perl -w
use Bio::SeqIO;
use Getopt::Long;
use lib 'lib';
use log qw(LOG);
#use SXP::Util qw(LOG);
use List::Util qw(max min sum);
use Statistics::Basic qw(:all);
use Bio::Tools::CodonTable;
use Statistics::R;
use Term::ANSIColor;
use File::Basename;
use Cwd;
use File::chdir;
use File::Spec;


######################################################
#	contact: Xuepeng Sun (xs57@cornell.edu)	     #		
#	date:	         08/28/2017                  #       
######################################################

my ($genome,$gff,$cds,$prefix,$type,$ext_len,$mapping_qual,$mismatch,$seedLen,$read,$threads,$bam,$cdsFile,$readLenRange,$program,$codonRange,$outputDir,$minCount);
my ($strand,$engine,$bamRFCnt,$bamRFTrt,$bamRSCnt,$bamRSTrt,$src,$disp,$rdmin,$padj,$plot,$minP,$rfLen,$rsLen,$tool);

$program = {'map'=>1,'meta'=>1,'plot'=>1,'diff'=>1};

GetOptions(
	'genome=s'    => \$genome,
	'gff=s'		  => \$gff,
	'prefix=s'    => \$prefix,
	'type=s'	  => \$type,
	'ext=i'		  => \$ext_len,
	'mapq=i'	  => \$mapping_qual,
	'mis=i'		  => \$mismatch,
	'seed=i'	  => \$seedLen,
	'read=s'	  => \$read,
	'threads=i'   => \$threads,
	'bam=s'       => \$bam,
	'cds=s'       => \$cdsFile,
	'rdRange=s'   => \$readLenRange,
	'cdRange=s'   => \$codonRange,
	'dir=s'       => \$outputDir,
	'minCount=i'  => \$minCount,
	'strand=i'	  => \$strand,
	'engine=s'	  => \$engine,
	'bamRFCnt=s'  => \$bamRFCnt,
	'bamRFTrt=s'  => \$bamRFTrt,
	'bamRSCnt=s'  => \$bamRSCnt,
	'bamRSTrt=s'  => \$bamRSTrt,
	'src=s'		  => \$src,
	'disp=i'	  => \$disp,
	'rdmin=i'     => \$rdmin,
	'padj=s'      => \$padj,
	'plot=i'      => \$plot,
	'minP=f'	  => \$minP,
	'rfLen=s'     => \$rfLen,
	'rsLen=s'	  => \$rsLen,
	'tool=s'	  => \$tool
);

unless(scalar @ARGV == 1 && $program->{$ARGV[0]}){ERROR_INFOR("Check the \"task\" command!")}
my $task = shift @ARGV;

$outputDir //= 'RF_output';
$prefix //= 'riboseq';
$threads //= 10;

my $logFile = $prefix.'.log';
my $cwd = getcwd();
my $outWD = $cwd.'/'.$outputDir;
unless (-d $outWD){mkdir $outWD}


if($task eq 'map'){
	unless($genome && $gff && $read){ERROR_INFOR("Check the files: $genome, $gff, $read!")}
	$ext_len //= 45;
	$mapping_qual //= 20;
	$mismatch //= 1;
	$seedLen //= 22;   # default in bowtie2, but higher value may increase unique mapping 
	$type //= 'cds';
	
	$genome = File::Spec->rel2abs($genome);
	$gff = File::Spec->rel2abs($gff);
	my @readFile = split(/,/,$read);
	@readFile = map{File::Spec->rel2abs($_)} @readFile;
	
	my $mapDir = $outWD.'/map';
	unless(-d $mapDir){mkdir $mapDir}
	chdir $mapDir;
	
	my $index = &prepareGenomeAndIndex($genome,$gff,$type,$prefix,$ext_len);
	LOG("Building bowti2 index");
	unless(-d 'index'){mkdir "index"}
	system("bowtie2-build -q --threads $threads $index index/$prefix");
	
	foreach my $f(@readFile){
		LOG("Mapping file: $f");
		my ($filename,$dir) = fileparse($f);
		$filename=~/\.fq|\.fastq/;
		my $bamOut = $`.'.sorted.bam';
		system("bowtie2 -L $seedLen -p $threads -x index/$prefix -U $f 2>>$logFile | sambamba view -t 20 -f bam -S -F \"mapping_quality >= $mapping_qual and [NM] <=$mismatch\" /dev/stdin | sambamba sort -t 20 -o $bamOut /dev/stdin");
	}
}

elsif($task eq 'meta'){
	unless($bam && $cdsFile){ERROR_INFOR("Check the files: $bam,$cdsFile!")}
	
	$cdsFile = File::Spec->rel2abs($cdsFile);
	my @bams = split(/,/,$bam);
	@bams = map{File::Spec->rel2abs($_)} @bams;
	
	$readLenRange //= '27-30';
	$codonRange //= '-15,15';
	$minCount //= 20 ;
	
	my $valDir = $outWD.'/meta';
	unless(-d $valDir){mkdir $valDir}
	chdir $valDir;
	
	my ($metagene,$geneLenMatrics) = &readCDSFile($cdsFile,$codonRange);
	my @countRange = split(/-/,$readLenRange);
	
	foreach $bf(@bams){
		LOG("Processing bam file: $bf");
		my $bam_prefix = fileparse($bf,".bam");
		my (%phase,%bodyCov,%codonCov,$psite);
		open (BAM,"sambamba view $bf | ") || die "cannot find $bf $!";
		while(<BAM>){
			my @s = split;
			my $frame = $metagene->{$s[2]}{$s[3]}->{'frame'};
			$phase{length($s[9])}{$frame}++;
			if($geneLenMatrics->{$s[2]}->{'cds'} >= 150){     ## use genes >= 50AA
				$bodyCov{$s[2]}{$s[3]}++;
				if(length($s[9])>=$countRange[0] && length($s[9])<=$countRange[1]){
					if($frame == 0 || $frame == 1){       # get codon position for P site
						$psite = $metagene->{$s[2]}{$s[3]}->{'num'} + 4;
					}
					else{
						$psite = $metagene->{$s[2]}{$s[3]}->{'num'} + 5;
					}
					$codonCov{$s[2]}{$psite}++;	
				}
			}
		}
		close BAM;
		
		# output frame info
		open FRAME,">$bam_prefix\.frame";
		print FRAME "Length\tframe0\tframe1\tframe2\n";
		foreach $L(sort {$a<=>$b} keys %phase){
			print FRAME $L;
			for(0..2){
				$phase{$L}{$_}//=0;
				print FRAME "\t",$phase{$L}{$_};
			}
			print FRAME "\n";
		}
		close FRAME;
		
		# output P site distribution
		my (%codonStart,%codonEnd,%codonStartNorm,%codonEndNorm);
		my $gNum = 0;
		my @range = split(/,/,$codonRange);
		foreach $g(%codonCov){
			next if scalar keys %{$codonCov{$g}} < 20 ;   ## remove gene with < 20 codon covered
			next if sum (values %{$codonCov{$g}}) < $minCount ;     ### remove gene with < $minCount countable reads
			my $max = max (values %{$codonCov{$g}});          ### switch for average
			$gNum ++;
			for ($range[0]..29){    				# -15 + 30 codon
				$codonCov{$g}{$_}//=0;
				push @{$codonStartNorm{$_}},$codonCov{$g}{$_}/$max;		### switch for average
				$codonStart{$_}+=$codonCov{$g}{$_};				    ### switch for sum
			}
			my $codonTotal = $geneLenMatrics->{$g}->{'cds'} / 3;
			for($i=$codonTotal+$range[1];$i>=$codonTotal-29;$i--){
				$codonCov{$g}{$i-1}//=0;
				push @{$codonEndNorm{$i-$codonTotal}}, $codonCov{$g}{$i-1}/$max;	   ### switch for average
				$codonEnd{$i-$codonTotal}+=$codonCov{$g}{$i-1};				   ### switch for sum
			}
		}
		
		LOG("Total genes used for codon coverage: $gNum.");
		
		open CODON,">$bam_prefix\.Psite.start";
		print CODON "Postion\tPercentage\n";
		foreach(sort{$a<=>$b} keys %codonStart){
			print CODON $_,"\t",$codonStart{$_},"\n";
		}
		close CODON;
		
		open CODON,">$bam_prefix\.Psite.start.norm";
		print CODON "Postion\tPercentage\n";
		foreach(sort{$a<=>$b} keys %codonStartNorm){
			print CODON $_,"\t",mean(@{$codonStartNorm{$_}}),"\n";
		}
		close CODON;	
		
		open CODON,">$bam_prefix\.Psite.end";
		print CODON "Postion\tPercentage\n";
		foreach(sort{$a<=>$b} keys %codonEnd){
			print CODON $_,"\t",$codonEnd{$_},"\n";
		}
		close CODON;
		
		open CODON,">$bam_prefix\.Psite.end.norm";
		print CODON "Postion\tPercentage\n";
		foreach(sort{$a<=>$b} keys %codonEndNorm){
			print CODON $_,"\t",mean(@{$codonEndNorm{$_}}),"\n";
		}
		close CODON;
		
		# output genebody distribution
		my %mat_tmp=();
		foreach $g(keys %bodyCov){
			my $CDSLen = $geneLenMatrics->{$g}->{'cds'};
			my $urt5Len = $geneLenMatrics->{$g}->{'utr5'};
			my @s = map{int($CDSLen*$_ / 100) + $urt5Len} (1..100);
			my @dp_tmp;
			foreach(@s){
				$bodyCov{$g}{$_} //= 0;
				push @dp_tmp, $bodyCov{$g}{$_};
			}
			my $mean = mean(@dp_tmp);     ## use max to switch to max value
			next unless $mean > 0;
			my @s2 = map{$_/$mean} @dp_tmp;
			for(1..100){
				push @{$mat_tmp{$_}},$s2[$_-1];
			}
		}
		open BODY,">$bam_prefix\.bodyCov";
		print BODY "Percentile\tPercentage\n";
		foreach(sort {$a<=>$b} keys %mat_tmp){
			my $mean = mean(@{$mat_tmp{$_}});
			print BODY $_,"\t",$mean,"\n";
		}
		close BODY;
	}
}

elsif($task eq 'plot'){
	unless($bam){ERROR_INFOR("Check the bam files!")}
	my @bams = split(/,/,$bam);
	@bams = map{File::Spec->rel2abs($_)} @bams;
	my $valDir = $outWD.'/plot';
	unless(-d $valDir){mkdir $valDir}
	chdir $valDir;
	
	foreach $bfile(@bams){
		my $name = fileparse($bfile,".bam");
		RPlot($name,$outWD);
	}
}

elsif($task eq 'diff'){
	unless($bamRFCnt && $bamRFTrt && $bamRSCnt && $bamRSTrt && $cdsFile){ERROR_INFOR("Check the imput bam files!")}
	$strand //= 2;
	$rdmin //= 10;
	$rfLen //= '27-30';
	$rsLen //= '27-30';
	$engine //= 'ribodiff';
	
	my @bamsRFCnt = split(/,/,$bamRFCnt);
	my @bamsRFTrt = split(/,/,$bamRFTrt);
	my @bamsRSCnt = split(/,/,$bamRSCnt);
	my @bamsRSTrt = split(/,/,$bamRSTrt);
	
	my @rfLens = split(/-/,$rfLen);
	my @rsLens = split(/-/,$rsLen);
	
	$cdsFile = File::Spec->rel2abs($cdsFile);
	@bamsRFCnt = map{File::Spec->rel2abs($_)} @bamsRFCnt;
	@bamsRFTrt = map{File::Spec->rel2abs($_)} @bamsRFTrt;
	@bamsRSCnt = map{File::Spec->rel2abs($_)} @bamsRSCnt;
	@bamsRSTrt = map{File::Spec->rel2abs($_)} @bamsRSTrt;
		
	my $valDir = $outWD.'/diff';
	unless(-d $valDir){mkdir $valDir}
	chdir $valDir;

	my (%stranBias, %countRF, %countRS, %countPos);
	
	LOG("Reading $cdsFile");
	open CDS,$cdsFile || die "$!";
	while(<CDS>){
		if(my ($gid,$gleft,$glen)=$_=~/^>(\S+).*left=(\d+).*cds_len=(\d+)/){
			$countPos{$gid} = {'min'=>$gleft + 1, 'max'=>$gleft + $glen};
		}
	}
	close CDS;
	
	my (@allRF,@allRS);
	push (@allRF,@bamsRFCnt);
	push (@allRF,@bamsRFTrt);
	push (@allRS,@bamsRSCnt);
	push (@allRS,@bamsRSTrt);
	
	foreach(@allRF){
		my $bam_prefix = fileparse($_,".bam");
		LOG("Counting $bam_prefix");
		open BAM,"sambamba view $_ |";
		while(<BAM>){
			my @s = split;
			next if (length($s[9]) > $rfLens[1] || length($s[9]) < $rfLens[0]);
			if($s[3] + length($s[9]) -1 >= $countPos{$s[2]}->{'min'} && $s[3] <= $countPos{$s[2]}->{'max'}){
				$countRF{$bam_prefix}{$s[2]}{$s[1]}++;
			}
		}
		close BAM;
	}
	
	foreach(@allRS){
		my $bam_prefix = fileparse($_,".bam");
		LOG("Counting $bam_prefix");
		open BAM,"sambamba view $_ |";
		while(<BAM>){
			my @s = split;
			next if (length($s[9]) > $rsLens[1] || length($s[9]) < $rsLens[0]);
			if($s[3] + length($s[9]) -1 >= $countPos{$s[2]}->{'min'} && $s[3] <= $countPos{$s[2]}->{'max'}){
				$countRS{$bam_prefix}{$s[2]}{$s[1]}++;
			}
		}
		close BAM;
	}
	
	my @tmpBams;
	push (@tmpBams,sort keys %countRF);
	push (@tmpBams,sort keys %countRS);
	
	open ALL,">All.count";
	print ALL "gene\t",join("\t",@tmpBams),"\n";
	foreach $g(keys %countPos){
		print ALL $g;
		foreach $b(@tmpBams){
			my ($pos,$neg);
			if($countRF{$b}){
				$countRF{$b}{$g}{0} //= 0;    ## positive strand
				$countRF{$b}{$g}{16} //= 0;   ## negative strand
				$stranBias{$b}{0}+=$countRF{$b}{$g}{0};
				$stranBias{$b}{16}+=$countRF{$b}{$g}{16};
				$pos = $countRF{$b}{$g}{0};
				$neg = $countRF{$b}{$g}{16};
			}
			else{
				$countRS{$b}{$g}{0} //= 0;    ## positive strand
				$countRS{$b}{$g}{16} //= 0;   ## negative strand
				$stranBias{$b}{0}+=$countRS{$b}{$g}{0};
				$stranBias{$b}{16}+=$countRS{$b}{$g}{16};
				$pos = $countRS{$b}{$g}{0};
				$neg = $countRS{$b}{$g}{16};			
			}
			
			if($strand == 0){  ## positive strand
				print ALL "\t",$pos;
			}
			elsif($strand == 1){
				print ALL "\t",$neg;
			}
			elsif($strand == 2){
				print ALL "\t",$pos+$neg;
			}
			else{ERROR_INFOR("Check strand value!")}
		}
		print ALL "\n";
	}
	close ALL;
	
	my $RFnum = scalar @allRF;
	open RF,">Riboseq.count";
	open RS,">mRNAseq.count";
	open IN,"All.count" || die "$!";
	while(<IN>){
		chomp;
		my @ss = split;
		my $rf = join("\t",@ss[0..$RFnum]);
		my $rs = join("\t",@ss[0,$RFnum+1..$#ss]);
		print RF $rf,"\n";
		print RS $rs,"\n";
	}
	close IN;
	close RF;
	close RS;

	foreach my $b(sort keys %stranBias){
		my $log = "Mapping strandness for ".$b." - Postive: ".$stranBias{$b}{0}." Negtitve: ".$stranBias{$b}{16};
		LOG($log);
	}
	
	open CON,">condition.csv";
	print CON "Samples,Data_Type,Conditions\n";
	foreach(@bamsRFCnt){
		my $name = fileparse($_,".bam");
		print CON $name,",Ribo-Seq,Control\n";
	}
	foreach(@bamsRFTrt){
		my $name = fileparse($_,".bam");
		print CON $name,",Ribo-Seq,Treated\n";
	}
	foreach(@bamsRSCnt){
		my $name = fileparse($_,".bam");
		print CON $name,",RNA-Seq,Control\n";
	}	
	foreach(@bamsRSTrt){
		my $name = fileparse($_,".bam");
		print CON $name,",RNA-Seq,Treated\n";
	}	
	close CON;
	
	LOG("Running DESeq2");
	DESeq2('mRNAseq.count','mRNA',scalar @bamsRSCnt,scalar @bamsRSTrt,$rdmin);
	DESeq2('Riboseq.count','Ribo',scalar @bamsRFCnt,scalar @bamsRFTrt,$rdmin);

	LOG("Running $engine");
	if($engine eq 'ribodiff'){
		$src //= 'utils/TE.py';
		$disp //= 0;
		$padj //= 'BH';
		$plot //= '0';
		$minP //= '0.1';
		system("python $src -e condition.csv -c All.count -o TE.ribodiff.out -d $disp -s $rdmin -m $padj -p $plot -q $minP 2>ribodiff.log");
	}
	elsif($engine eq 'riborex'){
		$tool //= 'DESeq2';
		riborex(scalar @bamsRSCnt,scalar @bamsRSTrt,scalar @bamsRFCnt,scalar @bamsRFTrt,$tool);
	}
	else{ERROR_INFOR("check diff options!")}
}



sub readCDSFile{
	LOG("Reading CDS file.");
	my ($cds_file,$codon_range) = @_;
	my (%hash,%gene_length);
	my @range = split(/,/,$codon_range);
	open (CDS,"< $cds_file") || die "cannot find $cds_file $!";
	my @content = <CDS>;
	close CDS;
	die "$cds_file is not complete at the end of file!" if (scalar (@content) % 2 != 0); 
	my $CodonTable = Bio::Tools::CodonTable->new();
	while(my @tp = splice(@content,0,2)){
		chomp(@tp);
		if(my ($geneID,$left,$right,$len) = $tp[0]=~/^>(\S+)\s+left=(\d+)\s+right=(\d+).*cds_len=(\d+)/){
			$gene_length{$geneID}={'cds'=>$len,'utr5'=>$left,'utr3'=>$right};
			my $count = 0;
			for(my $i=$left; $i>=3; $i-=3){
				$count--;
				last if $count < $range[0];
				my $codon = $CodonTable->translate(substr($tp[1],$i-3,3));
				$hash{$geneID}{$i} = {'frame'=>2, 'codon'=>$codon, 'num'=>$count};
				$hash{$geneID}{$i-1} = {'frame'=>1, 'codon'=>$codon, 'num'=>$count};
				$hash{$geneID}{$i-2} = {'frame'=>0, 'codon'=>$codon, 'num'=>$count};
			}
			$count = 0;
			for(my $i=$left+1; $i<=$left+$len+$right-2; $i+=3){
				last if ($count + 1 - $len / 3) > $range[1];
				my $codon = $CodonTable->translate(substr($tp[1],$i-1,3));
				$hash{$geneID}{$i} = {'frame'=>0, 'codon'=>$codon, 'num'=>$count};
				$hash{$geneID}{$i+1} = {'frame'=>1, 'codon'=>$codon, 'num'=>$count};
				$hash{$geneID}{$i+2} = {'frame'=>2, 'codon'=>$codon, 'num'=>$count};
				$count++;
			}
		}
	}
	return (\%hash,\%gene_length);
}

sub prepareGenomeAndIndex{
	my ($genome_file,$gff_file,$seq_type,$output,$extraSeqLen) = @_;
	LOG("Reading GFF file");
	my (%gene,%mRNA,%filter,%gseq,%strand);
	open(GFF,"<$gff_file") || die;
	while(<GFF>){
		last if /^##FASTA/;
		next if /^#/;
		my @s = split /\t/;
		if($s[2] eq 'mRNA'){
			my ($mm,$gg) = $s[8]=~/ID=(.*?);.*Parent=(.*?)[;\n]/;
			$gene{$gg}{$mm}=abs($s[4]-$s[3]);
			$strand{$mm} = $s[6];
		}
		elsif($s[2] eq 'CDS'){
			$s[8]=~/Parent=(.*?)[;\n]/;
			$mRNA{$s[0]}{$1}{min($s[3],$s[4])} = max($s[3],$s[4]);
		}
		else{next}
	}
	close GFF;
	foreach $g(keys %gene){my @key=sort{$gene{$g}{$b}<=>$gene{$g}{$a}} keys %{$gene{$g}};$filter{$key[0]}++}	
	LOG("Filtering isoforms: ".(scalar keys %filter)." retained");
	my $io = new Bio::SeqIO(-file=>$genome_file,-format=>'fasta');
	while(my $seq=$io->next_seq){$gseq{$seq->id}=$seq->seq}
	$io->close;
	my %cutRegion=();
	my $return = $output.'.genome.fa';
	open OUT,">$return";	
	foreach my $chrm(keys %mRNA){
		foreach my $rna(keys %{$mRNA{$chrm}}){
			next if !$filter{$rna};
			my $cds_tmp = '';
			foreach $left(sort{$a<=>$b} keys %{$mRNA{$chrm}{$rna}}){
				$cds_tmp.=substr($gseq{$chrm},$left-1,$mRNA{$chrm}{$rna}{$left}-$left+1);
			}
			my $start_codon = substr($cds_tmp,0,3);
			my $stop_codon = substr($cds_tmp,-3,3);
			my $nt_len = length($cds_tmp);
			my $st = min(keys %{$mRNA{$chrm}{$rna}})-$extraSeqLen < 1 ? 1 : min(keys %{$mRNA{$chrm}{$rna}})-$extraSeqLen;
			my $st_len = min(keys %{$mRNA{$chrm}{$rna}}) - $st;
			my $end = max(values %{$mRNA{$chrm}{$rna}})+$extraSeqLen > length($gseq{$chrm}) ? length($gseq{$chrm}) : max(values %{$mRNA{$chrm}{$rna}})+$extraSeqLen;
			my $end_len = $end - max(values %{$mRNA{$chrm}{$rna}});
			$cds_tmp = substr($gseq{$chrm},$st-1,$st_len).$cds_tmp.substr($gseq{$chrm},max(values %{$mRNA{$chrm}{$rna}}),$end_len);
			if($strand{$rna} eq '-'){
				$cds_tmp=join("",reverse split("",$cds_tmp));
				$cds_tmp=~tr/ATGCatgc/TACGtacg/;
				$st_tmp=join("",reverse split("",$start_codon));
				$st_tmp=~tr/ATGCatgc/TACGtacg/;
				$sp_tmp=join("",reverse split("",$stop_codon));
				$sp_tmp=~tr/ATGCatgc/TACGtacg/;
				($start_codon,$stop_codon)=($sp_tmp,$st_tmp);
				($st_len,$end_len) = ($end_len,$st_len);
			}
			print OUT '>',$rna," left=$st_len right=$end_len start_codon=$start_codon stop_codon=$stop_codon cds_len=$nt_len\n",$cds_tmp,"\n";
			$cutRegion{$chrm}{$st} = $end;
		}
	}
	if($seq_type eq 'genome'){
		foreach my $c(keys %gseq){
			my @start = sort{$a<=>$b} keys %{$cutRegion{$c}};
			my $plen = $start[0] - 1;
			if($plen > 0){
				my $string = substr($gseq{$c},0,$plen);
				my $string_ck = $string;
				$string_ck=~s/N//gi;
				if(length ($string_ck) >= 50 ){   ## removing gap-only sequences;
					print OUT ">",$c,':1-',$start[0]-1,"\n",$string,"\n";
				}
			}
			for(my $i=0;$i<$#start;$i++){
				next if $cutRegion{$c}{$start[$i]} >= $start[$i+1];
				my $len = $start[$i+1] -1 - $cutRegion{$c}{$start[$i]};
				my $string = substr($gseq{$c},$cutRegion{$c}{$start[$i]},$len);
				my $string_ck = $string;
				$string_ck=~s/N//gi;
				if(length ($string_ck) >= 50 ){   ## removing gap-only sequences;
					print OUT ">",$c,':',$cutRegion{$c}{$start[$i]}+1,'-',$start[$i+1]-1,"\n",$string,"\n";
				}
			}
			my $llen = length($gseq{$c}) - $cutRegion{$c}{$start[-1]};
			if($llen > 0){
				my $string = substr($gseq{$c},$cutRegion{$c}{$start[-1]},$llen);
				my $string_ck = $string;
				$string_ck=~s/N//gi;
				if(length ($string_ck) >= 50 ){   ## removing gap-only sequences;
					print OUT ">",$c,":",$cutRegion{$c}{$start[-1]}+1,"-",length($gseq{$c}),"\n",$string,"\n";
				}
			}
		}
	}
	close OUT;
	return $return;
}

sub ERROR_INFOR{
	my $print = shift;
	print STDERR "\n====================================================\n";
	print STDERR "Error:  ",color("green"),$print,color("reset"),"\n";
	print STDERR "====================================================\n";
	print STDERR <DATA>;
	exit;
}

sub RPlot{
	my ($file,$dir) = @_;
	my $R = Statistics::R->new( shared => 1);

my $cmds = <<EOF;
	suppressMessages(library(ggplot2))
	suppressMessages(library(reshape))
	suppressMessages(library(gridExtra))
	
	dframe = paste("$dir","/meta/","$file",".frame",sep="")
	bodycov = paste("$dir","/meta/","$file",".bodyCov",sep="")
	pfive = paste("$dir","/meta/","$file",".Psite.start",sep="")
	pfiveNorm = paste("$dir","/meta/","$file",".Psite.start.norm",sep="")
	pthree = paste("$dir","/meta/","$file",".Psite.end",sep="")
	pthreeNorm = paste("$dir","/meta/","$file",".Psite.end.norm",sep="")

	outFile = paste("$dir","/plot/","$file",".pdf",sep="")

	data <- read.table(dframe,header=T)
	sortData <- subset(data, Length>=27 & Length <=30)
	mdata <- melt(sortData, id=c("Length"))
	gframe <- ggplot(data=mdata, aes(x=Length, y=value, fill=variable)) +
		geom_bar(stat="identity", color="black", position=position_dodge())+
		theme_minimal() +
		xlab("Read length") +
		ylab("Count") +
		ggtitle("3nt Periodicity") +
		scale_fill_brewer(palette="Greens",name="Frame") + 
		theme(plot.title = element_text(hjust = 0.5))
		
	
	data <- read.table(bodycov,header=T)
	gbodycov <- ggplot(data=data, aes(x=Percentile, y=Percentage)) +
		geom_line(size=1.5,color="#009E73") +
		theme_minimal() +
		xlab("Position") +
		ylab("Relative abundance") +
		ggtitle("Gene body Coverage") +
		ylim(0,1.5*max(data\$Percentage)) + 
		theme(plot.title = element_text(hjust = 0.5))
		
	data <- read.table(pfive,header=T)
	gpFive <- ggplot(data=data, aes(x=Postion, y=Percentage)) +
		geom_line(size=1.5,color="#009E73") +
		theme_minimal() +
		xlab("Codon poistion (5')") +
		ylab("Count") +
		ggtitle("P site distribution") +
		ylim(0,1.5*max(data\$Percentage)) + 
		theme(plot.title = element_text(hjust = 0.5))	
		
	data <- read.table(pfiveNorm,header=T)
	gpfiveNorm <- ggplot(data=data, aes(x=Postion, y=Percentage)) +
		geom_line(size=1.5,color="#009E73") +
		theme_minimal() +
		xlab("Codon poistion (5')") +
		ylab("Normalized abundance") +
		ggtitle("P site distribution") +
		ylim(0,1) + 
		theme(plot.title = element_text(hjust = 0.5))
		
	data <- read.table(pthree,header=T)
	gpthree <- ggplot(data=data, aes(x=Postion, y=Percentage)) +
		geom_line(size=1.5,color="#009E73") +
		theme_minimal() +
		xlab("Codon poistion (3')") +
		ylab("Count") +
		ggtitle("P site distribution") +
		ylim(0,1.5*max(data\$Percentage)) + 
		theme(plot.title = element_text(hjust = 0.5))	
		
	data <- read.table(pthreeNorm,header=T)
	gpthreeNorm <- ggplot(data=data, aes(x=Postion, y=Percentage)) +
		geom_line(size=1.5,color="#009E73") +
		theme_minimal() +
		xlab("Codon poistion (3')") +
		ylab("Normalized abundance") +
		ggtitle("P site distribution") +
		ylim(0,1) + 
		theme(plot.title = element_text(hjust = 0.5))
	
	pdf(file=outFile)
	grid.arrange(gframe,gbodycov,gpFive,gpfiveNorm,gpthree,gpthreeNorm, ncol=2)
	dev.off()
EOF
	$R->run($cmds);
	$R->stop;
}

sub DESeq2{
	my ($file,$type,$ctrNum,$trtNum,$minReads) = @_;
	my $R = Statistics::R->new( shared => 1);
	my $stats_out = $type.'.DESeq2.out';
	my $pca = $type.'.PCA.pdf';
	
my $cmds = <<EOF;
	suppressMessages(library("DESeq2"))
	data<-read.table("$file",header=T,row.names=1)
	pasillaCountTable <- data[which(rowSums(data)>=$minReads),]

	pasillaDesign = data.frame(
		row.names = colnames(pasillaCountTable),
		condition = c(rep("control",$ctrNum),rep("treated",$trtNum)),
		libType = c(rep("SE",$ctrNum+$trtNum))
	)

	colData <- pasillaDesign[,c("condition","libType")]
	dds <- DESeqDataSetFromMatrix(countData = pasillaCountTable,
		colData = colData,
		design = ~ condition
	)

	dds <- DESeq(dds)
	res <- results(dds, contrast = c("condition", "treated","control"))
	write.table(res,file="$stats_out",quote=F,sep="\t",row.name=T)

	rld <- rlog(dds, blind = FALSE)

	pdf(file="$pca", useDingbats=FALSE) 
	plotPCA(rld, intgroup = c("condition"))
	dev.off()
EOF
	$R->run($cmds);
	$R->stop;
}

sub riborex{       
### minimal read count cutoff is not applied here, because riborex needs same gene number in RS and RF
	my ($RSctrNum,$RStrtNum,$RFctrNum,$RFtrtNum,$tool) = @_;
	my $R = Statistics::R->new( shared => 1);
	my $output = 'riborex.'.$tool.'.output';
	
my $cmds = <<EOF;
	suppressMessages(library("riborex"))

	RNACntTable <- read.table('mRNAseq.count',header=T,row.names=1)
	RiboCntTable <- read.table('Riboseq.count',header=T,row.names=1)

	riboCond <- c(rep("control",$RFctrNum),rep("treated",$RFtrtNum))
	rnaCond <- c(rep("control",$RSctrNum),rep("treated",$RStrtNum))

	res <- riborex(RNACntTable,RiboCntTable,rnaCond,riboCond,"$tool")
	write.table(res,file="$output",quote=F,sep="\t",row.name=T)
EOF

	$R->run($cmds);
	$R->stop;
}


__DATA__

Useage:
	perl Riboseq_tools.pl [task] [options]

Task: map meta plot diff

External tools required: bowtie2 sambamba 

options:

Global parameters:
	
	-dir      output dir, default: RF_output
	-prefix   prefix of output, default: riboseq.
	-threads  threads number, default: 10.
	
** map **
This command creates a modified cds or intro-free genome for reads mapping, 
but only cds sequences are used for downstream analyses.

	-read	  reads file,compressed or not,separated by comma, eg: read1.fq,read2,fq [required]   
	-genome   genome file in fasta format [required]
	-gff      gff file [required]
	-type     mappting to cds or intron-free genome, options: cds,genome. Default: cds
	-ext      N bp upstream and downstream sequences to be added to CDS, default: 45	
	-mapq     minimal mapping quality, default: 20
	-mis      maximal mismatch allowed, default: 1
	-seed     seed length for bowtie2, default: 22

** meta **
This command evaluates the patterns of phasing, reads coverage and codon enrichment.
Figures are generated automatically.

	-bam      sorted bam files, separated by comma, eg: bam1,bam2 [required]
	-cds      cds sequence file generated by 'map' command [required]
	-rdRange  read length ranges for codon analysis, format: num-num. Default: 27-30
	-cdRange  codon range to count, cannot > ("ext"/3), format: neg,pos. Default: -15,15
	-minCount minimal count of a gene to use, default: 20
	
** plot **
Plot data generated by meta command, one bam file one graph (.pdf)

	-bam      bam files used in meta command, separated by comma, eg: bam1,bam2 [required]

** diff **
Calculate differential expression/translation using DESeq2 (for mRNA/RFP) or Ribodiff (for TE)
TE can also be calculated using Riborex, which is faster but cannot be applied when the number
of samples are different for mRNA and RFP.

	-cds      cds sequence file generated by 'map' command [required]
	-rfLen    reads within this range to count for Riboseq, default: 27-30
	-rsLen    reads within this range to count for mRNAseq, default: 27-30
	-strand   mappng strand to count: +: 0; -: 1; both: 2. Default: 2
	-engine	  tool to be used for differential TE analysis, ribodiff or riborex. Default: ribodiff
	-bamRFCnt bam files of riboseq used as control, separated by comma. [required]
	-bamRFTrt bam files of riboseq used as treatmeat. [required]
	-bamRSCnt bam files of mRNAseq used as control. [required]
	-bamRSTrt bam files of mRNAseq used as treatmeat. [required]
	-rdmin	  minimal sum of normalized reads for the test. default: 10

options for ribodiff
	
	-src	  location for ribodiff script (TE.py), defualt: utils/TE.py
	-disp	  dispersion for riboseq or mRNAseq. 1: different dispersion; 0: same dispersion. default: 0
	-padj	  method for multiple test correction, BH, Bonferroni. Default: BH
	-plot 	  make plots to show the data and results. On: 1; Off: 0. Defulat: 0
	-minP     FDR cutoff for plot. Default: 0.1

options for riborex
		
	-tool	  tools for estimation, options: DESeq2, edgeR, edgeRD, Voom. Default: DESeq2.

