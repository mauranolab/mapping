#!/usr/bin/perl -w
#from /home/ehaugen/work/deepAIupdate/comprehensiveGenotyping4/src/parseSamtoolsGenotypesToBedFiles.pl

use strict;
use IO::File;

if (@ARGV < 5) {
    die "Usage:  $0 file.vcf[.gz] prefix_for_unsorted_output min_site_quality_threshold min_total_depth_threshold min_allele_depth_threshold\n";
}
#my $vcfInput = $ARGV[0];
my ($vcfInput, $outputPrefix, $min_site_quality_threshold, $min_total_depth_threshold, $min_allele_depth_threshold) = @ARGV;

my @sampleNames = undef;
my $NUM_SAMPLES = undef;
my %outs = ();

print "Parsing $vcfInput. min_site_quality_threshold=$min_site_quality_threshold, min_total_depth_threshold=$min_total_depth_threshold, min_allele_depth_threshold=$min_allele_depth_threshold\n";

open IN, "zcat -f $vcfInput |" or die "$!";
while (<IN>) {
    chomp;
    next if (/^\#\#/);
    my $inputLine = $_;
    if (/^#CHROM\tPOS/) {
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/AG10803-DS12374.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/GM06990-DS7748.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HAEpiC-DS12663.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCF-DS12501.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCM-DS12599.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCPEpiC-DS12447.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HIPEpiC-DS12684.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HRCE-DS10666.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/IMR90-DS13219.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/NHA-DS12800.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/SAEC-DS10518.hg19.bam
        my ($chromLABEL, $posLABEL, $idLABEL, $refLABEL, $altLABEL, $qualLABEL, $filterLABEL, $infoLABEL, $formatLABEL, @sampleNamesList ) = split /\t/;
        @sampleNames = sampleNamesFromFilenames( @sampleNamesList );  
        $NUM_SAMPLES = scalar @sampleNames;
        foreach my $name (@sampleNames) {
            $outs{$name} = new IO::File "> $outputPrefix.$name.txt";
        }
    } elsif (/^#/) {
        # header comment 
    } else {
        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @genotypes ) = split /\t/;
        if ($qual < $min_site_quality_threshold) {
            warn "Skipping $chrom $pos with site quality $qual < $min_site_quality_threshold\n";
            next;
        }
        my %code2base = ( 0 => $ref, 1 => $alt, '.' => 'N' );
#        if (length($ref) != 1 or length($alt) != 1) {
#            warn "Skipping $ref/$alt at $chrom $pos\n";
#            next;
#        }
        my $totalDP = undef;
        my($numRef, $numNonRef) = undef;
        foreach my $infoItem (split /\;/, $info) {
            my ($key,$value) = split /\=/, $infoItem;
            if ("DP" eq $key) {
                $totalDP = int($value);
            } elsif ("DP4" eq $key) {
                #NB summed across all indivs, so only useful for single-sample case
                #eg DP4=41,47,16,14
                my ($refFwd, $refRev, $nonrefFwd, $nonrefRev) = split /,/, $value;
                $numRef = $refFwd+$refRev;
                $numNonRef = $nonrefFwd+$nonrefRev;
            }
        }
        unless (defined($totalDP) and defined($numRef) and defined($numNonRef) ) {
            die "Failed to parse DP or DP4 from INFO ($info) in $_\n";
        }
        #if ($totalDP < $min_total_depth_threshold) {
        #    #warn "Skipping $chrom $pos with total depth $totalDP < $min_total_depth_threshold\n";
        #    next;
        #}
        if ($min_allele_depth_threshold > 0 and ($numRef < $min_allele_depth_threshold or $numNonRef < $min_allele_depth_threshold)) {
            #warn "Skipping $chrom $pos with ref/nonref depth $numRef/$numNonRef < $min_allele_depth_threshold\n";
            next;
        }
        my $min0 = $pos - 1;
        my $max1 = $pos - 1 + length($ref) ;
        unless ($chrom =~ /^chr/) {
            $chrom = "chr$chrom";
        }
        my $bed3 = "$chrom\t$min0\t$max1";
        for (my $i = 0; $i < $NUM_SAMPLES; ++$i) {
            my $sampleName = $sampleNames[$i];
            my ($betterGT, $GQ, $DP, $PL) = parseGenotype( $format, $genotypes[$i], \%code2base );
            next if ("N/N" eq $betterGT); # parsed from "./."
            next if ($DP < $min_total_depth_threshold); # Require minimum per-sample depth of 8
            my $outFH = $outs{$sampleName}; 
            print $outFH "$bed3\t$id\t$GQ\t$ref\t$alt\t$betterGT\t$PL\t$DP\t$numRef\t$numNonRef\n";
        }
    }
}
close IN;

# close outputs
foreach my $key (keys %outs) {
    delete $outs{$key};
}
%outs = ();
undef %outs;
warn "Output files closed safely?\n";
exit;

sub sampleNamesFromFilenames { 
    my @filenames = @_;
    my @result = ();
    foreach my $filename (@filenames) {
        #chomp( my $name = `basename $filename | cut -f1 -d .` );
        chomp( my $name = `basename $filename` );
        if ($name =~ /^([^\/\.]+)\.bam$/) {
            $name = $1;
        } else {
            $name = $filename;
        }
        #warn "$name\t$filename\n";
        push @result, $name;
    }
    return @result;
}

sub parseGenotype {
    my ($format, $genotype, $code2base) = @_;
    my ($GT, $PL, $qualityScore, $readDepth);
    if ($format eq "GT:PL:GQ") {
        #my ($GT,$PL,$GQ) = split /\:/, $genotype;
        # Skip this, now requiring read-depth DP field
        die "Unsupported genotype format ($format)\n"; 
    } elsif ($format eq "GT:PL:DP:AD:GQ") {
        my ($GT1,$PL1,$DP,$AD,$GQ) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = $PL1;
    } elsif ($format eq "GT:PL:DP:SP:GQ") {
        my ($GT1,$PL1,$DP,$SP,$GQ) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = $PL1;
    } elsif ($format eq "GT:PL:DP:GQ") {
        my ($GT1,$PL1,$DP,$GQ) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = $PL1;
    } elsif ($format eq "GT:PL:DP:DV:DP4:DPR:GQ") {
        my ($GT1,$PL1,$DP,$DV,$DP4,$DPR,$GQ) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = $PL1;
    } elsif ($format eq "GT:GT:PL:DP:SP:GQ:GQ") {
        #my ($GT1,$GT2,$PL,$DP,$SP,$GQ1,$GQ2) = split /\:/, $genotype;
        # oops it's actually mangled, with duplicates
        my ($GT1,$PL1,$DP,$SP,$GQ1,$GT2,$GQ2) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ1;
        $GT = $GT1;
        $PL = $PL1;
        if (($GT1 ne $GT2) or ($GQ1 != $GQ2)) {
            die "Wrong parsing assumptions for VCF genotype ($genotype)\n";
        }
    } else {
        die "Unsupported genotype format ($format)\n";
    }
    my @alleles = ();
    foreach my $code (split /\//, $GT) {
        my $base = $code2base->{$code};
        unless (defined($base)) {
            die "Failed to parse genotype ($GT)\n";
        }
        push @alleles, $base;
    }
    my $baseSlashBaseGenotype = join("/", @alleles);
    return ($baseSlashBaseGenotype, $qualityScore, $readDepth, $PL);
}



##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/AG10803-DS12374.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/GM06990-DS7748.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HAEpiC-DS12663.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCF-DS12501.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCM-DS12599.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HCPEpiC-DS12447.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HIPEpiC-DS12684.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/HRCE-DS10666.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/IMR90-DS13219.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/NHA-DS12800.hg19.bam	/net/glados/solexa_data/rsandstrom/all_stamlab_alignments/alignments/SAEC-DS10518.hg19.bam
#   10	61321	.	A	G	4.22	.	DP=11;VDB=2.880000e-02;RPB=8.249579e-01;AF1=0.1296;AC1=2;DP4=4,5,1,1;MQ=25;FQ=4.61;PV4=1,1,1,1	GT:PL:GQ	0/0:0,0,0:8	0/0:0,0,0:8	0/0:0,3,25:11	0/0:0,3,25:11	0/0:0,0,0:8	0/1:38,0,34:30	0/0:0,9,62:17	0/0:0,0,0:8	0/0:0,3,25:11	0/0:0,3,25:11	0/0:0,0,0:8
#   10	78106	.	C	A	145	.	DP=14;VDB=2.823745e-03;RPB=-1.200000e+00;AF1=0.5526;AC1=12;DP4=3,2,2,7;MQ=25;FQ=60.6;PV4=0.27,0.069,1,1	GT:PL:GQ	1/1:25,3,0:3	0/1:0,0,0:3	0/1:0,3,25:3	0/1:0,3,25:3	0/1:19,0,19:19	0/1:0,3,25:3	1/1:50,6,0:5	0/1:0,0,0:3	1/1:88,12,0:10	0/1:0,3,25:3	1/1:25,3,0:3
#   10	78118	.	C	T	77.3	.	DP=11;VDB=1.743144e-02;RPB=8.723915e-01;AF1=0.4463;AC1=10;DP4=0,5,2,4;MQ=25;FQ=58.8;PV4=0.45,0.31,1,1	GT:PL:GQ	0/0:0,3,25:4	0/1:0,0,0:3	0/1:0,0,0:3	0/0:0,3,25:4	0/1:25,3,0:4	0/0:0,3,25:4	0/1:19,0,19:19	0/1:0,0,0:3	0/1:59,0,13:18	0/1:0,0,0:3	0/1:25,3,0:4
#   10	95220	.	A	G	20.6	.	DP=9;VDB=6.737978e-02;RPB=-6.510653e-01;AF1=0.391;AC1=8;DP4=2,4,2,1;MQ=25;FQ=22.8;PV4=0.52,1,1,1	GT:PL:GQ	0/0:0,3,25:4	0/1:25,3,0:4	0/1:25,3,0:4	0/1:0,0,0:3	0/1:0,0,0:3	0/1:0,0,0:3	0/1:0,0,0:3	0/0:0,3,25:4	0/0:0,9,71:9	0/1:0,0,0:3	0/1:19,0,19:18
#   10	95289	.	C	A	145	.	DP=30;VDB=1.160764e-01;RPB=-1.740894e+00;AF1=0.4986;G3=5.303e-07,1,4.213e-07;HWE=0.00046;AC1=11;DP4=6,10,6,6;MQ=25;FQ=148;PV4=0.7,0.24,1,1	GT:PL:GQ	0/1:37,0,16:19	0/1:16,0,41:19	0/0:0,3,25:3	0/1:24,0,71:27	0/1:19,0,19:19	0/1:16,0,37:19	0/1:19,0,19:19	0/1:16,0,37:19	0/1:19,0,19:19	0/1:19,0,19:19	0/1:25,3,0:3
#   10	95301	.	G	T	22.8	.	DP=40;VDB=1.050387e-01;RPB=-2.653228e+00;AF1=0.1692;AC1=3;DP4=11,21,3,4;MQ=25;FQ=23.6;PV4=0.69,0.41,1,1	GT:PL:GQ	0/0:0,6,45:10	0/1:13,0,51:9	0/0:16,22,61:10	0/0:0,18,125:22	0/1:61,22,16:6	0/0:0,12,75:16	0/0:0,3,25:8	0/0:0,12,81:16	0/1:24,0,89:20	0/0:0,9,70:13	0/0:0,6,45:10
#   10	95306	.	C	A	145	.	DP=34;VDB=1.207853e-01;RPB=1.032875e+00;AF1=0.4166;AC1=9;DP4=4,15,7,4;MQ=25;FQ=148;PV4=0.047,0.12,1,1	GT:PL:GQ	0/1:25,3,0:4	0/0:0,12,75:11	0/0:0,6,50:6	0/1:31,0,56:32	0/0:0,9,70:8	0/0:0,6,45:6	0/1:25,3,0:4	0/1:19,0,19:19	0/1:53,0,46:49	0/1:19,0,19:19	1/1:45,6,0:4
#   10	97552	.	A	G	114	.	DP=36;VDB=2.029366e-03;RPB=-1.030255e+00;AF1=0.365;G3=0.04761,0.9524,2.26e-112;HWE=0.0103;AC1=8;DP4=13,12,7,4;MQ=25;FQ=116;PV4=0.72,1,1,1	GT:PL:GQ	0/0:0,6,50:7	0/1:16,0,37:17	0/0:0,3,25:4	0/1:13,0,51:14	0/0:0,3,25:4	0/1:34,0,38:34	0/1:16,0,41:17	0/1:41,0,13:19	0/1:10,0,57:11	0/1:47,0,86:48	0/0:0,6,50:7
#   10	97560	.	C	T	9.38	.	DP=31;VDB=2.649457e-02;RPB=5.597929e-01;AF1=0.1676;AC1=3;DP4=18,9,1,3;MQ=25;FQ=10.1;PV4=0.27,0.21,1,1	GT:PL:GQ	0/0:0,3,25:9	0/0:0,9,62:14	0/0:0,0,0:6	0/0:0,12,86:17	0/0:0,3,25:9	0/0:0,12,86:17	0/0:0,9,70:14	0/1:13,0,59:9	0/1:13,0,51:9	0/1:35,0,49:30	0/0:0,6,50:11
#   10	98706	.	A	G	18.3	.	DP=20;VDB=2.063840e-02;RPB=-8.468098e-01;AF1=0.1715;AC1=3;DP4=8,9,1,2;MQ=25;FQ=19.1;PV4=1,1,1,1	GT:PL:GQ	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,6,45:11	0/0:0,3,25:8	0/1:25,3,0:8	0/1:41,0,16:26	0/0:0,21,130:25	0/0:0,9,71:13	0/0:0,3,25:8	0/0:0,6,45:11
#   10	98717	.	G	A	15.9	.	DP=18;VDB=2.063840e-02;RPB=1.184698e-01;AF1=0.1778;AC1=4;DP4=6,9,1,2;MQ=25;FQ=16.7;PV4=1,1,1,1	GT:PL:GQ	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,3,25:8	0/1:19,0,19:15	0/1:41,0,16:26	0/0:0,18,120:22	0/0:0,9,71:13	0/0:0,3,25:8	0/0:0,6,45:10
#   10	101093	.	G	A	4.63	.	DP=13;VDB=2.845191e-02;RPB=-1.098701e+00;AF1=0.2602;AC1=6;DP4=4,6,1,2;MQ=25;FQ=5.52;PV4=1,0.07,1,0.16	GT:PL:GQ	0/0:0,6,45:11	0/0:0,0,0:5	0/0:0,3,25:8	0/0:0,9,56:13	0/1:19,0,19:15	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,0,0:5	0/1:12,0,19:8	0/1:19,0,19:15	0/0:0,3,25:8
#   10	101106	.	C	A	45.6	.	DP=17;VDB=1.025171e-02;RPB=1.317616e+00;AF1=0.3538;AC1=8;DP4=5,7,2,3;MQ=25;FQ=47.5;PV4=1,1,1,1	GT:PL:GQ	0/0:0,6,45:7	0/1:0,0,0:3	0/1:19,0,19:18	0/0:0,12,82:12	0/1:25,3,0:5	0/0:0,3,25:5	0/1:0,0,0:3	0/1:0,0,0:3	0/1:38,0,38:37	0/1:16,0,37:16	0/1:0,0,0:3
#   10	102534	.	C	T	74.1	.	DP=17;VDB=8.593325e-02;RPB=1.758816e+00;AF1=0.1922;AC1=4;DP4=4,7,2,4;MQ=25;FQ=75;PV4=1,0.25,1,1	GT:PL:GQ	0/0:0,6,45:10	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,0,0:5	0/0:0,3,25:7	0/0:0,6,45:10	0/0:0,3,25:7	0/1:59,0,6:16	0/0:0,9,71:13	0/1:59,0,13:22	0/0:0,0,0:5
#   10	103545	.	G	A	96.6	.	DP=8;VDB=3.564564e-02;RPB=1.114905e+00;AF1=0.808;AC1=18;DP4=1,1,1,5;MQ=25;FQ=11.2;PV4=0.46,0.0054,1,1	GT:PL:GQ	1/1:23,3,0:5	1/1:0,0,0:3	1/1:0,0,0:3	1/1:0,0,0:3	1/1:25,3,0:5	0/1:0,6,50:4	1/1:0,0,0:3	1/1:0,0,0:3	1/1:50,6,0:8	1/1:45,6,0:8	1/1:0,0,0:3
#   10	103556	.	T	A	118	.	DP=9;VDB=1.616396e-02;RPB=1.694823e+00;AF1=0.8551;AC1=19;DP4=0,2,2,5;MQ=25;FQ=4.33;PV4=1,0.26,1,1	GT:PL:GQ	1/1:25,3,0:6	1/1:0,0,0:4	1/1:0,0,0:4	1/1:0,0,0:4	1/1:25,3,0:6	0/1:16,0,37:24	1/1:0,0,0:4	1/1:0,0,0:4	1/1:50,6,0:9	1/1:45,6,0:9	1/1:0,0,0:4
#   10	107710	.	G	A	19.8	.	DP=10;VDB=2.309944e-02;RPB=7.878529e-01;AF1=0.2655;AC1=6;DP4=2,5,2,1;MQ=25;FQ=21.1;PV4=0.5,0.24,1,0.37	GT:PL:GQ	0/0:0,3,25:6	0/0:0,0,0:4	0/0:0,3,25:6	0/0:0,0,0:4	0/0:0,0,0:4	0/0:0,3,25:6	0/1:50,6,0:4	0/1:16,0,33:14	0/0:0,3,25:6	0/0:0,0,0:4	0/0:0,3,25:6
#   10	107938	.	A	G	67.3	.	DP=8;VDB=7.681113e-02;RPB=1.551181e+00;AF1=0.8152;AC1=18;DP4=2,1,3,2;MQ=25;FQ=15.7;PV4=1,1,1,1	GT:PL:GQ	1/1:0,0,0:3	1/1:0,0,0:3	1/1:0,0,0:3	1/1:45,6,0:8	1/1:0,0,0:3	1/1:0,0,0:3	1/1:25,3,0:5	1/1:0,0,0:3	0/1:13,0,49:20	1/1:0,0,0:3	1/1:25,3,0:5
#   10	107940	.	T	G	65.3	.	DP=8;VDB=7.777030e-02;RPB=1.551181e+00;AF1=0.8151;AC1=18;DP4=2,1,3,2;MQ=25;FQ=25.6;PV4=1,0.14,1,1	GT:PL:GQ	1/1:0,0,0:3	1/1:0,0,0:3	1/1:0,0,0:3	1/1:43,6,0:8	1/1:0,0,0:3	1/1:0,0,0:3	1/1:25,3,0:5	1/1:0,0,0:3	0/1:13,0,59:20	1/1:0,0,0:3	1/1:25,3,0:5
#   10	111174	.	A	G	18.9	.	DP=19;VDB=4.644186e-02;RPB=-1.677051e-01;AF1=0.09574;AC1=2;DP4=4,12,1,2;MQ=25;FQ=19.3;PV4=1,0.24,1,1	GT:PL:GQ	0/0:0,0,0:8	0/0:0,3,25:10	0/0:0,0,0:8	0/0:0,12,86:19	0/1:56,0,35:45	0/0:0,3,25:10	0/0:0,0,0:8	0/0:0,9,62:16	0/0:0,12,75:19	0/0:0,3,25:10	0/0:0,0,0:8
#   10	119590	.	A	C	5.8	.	DP=8;VDB=6.240000e-02;RPB=1.621514e+00;AF1=0.1742;AC1=4;DP4=5,1,2,0;MQ=25;FQ=6.42;PV4=1,0.34,1,1	GT:PL:GQ	0/0:0,6,45:12	0/0:0,0,0:7	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,0,0:7	0/0:0,0,0:7	0/0:0,0,0:7	0/1:45,6,0:7	0/0:0,0,0:7
#   10	119893	.	G	T	10.5	.	DP=100;VDB=7.681113e-02;RPB=-8.431514e-01;AF1=0.1642;AC1=4;DP4=39,49,2,3;MQ=25;FQ=11.2;PV4=1,0.21,1,1	GT:PL:GQ	0/0:4,0,103:3	0/1:25,3,0:8	0/0:0,24,152:29	0/0:0,51,217:56	0/0:0,18,123:23	0/0:0,39,202:44	0/0:0,33,183:38	0/0:0,24,131:29	0/0:0,33,183:38	0/1:23,0,116:18	0/1:19,0,19:14
#   10	119894	.	A	T	10.5	.	DP=100;VDB=7.681113e-02;RPB=-6.046844e-01;AF1=0.1642;AC1=4;DP4=40,48,2,3;MQ=25;FQ=11.2;PV4=1,0.2,1,1	GT:PL:GQ	0/0:4,0,101:3	0/1:25,3,0:8	0/0:0,24,152:29	0/0:0,48,220:53	0/0:0,18,123:23	0/0:0,39,202:44	0/0:0,33,183:38	0/0:0,24,126:29	0/0:0,36,188:41	0/1:23,0,111:18	0/1:19,0,19:14
#   10	120211	.	C	G	30.6	.	DP=96;VDB=1.009663e-01;RPB=1.188167e+00;AF1=0.0919;AC1=2;DP4=33,57,2,4;MQ=25;FQ=31;PV4=1,0.11,1,1	GT:PL:GQ	0/0:0,24,147:31	0/0:0,15,109:22	0/0:0,30,159:37	0/0:0,42,202:49	0/0:0,24,152:31	0/0:0,11,156:18	0/0:0,24,145:31	0/1:52,0,53:45	0/0:0,36,194:43	0/1:29,0,89:22	0/0:0,18,125:25
#   10	122613	.	C	G	33.2	.	DP=19;VDB=1.725012e-02;RPB=1.465882e+00;AF1=0.3002;AC1=6;DP4=4,5,2,2;MQ=25;FQ=34.7;PV4=1,1,1,1	GT:PL:GQ	0/1:50,6,0:4	0/0:0,0,0:3	0/0:0,3,25:5	0/0:0,6,45:8	0/0:0,3,25:5	0/0:0,3,25:5	0/1:16,0,37:15	0/0:0,3,25:5	0/0:0,0,0:3	0/1:19,0,19:17	0/0:0,0,0:3
#   10	126070	.	C	T	51	.	DP=18;VDB=5.770788e-02;RPB=-9.857281e-02;AF1=0.3829;AC1=8;DP4=6,7,2,3;MQ=25;FQ=53.1;PV4=1,1,1,0.42	GT:PL:GQ	0/0:0,9,71:9	0/1:0,0,0:3	0/1:25,3,0:4	0/1:25,3,0:4	0/0:0,6,50:7	0/1:0,0,0:3	0/1:10,0,77:11	1/1:50,6,0:3	0/0:0,9,70:9	0/0:0,3,25:4	0/1:0,0,0:3
#   10	126402	.	T	C	81.1	.	DP=13;VDB=1.045862e-01;RPB=-8.525197e-01;AF1=0.4945;AC1=11;DP4=2,5,1,5;MQ=25;FQ=84;PV4=1,0.48,1,1	GT:PL:GQ	0/1:16,0,41:18	0/1:0,0,0:3	0/1:0,0,0:3	0/0:0,6,45:5	0/1:0,0,0:3	0/0:0,6,50:5	0/1:25,3,0:3	0/1:0,0,0:3	0/0:0,3,23:3	0/1:25,3,0:3	1/1:71,9,0:6
#   10	127636	.	G	C	21.4	.	DP=16;VDB=6.621023e-02;RPB=-5.381382e-01;AF1=0.08591;AC1=1;DP4=7,6,1,2;MQ=25;FQ=21.8;PV4=1,0.29,1,1	GT:PL:GQ	0/0:0,3,25:11	0/0:0,3,25:11	0/0:0,3,25:11	0/0:0,9,65:17	0/0:0,3,25:11	0/0:0,0,0:8	0/1:59,0,13:26	0/0:0,0,0:8	0/0:0,3,25:11	0/0:0,6,50:14	0/0:0,6,50:14
#   10	127924	.	A	G	14.4	.	DP=12;VDB=6.186179e-02;RPB=3.698001e-01;AF1=0.3063;AC1=7;DP4=7,2,2,1;MQ=25;FQ=15.9;PV4=1,0.051,1,1	GT:PL:GQ	0/0:0,0,0:3	0/0:0,0,0:3	0/0:0,3,25:6	0/0:0,3,25:6	0/1:25,3,0:6	0/0:0,9,70:11	0/1:16,0,41:15	0/0:0,3,25:6	0/0:0,3,25:6	0/1:25,3,0:6	0/0:0,0,0:3
#   10	129508	.	C	T	3.07	.	DP=20;VDB=6.621023e-02;RPB=-6.149187e-01;AF1=0.2002;AC1=4;DP4=10,6,2,1;MQ=25;FQ=3.58;PV4=1,0.13,1,1	GT:PL:GQ	0/1:19,0,19:13	0/0:0,9,70:16	0/0:0,9,62:16	0/0:0,3,25:10	0/0:0,3,25:10	0/0:0,3,25:10	0/0:0,3,25:10	0/1:16,0,41:10	0/0:0,0,0:7	0/1:16,0,41:10	0/0:0,3,25:10
#   10	129509	.	G	A	3.58	.	DP=19;VDB=6.621023e-02;RPB=-4.738791e-01;AF1=0.2161;AC1=5;DP4=9,6,2,1;MQ=25;FQ=4.19;PV4=1,0.043,1,1	GT:PL:GQ	0/1:19,0,19:13	0/0:0,9,70:15	0/0:0,9,62:15	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,3,25:9	0/1:16,0,41:11	0/0:0,0,0:7	0/1:16,0,41:11	0/0:0,0,0:7
#   10	131337	.	A	G	999	.	DP=18;VDB=7.560302e-02;AF1=1;AC1=22;DP4=0,0,6,9;MQ=25;FQ=-31.8	GT:PL:GQ	1/1:0,0,0:11	1/1:45,6,0:17	1/1:25,3,0:14	1/1:0,0,0:11	1/1:45,6,0:17	1/1:25,3,0:14	1/1:71,9,0:20	1/1:50,6,0:17	1/1:25,3,0:14	1/1:45,6,0:17	1/1:25,3,0:14
#   10	131636	.	T	C	50.5	.	DP=15;VDB=1.045599e-01;RPB=-5.892557e-02;AF1=0.4018;AC1=9;DP4=7,2,1,5;MQ=25;FQ=52.8;PV4=0.041,0.13,1,1	GT:PL:GQ	0/1:29,0,16:21	0/1:19,0,19:19	0/1:0,0,0:3	0/0:0,6,50:6	0/0:0,3,21:4	0/0:0,3,25:4	0/1:16,0,37:17	0/1:25,3,0:4	0/0:0,3,25:4	0/1:19,3,0:4	0/1:0,0,0:3
#   10	133273	.	G	C	13.2	.	DP=19;VDB=6.086650e-02;RPB=1.066228e+00;AF1=0.1583;AC1=3;DP4=10,5,1,2;MQ=25;FQ=13.9;PV4=0.53,1,1,1	GT:PL:GQ	0/0:0,9,62:14	0/0:0,3,25:9	0/1:13,0,59:9	0/0:0,3,25:9	0/0:0,6,50:11	0/0:0,0,0:6	0/0:0,3,25:9	0/1:50,6,0:6	0/0:0,3,25:9	0/0:0,3,25:9	0/0:0,6,50:11
#   10	133402	.	G	C	3.01	.	DP=33;VDB=6.080000e-02;RPB=-6.413192e-01;AF1=0.1828;AC1=4;DP4=15,16,2,0;MQ=25;FQ=3.47;PV4=0.48,0.36,1,1	GT:PL:GQ	0/0:0,0,0:8	0/1:25,3,0:9	0/0:0,54,232:61	0/0:0,15,109:22	0/0:0,0,0:8	0/0:0,9,71:16	0/0:0,6,50:13	0/0:0,9,62:16	0/0:0,0,0:8	0/1:25,3,0:9	0/0:0,0,0:8
#   10	134767	.	A	G	19.8	.	DP=12;VDB=6.442004e-02;RPB=5.547002e-01;AF1=0.2017;AC1=4;DP4=3,6,2,1;MQ=25;FQ=20.7;PV4=0.52,0.16,1,1	GT:PL:GQ	0/0:0,3,25:7	0/0:0,3,25:7	0/0:0,3,25:7	0/0:0,0,0:5	0/0:0,6,43:10	0/0:0,0,0:5	0/1:41,0,16:25	0/0:0,3,25:7	0/0:0,6,50:10	0/0:0,0,0:5	0/1:25,3,0:7
#   10	135237	.	A	G	999	.	DP=12;VDB=5.934750e-02;RPB=1.423025e+00;AF1=0.8939;AC1=20;DP4=0,1,3,7;MQ=25;FQ=-14.2;PV4=1,1,1,1	GT:PL:GQ	1/1:0,0,0:5	1/1:0,0,0:5	1/1:25,3,0:8	1/1:0,0,0:5	1/1:25,3,0:8	1/1:50,6,0:11	1/1:25,3,0:8	1/1:25,3,0:8	0/1:0,3,25:8	1/1:70,9,0:13	1/1:25,3,0:8
#   10	135656	.	T	G	3.41	.	DP=12;VDB=5.760000e-02;RPB=-9.139077e-01;AF1=0.1184;AC1=2;DP4=3,5,2,0;MQ=25;FQ=3.72;PV4=0.44,0.31,1,1	GT:PL:GQ	0/0:0,3,25:12	0/0:0,0,0:9	0/0:0,3,25:12	0/0:0,3,25:12	0/0:0,0,0:9	0/0:0,0,0:9	0/0:0,6,43:15	0/0:0,3,25:12	0/0:0,0,0:9	0/0:0,3,25:12	0/1:37,0,16:26
#   10	135708	.	G	A	90.3	.	DP=18;VDB=1.167705e-01;RPB=-5.775402e-01;AF1=0.3842;G3=4.252e-05,1,9.799e-25;HWE=0.0132;AC1=8;DP4=6,4,7,1;MQ=25;FQ=92.4;PV4=0.31,1,1,1	GT:PL:GQ	0/0:0,3,25:4	0/1:0,0,0:3	0/1:37,0,16:21	0/0:0,3,25:4	0/0:0,3,24:4	0/1:0,0,0:3	0/1:16,0,37:17	0/1:37,0,16:21	0/1:0,0,0:3	0/1:19,0,19:18	0/1:34,0,38:34
#   10	137306	.	G	C	3.48	.	DP=16;VDB=4.160000e-02;RPB=-4.246039e-01;AF1=0.09683;AC1=2;DP4=4,9,1,1;MQ=25;FQ=3.73;PV4=1,1,1,1	GT:PL:GQ	0/0:0,3,25:13	0/0:0,0,0:10	0/0:0,0,0:10	0/0:0,3,25:13	0/0:0,3,25:13	0/0:0,6,45:16	0/1:38,0,34:28	0/0:0,0,0:10	0/0:0,3,25:13	0/0:0,3,25:13	0/0:0,12,91:22
#    
#    
#    GT:GT:PL:DP:SP:GQ:GQ	
#    0/0:0,0,0:0:0:9:0/0:9	
#    0/0:0,3,25:1:0:12:0/0:12	
#    0/0:0,0,0:0:0:9:0/0:9	
#    0/0:0,0,0:0:0:9:0/0:9	
#    0/0:0,3,25:1:0:12:0/0:12	
#    0/0:0,0,0:0:0:9:0/0:9	
#    0/0:0,3,25:1:0:12:0/0:12	
#    0/0:0,0,0:0:0:9:0/0:9	
#    0/0:0,0,0:0:0:9:0/0:9	
#    
#    GT	GT	PL	DP	SP	GQ	GQ	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,3,25	1	0	12	0/0	12	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,3,25	1	0	12	0/0	12	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,3,25	1	0	12	0/0	12	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,6,45	2	0	15	0/0	15	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,3,25	1	0	12	0/0	12	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,12,75	4	0	21	0/0	21	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,0,0	0	0	9	0/0	9	
#    0/0	0,3,25	1	0	12	0/0	12	
