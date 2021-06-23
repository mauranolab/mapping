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
        my %code2base = ( 0 => $ref, '.' => 'N' );
        #Add multiallelic alleles
        my @altalleles = split(",", $alt);
#        warn "\naltalleles $ref / $alt: @altalleles\n";
        for (my $i = 0; $i < @altalleles; $i++) {
            $code2base{$i+1} = $altalleles[$i];
        }
#        warn "code2base $ref / $alt: $_ $code2base{$_}\n" for keys %code2base;
        
        #set bed coords
        my $min0 = $pos - 1;
        my $max1 = $pos - 1 + length($ref);
        
        my $totalDP = undef;
        my $numRef = "NA";
        my $numNonRef = "NA";
        my $PE = undef;
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
            } elsif ("END" eq $key) {
                # structural variant coords (i.e. Delly) require overriding the prior coordinates using END
                $max1 = int($value) - 1;
            } elsif ("PE" eq $key) {
                $PE = int($value);
            }
        }
        if ($min_total_depth_threshold > 0 and ($totalDP ne undef and $totalDP < $min_total_depth_threshold)) {
            #warn "Skipping $chrom $pos with total depth $totalDP < $min_total_depth_threshold\n";
            next;
        }
        if ($min_allele_depth_threshold > 0 and (($numRef ne "NA" and $numRef < $min_allele_depth_threshold) or ($numNonRef ne "NA" and $numNonRef < $min_allele_depth_threshold))) {
            #warn "Skipping $chrom $pos with ref/nonref depth $numRef/$numNonRef < $min_allele_depth_threshold\n";
            next;
        }
        
        my $bed3 = "$chrom\t$min0\t$max1";
        for (my $i = 0; $i < $NUM_SAMPLES; ++$i) {
            my $sampleName = $sampleNames[$i];
            my ($betterGT, $GQ, $DP, $PL) = parseGenotype( $format, $genotypes[$i], \%code2base );
            next if ("N/N" eq $betterGT); # parsed from "./."
            next if ($DP ne "NA" and $DP < $min_total_depth_threshold); # Require minimum per-sample depth of 8
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
    #TODO generalize this rather than hardcoding format strings
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
    } elsif ($format eq "GT:PL:DP:AD") {
        my ($GT1,$PL1,$DP,$AD,$GQ) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = "NA";
        $GT = $GT1;
        $PL = $PL1;
    } elsif ($format eq "GT:AD:DP:GQ:PP") {
        my ($GT1,$AD,$DP,$GQ,$PP) = split /\:/, $genotype;
        $readDepth = $DP;
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = "NA";
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
    } elsif ($format eq "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV") {
        #Format for DELLY
        my ($GT1,$GL,$GQ,$FT,$RCL,$RC,$RCR,$RDCN,$DR,$DV,$RR,$RV) = split /\:/, $genotype;
        $readDepth = "NA";
        $qualityScore = $GQ;
        $GT = $GT1;
        $PL = "NA";
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
