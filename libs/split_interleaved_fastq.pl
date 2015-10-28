#!/usr/bin/perl-w
 
use strict;
 
use warnings;
 
use Getopt::Long;
 
use Pod::Usage;
 
use File::Basename;
 
# Date: 14-05-2010
 
# This program takes a fastq file containing paired-end reads in interleaved format as input and returns two separate files containing read 1 and read 2 in the correct order.
 
# Author: Ram Vinay Pandey 
 

 
# Define the variables
 
my $input="";
 
my $output="";
 
my $help=0;
 
my $test=0;
 
my $verbose=1;
 

 
my $usage="perl $0 --input interleaved_fastq_file.fastq --output output.fastq\n";
 

 
GetOptions(
 
    "input=s"       =>\$input,
 
    "output=s"      =>\$output,
 
    "test"          =>\$test,
 
    "help"          =>\$help
 
) or pod2usage(-msg=>"Wrong options",-verbose=>1);
 

 
pod2usage(-verbose=>2) if $help;
 
Test::runTests() if $test;
 

 
pod2usage(-msg=>"\n\tProvide an input file!!\n\n\t\t$usage\n\n",-verbose=>1) unless -e $input;
 
pod2usage(-msg=>"\n\tProvide an output file!!\n\n\t\t$usage\n\n",-verbose=>1) unless $output;
 

 

 
my ( $name, $path, $extension ) = File::Basename::fileparse ( $output, '\..*' );
 
my $output1 = $name."-read1.fastq";
 
my $output2 = $name."-read2.fastq";
 

 
open my $ofh1, ">$output1" or die "Could not open output file";
 
open my $ofh2, ">$output2" or die "Could not open output file";
 

 

 
open (IN, "<$input") or die ("Could not open file $input for reading\n");
 

 
while (<IN>) {
 
    chomp;
 
    s/\r/\n/g;
 
    # discard blank line
 
    if (m/^\s*$/g) {
 
        next;
 
    }
 
    else {
 
        # Reading all lines for read 1
 
        if (m/^\s*\@.*1$/) {
 

 
            print $ofh1 "$_\n";
 
            my $ct=0;
 
            while(my $l = <IN>) {
 
                
 
                $ct++;
 
                chomp $l;
 
                s/\r/\n/g;
 
                print $ofh1 "$l\n";
 
                
 
                last if($ct ==3);
 
                
 
            }
 
            
 
        }
 
        # Reading all lines for read 2
 
        if (m/^\s*\@.*2$/) {
 

 
            print $ofh2 "$_\n";
 
            my $ct=0;
 
            while(my $l = <IN>) {
 
                
 
                $ct++;
 
                chomp $l;
 
                s/\r/\n/g;
 
                print $ofh2 "$l\n";
 
                
 
                last if($ct ==3);
 
            }
 
            
 
            
 
        }
 
        
 
    }
 
}
 

 
close IN;
 
close $ofh1;
 
close $ofh2;
 

 

 

 
=head1 NAME
 

 
split-interleaved-fastq.pl - TThis program takes a fastq file containing paired-end reads in interleaved format as input and returns two separate files containing read 1 and read 2 in the correct order. 
 

 
=head1 SYNOPSIS
 

 
 perl split-interleaved-fastq.pl --input interleaved_fastq_file.fastq --output output.fastq
 

 
=head1 OPTIONS
 

 
=over 4
 

 
=item B<--input>
 

 
The input file which contains read1 and read2 in a single file in FASTQ format. Mandatory parameter
 

 
=item B<--output>
 

 
The output file. Mandatory parameter
 

 
=item B<--help>
 

 
Display help for this script
 

 
=back
 

 
=head1 Details
 

 
=head2 Input
 

 
The paired-end reads in interleaved format; input file looks like following:
 

 
 @fc_HWUSI-EAS613R:1:1:4:682#CATA/1
 
 TTGTANGATTTCGTCCAGACTTATCTGGAGCATCCGGACGGTCGGGTGAAGCTCAATCCTCAGCTGGTGTTG
 
 +fc_HWUSI-EAS613R:1:1:4:682#CATA/1
 
 baaa\DVbbbbaaaaa`[`abaaaab`b`]aab]_aaa``Z^`a[SN^QR^`]___aXK[a\T\[UTWWMZV
 
 @fc_HWUSI-EAS613R:1:1:4:682#CATA/2
 
 TCGACAGCTGCTGCTCCGTATTGAGGTACGGATCGTTCACGATCATATACGCCCTCTCTTTCAAAAACCTCA
 
 +fc_HWUSI-EAS613R:1:1:4:682#CATA/2
 
 bbbY`[`a\S_Y][XPaUDZ__LLL]TZPWXBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
 
 
 
=head2 Output
 

 
The output looks like as following for reads 1:
 

 
 @fc_HWUSI-EAS613R:1:1:4:682#CATA/1
 
 TTGTANGATTTCGTCCAGACTTATCTGGAGCATCCGGACGGTCGGGTGAAGCTCAATCCTCAGCTGGTGTTG
 
 +fc_HWUSI-EAS613R:1:1:4:682#CATA/1
 
 baaa\DVbbbbaaaaa`[`abaaaab`b`]aab]_aaa``Z^`a[SN^QR^`]___aXK[a\T\[UTWWMZV
 

 
The output looks like as following for reads 2:
 

 
 @fc_HWUSI-EAS613R:1:1:4:682#CATA/2
 
 TCGACAGCTGCTGCTCCGTATTGAGGTACGGATCGTTCACGATCATATACGCCCTCTCTTTCAAAAACCTCA
 
 +fc_HWUSI-EAS613R:1:1:4:682#CATA/2
 
 bbbY`[`a\S_Y][XPaUDZ__LLL]TZPWXBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
 

 

 
=head1 AUTHORS
 

 
Ram vinay pandey
 

 
=cut 

 