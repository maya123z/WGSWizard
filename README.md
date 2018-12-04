# WGS Wizard

WGS Wizard allows you to perform a BLAST search against a WGS database and use the results to find the complete DNA and protein sequence for your hits from the original scaffold.


## Getting Started

###Prerequisites

**Required Python packages:** Bio, os, re

**Required external programs:** taxid2wgs, sra-blastn (see below)

In order to convert species taxonomic ID's to the names of their WGS databases, you must download the taxid2wgs.pl script from <ftp://ftp.ncbi.nlm.nih.gov/blast/WGS_TOOLS>.
The taxid2wgs script also requires that you have the sra-blastn program installed: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software.

**Optional extras:** To follow the included script and generate multiple sequence alignments, as well as phylogenetic trees, you will need to install two command line tools.

MUSCLE: https://www.drive5.com/muscle/

PhyML: http://www.atgc-montpellier.fr/phyml/

###Installation

You can install WGS Wizard directly from Github.

##Usage

###Defining the WGS Databases: getdb

The `getdb` function converts the taxonomic ID for your species of interest into the names of its WGS databases using the taxid2wgs.pl script. (See Required external programs).

After inputting the taxonomic ID and the path to where you have taxid2wgs.pl installed, the function creates a text file in the current directory that contains the database names.

You only need to do this once per species.

###Performing the WGS BLAST: wgsblast

The `wgsblast` function performs a tBLASTn search against the NCBI WGS Database based on the GI number of a protein sequence. Note that your query species must be formatted as follows: "txid7244[orgn]"

This step often takes 30 seconds or more, depending on the database.

###Finding the full gene sewuence: wgsseq

The `wgsseq` function is the meat of this module! Let's take a look at the steps that this function takes in order to track down the full DNA and protein sequence for a BLAST query.

####Step 1: Get the scaffold sequence.

Searches Entrez for the scaffold using its accession ID.

####Step 2: Put the exons in the proper order.

Uses the starting and ending points to determine which exon goes first. (Note: This program currently works with only two exons. See Limitations.)

####Step 3: Find the first in-frame upstream stop codon.

Searches the region upstream of the hit start site for the nearest in-frame stop codon.

####Step 4: Find the possible starting methionine.

Identifies the first in-frame ATG codon after the upstream stop.
Note that this may result in extra nucleotides being included if the real start codon is farther downstream.
However, this method ensures that none of the sequence is missed.

####Step 5: Find the exon 1 splice site.

Searches for a 'GTAAGT' or 'GTGAGT' sequence downstream of the exon 1 hit end, indicating a splice site, and determines where the exon 1 sequence ends.

####Step 6: Find the exon 2 splice site.

Searches for a (C/T)AG/(A/G) sequence upstream of the exon 2 hit beginning. indicating a splice site, and joins exons 1 and 2 together accordingly.

####Step 7: Find the stop codon.

Searches the region downstream of exon 2 for the first in-frame stop codon.

####Step 8: Return the complete gene sequence.

The final sequence of the gene is returned as a Bio.Seq.Seq object.

##Current Limitations

WGS Wizard is still a work in progress and has many areas for improvement! Currently, in order for this module to work, it makes the following assumptions:

- The first BLAST result is the correct one.
- The gene being BLASTed contains exactly two exons.
- The species being BLASTed against also has a homolog with exactly two exons.
- The two exons fall on the same scaffold and appear in the same BLAST result.
- The query and subject species are all Drosophila (and use the typical Drosophila splice site motifs).

Future versions of this program will address these issues to make WGS Wizard useful across many more situations!
