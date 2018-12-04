#NCBI tBLASTn from WGS database

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
import os
import re

def getdb(species: str, filepath: str) -> None:
    """
    Accesses a perl script to find the database names for a given species within the WGS database. See ftp://ftp.ncbi.nlm.nih.gov/blast/WGS_TOOLS
    :param query_species: The taxonomic ID for the species to BLAST against
    :param filepath: Path to the perl file
    :return: None. Creates a file containing a list of WGS databases.
    """
    id = re.findall('\d+', species)[0]
    command = "perl " + filepath + " url_api_ready " + id + " > wgs_dbs_tax_" + id
    os.popen(command)


def wgsblast(query_seq: str, species: str) -> list:
    """
    Performs a tblastn search for an inputted sequence ID against a target species from a Whole-Genome Shotgun (WGS) database.
    :param query_seq: GI number of a protein sequence.
    :param species: The taxonomic ID for the species to BLAST against. Must be formatted like this: "txid7244[orgn]".
    :return: A list of BLAST results.
    """
    id = re.findall('\d+', species)[0]
    path = "wgs_dbs_tax_" + id
    db = open(path).read()
    result_handle = NCBIWWW.qblast(program="tblastn", database=db, sequence=query_seq, entrez_query=species, filter=True,
                               expect=10.0, word_size=6, matrix_name='BLOSUM62', gapcosts="11 1")
    blast_record = list(NCBIXML.parse(result_handle))
    return blast_record

#Find the complete gene sequence
def wgsseq(email: str, blast: list):
    """
    :param email: Your email address (required for Entrez queries).
    :param blast: A list of BLAST results.
    :return: A Bio.Seq.Seq object containing the complete gene sequence.
    """
    #Get scaffold sequence
    Entrez.email = email
    accession = blast[0].alignments[0].accession
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=accession) as handle:
        seq_record = SeqIO.read(handle, "gb")

    #Put exons in proper order
    if blast[0].alignments[0].hsps[0].sbjct_start < blast[0].alignments[0].hsps[1].sbjct_start:
        exon1 = blast[0].alignments[0].hsps[0]
        exon2 = blast[0].alignments[0].hsps[1]
    else:
        exon1 = blast[0].alignments[0].hsps[1]
        exon2 = blast[0].alignments[0].hsps[0]

    #Find the first upstream in-frame stop codon
    e1upstream = seq_record.seq[0 : exon1.sbjct_start + 4]
    e1stop = e1upstream.rfind('TAA' or 'TAG' or 'TGA')
    while (exon1.sbjct_start - 1 - e1stop) % 3 != 0 and e1stop != exon1.sbjct_start:
        e1stop = e1upstream[0 : e1stop].rfind('TAA' or 'TAG' or 'TGA')

    #Find the first in-frame start codon after the stop
    e1start = e1upstream[e1stop : ].find('ATG')
    if e1start != -1:
        e1start = e1start + e1stop
        while (exon1.sbjct_start - 1 - e1start) % 3 != 0 and e1start != exon1.sbjct_start:
            e1start = e1upstream[e1stop + e1start + 1 : ].find('ATG') + e1stop + e1start + 1
    else:
        e1start = e1stop + 3

    #Find the exon1 splice site
    e1downstream = seq_record.seq[exon1.sbjct_end - 6 :]
    e1splice = e1downstream.find('GTAAGT' or 'GTGAGT')
    e1splice = e1splice + exon1.sbjct_end - 6
    exon1_complete = seq_record.seq[e1start : e1splice]

    #Find the exon2 splice site
    e2upstream = seq_record.seq[0 : exon2.sbjct_start + 5]
    e2splice = e2upstream.rfind('CAGA' or 'CAGG' or 'TAGA' or 'TAGG') + 3
    geneseq = exon1_complete + seq_record.seq[e2splice :]

    #Find the stop codon
    length = len(geneseq.translate(to_stop = True)) * 3
    geneseq = geneseq[0 : length + 3]
    return geneseq