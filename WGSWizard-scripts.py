
#Get the original sequence for comparison

from Bio import Entrez, SeqIO

results = []
Entrez.email = email
cds = "1280421"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=cds) as handle:
    seq_record = SeqIO.read(handle, "gb")
for feature in seq_record.features:
    if feature.type == "CDS":
        dmel_cds = feature.location.extract(seq_record).seq
results.append({'DNA': str(dmel_cds), 'Protein': str(dmel_cds.translate(to_stop=True))})


#BLAST against other species and find complete sequences

import WGSWizard as wiz

dmel_seq = "442632573"  #D melanogaster  Commissureless
species_dict = ["txid7244[orgn]", "txid7260[orgn]", "txid7291[orgn]"]  #D virilis, D willistoni, D albomicans
email = "maya.gosztyla@gmail.com"
path = "/Users/mayagosztyla/taxid2wgs.pl"

for species in species_list:
    wiz.getdb(species=species, filepath=path)
#Wait for them to pop up in PyCharm! Takes a few seconds.

for species in species_list:
    print("Blasting " + species)
    blast_record = wiz.wgsblast(query_seq=dmel_seq, species=species)
    seq = wiz.wgsseq(email=email, blast=blast_record)
    results.append({'DNA': str(seq), 'Protein': str(seq.translate(to_stop=True))})
    print("Finished " + species)


#Generate multiple-sequence alignment

import pandas as pd
from Bio.Align.Applications import MuscleCommandline

results_df = pd.DataFrame(results, index=["Dmelanogaster", "Dvirilis", "Dwillisoni", "Dalbomicans"])
results_df['Protein'].to_csv("My BLAST results.tab", sep="\t", header=False)
SeqIO.convert("My BLAST results.tab", "tab", "My BLAST results.fasta", "fasta")

muscle = r"/Users/mayagosztyla/muscle3.8.31_i86darwin64"
cmdline = MuscleCommandline(muscle, input= "My BLAST results.fasta", out= "My BLAST results.aln", clw=True)
cmdline()


#Generate phylogenetic tree

from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline

phyml = r"/Users/mayagosztyla/PhyML-3.1_macOS-MountainLion"
AlignIO.convert("My BLAST results.aln", "clustal", "My BLAST results.phy", "phylip-relaxed")
cmdline2 = PhymlCommandline(phyml, input='My BLAST Results.phy', datatype='aa', model='WAG', alpha='e', bootstrap=100)
out_log, err_log = cmdline2()

my_tree = Phylo.read("My BLAST results.phy_phyml_tree.txt", "newick")
Phylo.draw_ascii(my_tree)
Phylo.write(my_tree, "My BLAST results tree.xml", 'phyloxml')