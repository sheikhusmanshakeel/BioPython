# do not copy these lines if you use python/ipython
# this is only used to plot in ipython notebooks
# %pylab inline
# %load_ext autoreload
# %autoreload 2

# import required Biophython functions
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SubsMat

Entrez.email = 'A.N.Other@example.com'
handle = Entrez.efetch(db="nucleotide", id="NP_000042.3", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

# show the sequence record
# here we have choosen the human Thy1 gene
print(record)

# run a BLAST via the web
result_handle = NCBIWWW.qblast('blastp', 'swissprot', record.seq)
# parse the results
blast_record = NCBIXML.read(result_handle)  # This takes some time to run

# show all matches returned with names
print('===============================================================================================')
print('===============================================================================================')
print('========================== BLAST P Allignments =======================================')
a = blast_record.alignments
for aa in a:
    expectedValue = 0
    print('Length of expected values', len(aa.hsps))
    for h in aa.hsps:
        expectedValue = h.expect
    print(aa.title.split('|')[4].split(' ')[0], expectedValue)

# print all results with e-value below this value:
E_VALUE_THRESH = 5e-100
# lower this threshold to also see sequences with poor match

# now print, for each match:

# name of alignment
# length of alignment
# e-value
# Query sequence
# Matching sequence
# Alignment info

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:75] + '...')
            print(hsp.sbjct[0:75] + '...')
            print(hsp.match[0:75] + '...')

# now work with all results with e-value below this value:
E_VALUE_THRESH = 5e-100


def get_seqrecs(alignments, threshold):
    # a little helper function to get the sequence records
    for aln in alignments:
        for hsp in aln.hsps:
            if hsp.expect < threshold:
                id = aln.title.split('|')[4].split(' ')[0]
                # id = aln.accession
                print(id)
                yield SeqRecord(Seq(hsp.sbjct), id=id)  # ,description=str(aln.title.split('|')[4]))
                break


best_seqs = get_seqrecs(blast_record.alignments, E_VALUE_THRESH)
# write out to a fasta file
SeqIO.write(best_seqs, 'family_alignment.fasta', 'fasta')

from Bio.Align.Applications import MuscleCommandline

# run Muscle MSA
cmdline = MuscleCommandline('./muscle3.8.31_i86linux64', input='family_alignment.fasta', out='family_alignment.aln',
                            clw=True)
cmdline()

alignment = AlignIO.read('family_alignment.aln', 'clustal')
print(alignment)

summary_align = AlignInfo.SummaryInfo(alignment)

# compute a consensus sequence by taking the most frequent letter
# positions below a thresold similarity are shown as 'X'
# threshold can be adjusted by adding e.g. threshold=0.5
print('Consensus sequence without gaps:')
print(summary_align.dumb_consensus())
print('Consensus sequence with gaps:')
print(summary_align.gap_consensus())

# print a Position Specific Score Matrix (PSSM)
# this shows the number of letters counted at each location
# in tyhe sequence, which is shown in vertical along the left
pssm = summary_align.pos_specific_score_matrix(summary_align.dumb_consensus(), chars_to_ignore=['X'])
print(pssm)

# How to construct a substitution matrix from the alignment
summary_align = AlignInfo.SummaryInfo(alignment)
# done here only for charged amino acids
replace_info = summary_align.replacement_dictionary(["G", "A", "V", "L", "I", "M", "P", "F",
                                                     "W", "S", "T", "N", "Q", "Y", "C"])
my_arm = SubsMat.SeqMat(replace_info)
my_lom = SubsMat.make_log_odds_matrix(my_arm)
my_lom.print_mat()

# you notice that in this example many entries are not defined
# this is because the alignment is too short and does not have all combinations

# you can use a different gene, and/or reduce the e-value threshold to include more data

# from Bio.Phylo.TreeConstruction import DistanceCalculator
from TreeConstruction import DistanceCalculator

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
print(dm)

from TreeConstruction import DistanceTreeConstructor

# here supply the keyword upgma or nj
# compare the trees you get from both methods
constructor = DistanceTreeConstructor(calculator, 'upgma')
tree = constructor.build_tree(alignment)
print(tree)

from Bio import Phylo

# now draw the tree, try out these three methods:

# Phylo.draw_ascii(tree)
# Phylo.draw_graphviz(tree)
Phylo.draw(tree)
