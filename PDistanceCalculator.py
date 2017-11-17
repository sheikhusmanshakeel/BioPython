# import required Biophython functions
from Bio import Entrez
from Bio import SeqIO

humanSequence = "NP_000042.3"
sequence1Reference = "Q6PQD5.2"  # pig
sequence2Reference = "Q62388.2"  # house mouse
sequence3Reference = "Q8WN22.1"  # dog
sequence4Reference = "P32871.1"  # cattle

Entrez.email = 'A.N.Other@example.com'
handle = Entrez.efetch(db="nucleotide", id=humanSequence, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
sequence1 = record._seq;

handle = Entrez.efetch(db="nucleotide", id=sequence1Reference, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
sequence2 = record._seq;

handle = Entrez.efetch(db="nucleotide", id=sequence2Reference, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
sequence3 = record._seq;

handle = Entrez.efetch(db="nucleotide", id=sequence3Reference, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
sequence4 = record._seq;

handle = Entrez.efetch(db="nucleotide", id=sequence4Reference, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
sequence5 = record._seq;

handle.close()

print(sequence1)  # human sequence
print(sequence2)
print(sequence3)
print(sequence4)
print(sequence5)

lengthOfSequence1 = len(sequence1)
lengthOfSequence2 = len(sequence2)
lengthOfSequence3 = len(sequence3)
lengthOfSequence4 = len(sequence4)
lengthOfSequence5 = len(sequence5)

# comparing sequence 1 with sequence 2
distance12 = 0;
distance13 = 0;
distance14 = 0;
distance15 = 0;

for i in range(0, lengthOfSequence1):
    if sequence1[i] != sequence2[i]:
        distance12 += 1
print("Distance 1-2 :",distance12)

for i in range(0, lengthOfSequence1):
    if sequence1[i] != sequence3[i]:
        distance13 += 1
print("Distance 1-3 :", distance13)

for i in range(0, lengthOfSequence1):
    if sequence1[i] != sequence4[i]:
        distance14 += 1
print("Distance 1-4 :",distance14)

for i in range(0, lengthOfSequence1):
    if i < lengthOfSequence5:
        if sequence1[i] != sequence5[i]:
            distance15 += 1
    else:
        distance15 +=1
print("Distance 1-5 :", distance15)


distance23 = 0
distance24 = 0;
distance25 = 0;


for i in range(0, lengthOfSequence2):
    if sequence2[i] != sequence3[i]:
        distance23 += 1
print("Distance 2-3 :", distance23)

for i in range(0, lengthOfSequence2):
    if sequence2[i] != sequence4[i]:
        distance24 += 1
print("Distance 2-4 :", distance24)

for i in range(0, lengthOfSequence2):
    if i < lengthOfSequence5:
        if sequence2[i] != sequence5[i]:
            distance25 += 1
    else:
        distance25 +=1
print("Distance 2-5 :", distance25)


distance34 = 0
distance35 = 0

for i in range(0, lengthOfSequence3):
    if sequence3[i] != sequence4[i]:
        distance34 += 1
print("Distance 3-4 :", distance34)

for i in range(0, lengthOfSequence3):
    if i < lengthOfSequence5:
        if sequence3[i] != sequence5[i]:
            distance35 += 1
    else:
        distance35 +=1
print("Distance 3-5 :", distance35)

distance45 = 0
for i in range(0, lengthOfSequence4):
    if i < lengthOfSequence5:
        if sequence4[i] != sequence5[i]:
            distance45 += 1
    else:
        distance45 +=1
print("Distance 4-5 :", distance45)

