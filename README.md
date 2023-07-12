# Elli-Houser-PythonII-Portfolio
This is the code that I wrote for my BISC-450C, part II, Python class.


Sequence Objects Notebook
```python
# In this notebook, we took sequence objects, otherwise known as a collection or series of similar aspects, and we combined them in order to
# manipulate and perform various functions with them. By performing these specific commands on various sequences, it allows us to relate
# them to computational biology.
```


```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
my_seq[2:3]
```




    Seq('T')




```python
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N"*10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
coding_dna.translate(table = "Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table=2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table=2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None,10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAAGGATAGT"}, length = 159345973)
```

```python
seq[1000:1020]
```

    Seq(None, length=20)

```python
seq[117512690:117512700]
```

    Seq('CCTGAATGTG')

```python
seq[117512670:]
```

    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAAGGATAGT'}, length=41833303)

```python
seq = Seq("ACGT")
```

```python
undefined_seq = Seq(None, length = 10)
```

```python
seq + undefined_seq + seq
```

    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)

```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```

```python
from Bio.Seq import MutableSeq
```

```python
mutable_seq = MutableSeq(my_seq)
```

```python
mutable_seq
```
    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq[5] = "C"
```

```python
mutable_seq
```

    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq.remove("T")
```

```python
mutable_seq
```
    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')

```python
mutable_seq.reverse()
```

```python
mutable_seq
```

    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')

```python
new_seq = Seq(mutable_seq)
```

```python
new_seq
```
    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')

```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```

```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```

```python
reverse_complement(my_string)
```

    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'
    
```python
transcribe(my_string)
```

    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'

```python
back_transcribe(my_string)
```

    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'

```python
translate(my_string)
```

    'AVMGRWKGGRAAG*'


Sequence Annotation Notebook
```python
# In this notebook, we performed various commands and functions on numerous sequences and manipulated them through various annotation commands.
# By doing so, this allows us to view specific details about each sequence chosen, as well as details that would have been unknown without
# the annotation commands.
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq ("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example."
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example.



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
# Start of Video 4 beyond this point; Part 3 URL was unavailable!
```


```python
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id = "Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(record), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



Sequence Input & Output Notebook
```python
# In this notebook, we worked with various seek files in order to view the results of various commands on different types of sequences. The file in
# which was used for this notebook was called ls_orchid, and we were able to demonstrate it as a fasta file, as well as a genbank file throughout
# the commands in the notebook. By doing so, we were able to view specific details regarding annotation commands, as well as listing commands.
```


```python
# https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk
```


```python
from Bio import SeqIO
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records [-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "desired new description"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new_id desired new description
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVSEITAVTFRGPSETHLDSMVGQALFGD"
   "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
   "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
   "SSAC",
   ),
id = "gi|14150838|gb|AAK54648.1|AF376133_1",
description = "chalcone synthase [Cucumis sativus]",)
```


```python
rec2 = SeqRecord(
Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRMYLTEEILKENPSMCEMAPSLDARQ"
   "DMVVVEIPKLGKEAAVKAIKEWGQ",
   ),
    id = "gi | 13919613 | gb| AAK33142.1|",
    description = "chalcone synthase [Fragaria vesca subsp. bracteata]",
)
```


```python
rec3 = SeqRecord(
Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
   "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
   "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
   "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
   "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
   "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
   "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
   ),
    id = "gi | 13925890 | gb| AAK49457.1|",
    description = "chalcone synthase [Nocotiana tabacum]",
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.fasta", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement)
```

    Z78533.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')>
    Z78532.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')>
    Z78531.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')>
    Z78530.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')>
    Z78529.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')>
    Z78527.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')>
    Z78526.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')>
    Z78525.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')>
    Z78524.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')>
    Z78523.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')>
    Z78522.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')>
    Z78521.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')>
    Z78520.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')>
    Z78519.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')>
    Z78518.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')>
    Z78517.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')>
    Z78516.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')>
    Z78515.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')>
    Z78514.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')>
    Z78513.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')>
    Z78512.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')>
    Z78511.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')>
    Z78510.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')>
    Z78509.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')>
    Z78508.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')>
    Z78507.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')>
    Z78506.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')>
    Z78505.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')>
    Z78504.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')>
    Z78503.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')>
    Z78502.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')>
    Z78501.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')>
    Z78500.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')>
    Z78499.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')>
    Z78498.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')>
    Z78497.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78496.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78495.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')>
    Z78494.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')>
    Z78493.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')>
    Z78492.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')>
    Z78491.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78490.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78489.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')>
    Z78488.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')>
    Z78487.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78486.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')>
    Z78485.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')>
    Z78484.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78483.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78482.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')>
    Z78481.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78480.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78479.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')>
    Z78478.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')>
    Z78477.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78476.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78475.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')>
    Z78474.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')>
    Z78473.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78472.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78471.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78470.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78469.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78468.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')>
    Z78467.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78466.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78465.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78464.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')>
    Z78463.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')>
    Z78462.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')>
    Z78461.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78460.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')>
    Z78459.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78458.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')>
    Z78457.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')>
    Z78456.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78455.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')>
    Z78454.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')>
    Z78453.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78452.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78451.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78450.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')>
    Z78449.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78448.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78447.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')>
    Z78446.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78445.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')>
    Z78444.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')>
    Z78443.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78442.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')>
    Z78441.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')>
    Z78440.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')>
    Z78439.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')>



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
print(records)
```

    [SeqRecord(seq=Seq('GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGT...ACG'), id='rc_gi|2765658|emb|Z78533.1|CIZ78533', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAAT...ACG'), id='rc_gi|2765657|emb|Z78532.1|CCZ78532', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCC...ACG'), id='rc_gi|2765656|emb|Z78531.1|CFZ78531', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGG...ACG'), id='rc_gi|2765655|emb|Z78530.1|CMZ78530', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAAC...CGT'), id='rc_gi|2765654|emb|Z78529.1|CLZ78529', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGT...ACG'), id='rc_gi|2765652|emb|Z78527.1|CYZ78527', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAAT...ACG'), id='rc_gi|2765651|emb|Z78526.1|CGZ78526', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATC...ACA'), id='rc_gi|2765650|emb|Z78525.1|CAZ78525', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAAC...ACG'), id='rc_gi|2765649|emb|Z78524.1|CFZ78524', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCA...ACG'), id='rc_gi|2765648|emb|Z78523.1|CHZ78523', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCC...ACG'), id='rc_gi|2765647|emb|Z78522.1|CMZ78522', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAG...TAC'), id='rc_gi|2765646|emb|Z78521.1|CCZ78521', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACC...ACG'), id='rc_gi|2765645|emb|Z78520.1|CSZ78520', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAA...TAT'), id='rc_gi|2765644|emb|Z78519.1|CPZ78519', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATC...ACG'), id='rc_gi|2765643|emb|Z78518.1|CRZ78518', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAAT...ACG'), id='rc_gi|2765642|emb|Z78517.1|CFZ78517', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCA...ACG'), id='rc_gi|2765641|emb|Z78516.1|CPZ78516', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAAT...ACG'), id='rc_gi|2765640|emb|Z78515.1|MXZ78515', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCC...ACG'), id='rc_gi|2765639|emb|Z78514.1|PSZ78514', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGG...ACG'), id='rc_gi|2765638|emb|Z78513.1|PBZ78513', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGT...ACG'), id='rc_gi|2765637|emb|Z78512.1|PWZ78512', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGA...ACG'), id='rc_gi|2765636|emb|Z78511.1|PEZ78511', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...TAG'), id='rc_gi|2765635|emb|Z78510.1|PCZ78510', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACC...ACG'), id='rc_gi|2765634|emb|Z78509.1|PPZ78509', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGG...ACG'), id='rc_gi|2765633|emb|Z78508.1|PLZ78508', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGG...ACG'), id='rc_gi|2765632|emb|Z78507.1|PLZ78507', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCC...ACG'), id='rc_gi|2765631|emb|Z78506.1|PLZ78506', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765630|emb|Z78505.1|PSZ78505', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCA...ACG'), id='rc_gi|2765629|emb|Z78504.1|PKZ78504', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765628|emb|Z78503.1|PCZ78503', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_gi|2765627|emb|Z78502.1|PBZ78502', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGT...ACG'), id='rc_gi|2765626|emb|Z78501.1|PCZ78501', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTC...ACG'), id='rc_gi|2765625|emb|Z78500.1|PWZ78500', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGT...ACG'), id='rc_gi|2765624|emb|Z78499.1|PMZ78499', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAAC...ACG'), id='rc_gi|2765623|emb|Z78498.1|PMZ78498', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765622|emb|Z78497.1|PDZ78497', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCA...ACG'), id='rc_gi|2765621|emb|Z78496.1|PAZ78496', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGA...ACG'), id='rc_gi|2765620|emb|Z78495.1|PEZ78495', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATT...ACG'), id='rc_gi|2765619|emb|Z78494.1|PNZ78494', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATG...ACG'), id='rc_gi|2765618|emb|Z78493.1|PGZ78493', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAAT...ACG'), id='rc_gi|2765617|emb|Z78492.1|PBZ78492', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAA...ACG'), id='rc_gi|2765616|emb|Z78491.1|PCZ78491', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCA...ACG'), id='rc_gi|2765615|emb|Z78490.1|PFZ78490', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765614|emb|Z78489.1|PDZ78489', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGG...CAG'), id='rc_gi|2765613|emb|Z78488.1|PTZ78488', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATC...ACG'), id='rc_gi|2765612|emb|Z78487.1|PHZ78487', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGAT...ACG'), id='rc_gi|2765611|emb|Z78486.1|PBZ78486', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACT...CAG'), id='rc_gi|2765610|emb|Z78485.1|PHZ78485', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCA...ACG'), id='rc_gi|2765609|emb|Z78484.1|PCZ78484', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765608|emb|Z78483.1|PVZ78483', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...AGA'), id='rc_gi|2765607|emb|Z78482.1|PEZ78482', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_gi|2765606|emb|Z78481.1|PIZ78481', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGAT...ACG'), id='rc_gi|2765605|emb|Z78480.1|PGZ78480', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTG...ACG'), id='rc_gi|2765604|emb|Z78479.1|PPZ78479', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAA...ACG'), id='rc_gi|2765603|emb|Z78478.1|PVZ78478', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_gi|2765602|emb|Z78477.1|PVZ78477', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGG...ACG'), id='rc_gi|2765601|emb|Z78476.1|PGZ78476', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCC...ACG'), id='rc_gi|2765600|emb|Z78475.1|PSZ78475', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAAC...ACG'), id='rc_gi|2765599|emb|Z78474.1|PKZ78474', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAAT...ACG'), id='rc_gi|2765598|emb|Z78473.1|PSZ78473', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765597|emb|Z78472.1|PLZ78472', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765596|emb|Z78471.1|PDZ78471', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTAT...ACG'), id='rc_gi|2765595|emb|Z78470.1|PPZ78470', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTA...ACG'), id='rc_gi|2765594|emb|Z78469.1|PHZ78469', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAA...ACG'), id='rc_gi|2765593|emb|Z78468.1|PAZ78468', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATC...ACG'), id='rc_gi|2765592|emb|Z78467.1|PSZ78467', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGG...ACG'), id='rc_gi|2765591|emb|Z78466.1|PPZ78466', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCA...ACG'), id='rc_gi|2765590|emb|Z78465.1|PRZ78465', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCA...ACG'), id='rc_gi|2765589|emb|Z78464.1|PGZ78464', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTC...ACG'), id='rc_gi|2765588|emb|Z78463.1|PGZ78463', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAAT...ACG'), id='rc_gi|2765587|emb|Z78462.1|PSZ78462', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCA...ACG'), id='rc_gi|2765586|emb|Z78461.1|PWZ78461', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCA...ACG'), id='rc_gi|2765585|emb|Z78460.1|PCZ78460', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_gi|2765584|emb|Z78459.1|PDZ78459', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAAC...ACG'), id='rc_gi|2765583|emb|Z78458.1|PHZ78458', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGA...ACG'), id='rc_gi|2765582|emb|Z78457.1|PCZ78457', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAA...ACG'), id='rc_gi|2765581|emb|Z78456.1|PTZ78456', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCA...ACG'), id='rc_gi|2765580|emb|Z78455.1|PJZ78455', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGG...ACG'), id='rc_gi|2765579|emb|Z78454.1|PFZ78454', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765578|emb|Z78453.1|PSZ78453', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATC...ACG'), id='rc_gi|2765577|emb|Z78452.1|PBZ78452', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCA...ACG'), id='rc_gi|2765576|emb|Z78451.1|PHZ78451', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_gi|2765575|emb|Z78450.1|PPZ78450', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCAT...ACG'), id='rc_gi|2765574|emb|Z78449.1|PMZ78449', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTG...ACG'), id='rc_gi|2765573|emb|Z78448.1|PAZ78448', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765572|emb|Z78447.1|PVZ78447', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGG...ACG'), id='rc_gi|2765571|emb|Z78446.1|PAZ78446', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCA...ACG'), id='rc_gi|2765570|emb|Z78445.1|PUZ78445', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCT...ACG'), id='rc_gi|2765569|emb|Z78444.1|PAZ78444', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGG...ACG'), id='rc_gi|2765568|emb|Z78443.1|PLZ78443', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATC...TAC'), id='rc_gi|2765567|emb|Z78442.1|PBZ78442', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCA...TCC'), id='rc_gi|2765566|emb|Z78441.1|PSZ78441', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATC...ACG'), id='rc_gi|2765565|emb|Z78440.1|PPZ78440', name='<unknown name>', description='reverse complement', dbxrefs=[]), SeqRecord(seq=Seq('GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATA...ATG'), id='rc_gi|2765564|emb|Z78439.1|PBZ78439', name='<unknown name>', description='reverse complement', dbxrefs=[])]



```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18



Sequence Alignment Notebook
```python
# In this notebook, we compared and contrasted multiple sequences / gene alleles across various organisms. By doing so, this allows us to compare
# and contrast how closely these organisms can be related, or even allow us to view their genetic history from a different point of view.
```


```python
# https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05361_seed.sth
```


```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
>>> align1 = MultipleSeqAlignment(
    [
...       SeqRecord(Seq("ACTGCTAGCTAG"), id = "Alpha"),
...       SeqRecord(Seq("ACT-CTAGCTAG"), id = "Beta"),
...       SeqRecord(Seq("ACTGCTAGDTAG"), id = "Gamma"),
...   ]
... )

>>> align2 = MultipleSeqAlignment(
    [
...        SeqRecord(Seq("GTCAGC-AG"), id = "Delta"),
...        SeqRecord(Seq("GACAGCTAG"), id = "Epsilon"),
...        SeqRecord(Seq("GTCAGCTAG"), id = "Zeta"),
...     ]
... )

>>> align3 = MultipleSeqAlignment(
    [
...        SeqRecord(Seq("ACTAGTACAGCTG"), id = "Eta"),
...        SeqRecord(Seq("ACTAGTACAGCT-"), id = "Theta"),
...        SeqRecord(Seq("-CTACTACAGGTG"), id = "Iota"),
...     ]
... )
```


```python
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fc4bcc90390>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fc4bcc90490>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7fc4bcc90590>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7fc4bcc90390>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7fc4bcc90490>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7fc4bcc90590>]



```python
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA seq0
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA seq1
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA seq2
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq3
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq4
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq5
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA seq6



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2,6])
```

    T



```python
print(alignment[2].seq[6])
```

    T



```python
print(alignment[:,6])
```

    TTT---T



```python
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.SeqRecord import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id = "seq1"),
        SeqRecord(Seq("AAT-CTA"), id = "seq2"),
        SeqRecord(Seq("CCTACT-"), id = "seq3"),
        SeqRecord(Seq("TCTCCTC"), id = "seq4"),
     ]
)
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
substitutions = alignment.substitutions
```


```python
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
help(ClustalwCommandline)
```

    Help on class ClustalwCommandline in module Bio.Align.Applications._Clustalw:
    
    class ClustalwCommandline(Bio.Application.AbstractCommandline)
     |  ClustalwCommandline(cmd='clustalw', **kwargs)
     |  
     |  Command line wrapper for clustalw (version one or two).
     |  
     |  http://www.clustal.org/
     |  
     |  Notes
     |  -----
     |  Last checked against versions: 1.83 and 2.1
     |  
     |  References
     |  ----------
     |  Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA,
     |  McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD,
     |  Gibson TJ, Higgins DG. (2007). Clustal W and Clustal X version 2.0.
     |  Bioinformatics, 23, 2947-2948.
     |  
     |  Examples
     |  --------
     |  >>> from Bio.Align.Applications import ClustalwCommandline
     |  >>> in_file = "unaligned.fasta"
     |  >>> clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
     |  >>> print(clustalw_cline)
     |  clustalw2 -infile=unaligned.fasta
     |  
     |  You would typically run the command line with clustalw_cline() or via
     |  the Python subprocess module, as described in the Biopython tutorial.
     |  
     |  Method resolution order:
     |      ClustalwCommandline
     |      Bio.Application.AbstractCommandline
     |      builtins.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self, cmd='clustalw', **kwargs)
     |      Initialize the class.
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Bio.Application.AbstractCommandline:
     |  
     |  __call__(self, stdin=None, stdout=True, stderr=True, cwd=None, env=None)
     |      Execute command, wait for it to finish, return (stdout, stderr).
     |      
     |      Runs the command line tool and waits for it to finish. If it returns
     |      a non-zero error level, an exception is raised. Otherwise two strings
     |      are returned containing stdout and stderr.
     |      
     |      The optional stdin argument should be a string of data which will be
     |      passed to the tool as standard input.
     |      
     |      The optional stdout and stderr argument may be filenames (string),
     |      but otherwise are treated as a booleans, and control if the output
     |      should be captured as strings (True, default), or ignored by sending
     |      it to /dev/null to avoid wasting memory (False). If sent to a file
     |      or ignored, then empty string(s) are returned.
     |      
     |      The optional cwd argument is a string giving the working directory
     |      to run the command from. See Python's subprocess module documentation
     |      for more details.
     |      
     |      The optional env argument is a dictionary setting the environment
     |      variables to be used in the new process. By default the current
     |      process' environment variables are used. See Python's subprocess
     |      module documentation for more details.
     |      
     |      Default example usage::
     |      
     |          from Bio.Emboss.Applications import WaterCommandline
     |          water_cmd = WaterCommandline(gapopen=10, gapextend=0.5,
     |                                       stdout=True, auto=True,
     |                                       asequence="a.fasta", bsequence="b.fasta")
     |          print("About to run: %s" % water_cmd)
     |          std_output, err_output = water_cmd()
     |      
     |      This functionality is similar to subprocess.check_output(). In general
     |      if you require more control over running the command, use subprocess
     |      directly.
     |      
     |      When the program called returns a non-zero error level, a custom
     |      ApplicationError exception is raised. This includes any stdout and
     |      stderr strings captured as attributes of the exception object, since
     |      they may be useful for diagnosing what went wrong.
     |  
     |  __repr__(self)
     |      Return a representation of the command line object for debugging.
     |      
     |      e.g.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
     |      >>> cline.asequence = "asis:ACCCGGGCGCGGT"
     |      >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
     |      >>> cline.outfile = "temp_water.txt"
     |      >>> print(cline)
     |      water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
     |      >>> cline
     |      WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5)
     |  
     |  __setattr__(self, name, value)
     |      Set attribute name to value (PRIVATE).
     |      
     |      This code implements a workaround for a user interface issue.
     |      Without this __setattr__ attribute-based assignment of parameters
     |      will silently accept invalid parameters, leading to known instances
     |      of the user assuming that parameters for the application are set,
     |      when they are not.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5, stdout=True)
     |      >>> cline.asequence = "a.fasta"
     |      >>> cline.bsequence = "b.fasta"
     |      >>> cline.csequence = "c.fasta"
     |      Traceback (most recent call last):
     |      ...
     |      ValueError: Option name csequence was not found.
     |      >>> print(cline)
     |      water -stdout -asequence=a.fasta -bsequence=b.fasta -gapopen=10 -gapextend=0.5
     |      
     |      This workaround uses a whitelist of object attributes, and sets the
     |      object attribute list as normal, for these.  Other attributes are
     |      assumed to be parameters, and passed to the self.set_parameter method
     |      for validation and assignment.
     |  
     |  __str__(self)
     |      Make the commandline string with the currently set options.
     |      
     |      e.g.
     |      
     |      >>> from Bio.Emboss.Applications import WaterCommandline
     |      >>> cline = WaterCommandline(gapopen=10, gapextend=0.5)
     |      >>> cline.asequence = "asis:ACCCGGGCGCGGT"
     |      >>> cline.bsequence = "asis:ACCCGAGCGCGGT"
     |      >>> cline.outfile = "temp_water.txt"
     |      >>> print(cline)
     |      water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5
     |      >>> str(cline)
     |      'water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5'
     |  
     |  set_parameter(self, name, value=None)
     |      Set a commandline option for a program (OBSOLETE).
     |      
     |      Every parameter is available via a property and as a named
     |      keyword when creating the instance. Using either of these is
     |      preferred to this legacy set_parameter method which is now
     |      OBSOLETE, and likely to be DEPRECATED and later REMOVED in
     |      future releases.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Bio.Application.AbstractCommandline:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Bio.Application.AbstractCommandline:
     |  
     |  parameters = None
    



```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln.txt", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd.txt", "newick")
```


```python
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



Open Computer Vision Basics Notebook
```python
# In this notebook, we used a downloadable package called, Open Computer Vision to perform basic commands concerning phenomics and photo manipulation.
# By using this specific package, we were able to perform comands that changed the colors, shape, and sizes of various images.
```


```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("sunset.jpg")
```


```python
type(img)
```




    numpy.ndarray




```python
img_wrong = cv2.imread('wrong/path/doesnot/abcdefgh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f954d188bd0>




![png](output_7_1.png)



```python
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7f95324ddc50>




![png](output_9_1.png)



```python
img_gray = cv2.imread("sunset.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (190, 265)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7f9540021e10>




![png](output_11_1.png)



```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f953242fc90>




![png](output_12_1.png)



```python
fix_img.shape
```




    (190, 265, 3)




```python
new_img = cv2.resize(fix_img, (1000, 400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f953239fc90>




![png](output_14_1.png)



```python
new_img.shape
```




    (400, 1000, 3)




```python
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7f9530b0dd90>




![png](output_17_1.png)



```python
new_img.shape
```




    (95, 132, 3)




```python
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7f9530af8290>




![png](output_19_1.png)



```python
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7f95309cbed0>




![png](output_20_1.png)



```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite('sunset_fixed_image.jpg', fix_img)
```




    True




```python
img = cv2.imread("sunset.jpg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f95309c1f90>




![png](output_24_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f9530930d50>




![png](output_26_1.png)



```python
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f9530897d90>




![png](output_28_1.png)



```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f9530882cd0>




![png](output_30_1.png)



```python
img1 = cv2.imread("tech.jpg")
img2 = cv2.imread("sunset.jpg")
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f95305c0150>




![png](output_32_1.png)



```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f953051eed0>




![png](output_34_1.png)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f9530488bd0>




![png](output_35_1.png)



```python
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200,1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f9530475a50>




![png](output_39_1.png)



```python
alpha = 0.8
beta = 0.2

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7f9522b410d0>




![png](output_40_1.png)



```python
img1 = cv2.imread('tech.jpg')
img2 = cv2.imread('sunset.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (200,500))
```


```python
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + large_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-116-ef99da0658f4> in <module>
          8 y_end = y_offset + large_img.shape[0]
          9 
    ---> 10 large_img[y_offset:y_end, x_offset:x_end] = small_img
         11 
         12 plt.imshow(large_img)


    ValueError: could not broadcast input array from shape (500,200,3) into shape (190,200,3)



```python

```


```python

```




Corner Detection Notebook
```python
# In this notebook, we used the downloadable package, Open Computer Vision, to detect corners on various images through the performance of specific
# commands regarding the images. We also were able to manipulate the images into a gray scale color scheme while working.
```


```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
flat_chess = cv2.imread('checkered.jpg')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ae591150>




![png](output_2_1.png)



```python
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7ff0ac8c7ed0>




![png](output_3_1.png)



```python
real_chess = cv2.imread("chess.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ac8b52d0>




![png](output_5_1.png)



```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7ff0ac823110>




![png](output_6_1.png)



```python
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ac787e50>




![png](output_8_1.png)



```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ac76bb10>




![png](output_9_1.png)



```python
# Shi-Tomasi Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y), 3, (255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ac6de610>




![png](output_11_1.png)



```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x, y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7ff0ac6b9d90>




![png](output_12_1.png)





Edge Detection Notebook
```python
# In this notebook, we used the downloadable package, Open Computer Vision, to detect edges of objects on various images through the performance of specific
# commands regarding the images. We also were able to manipulate the images into a gray scale color scheme while working, as well as map out edges through
# the view of specific lines and drawings to demonstrate where specific aspects were present in the images chosen.
```


```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("sunset.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f2bd828a090>




![png](output_4_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd51c5210>




![png](output_5_1.png)



```python
med_value = np.median(img)
med_value
```




    115.0




```python
lower = int(max(0, 0.7*med_value))
upper = int(min(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd5138b10>




![png](output_7_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd50ae350>




![png](output_8_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd509c310>




![png](output_9_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd500c090>




![png](output_10_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd4f6de50>




![png](output_11_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd4f5eb90>




![png](output_12_1.png)



```python
blurred_img = cv2.blur(img, ksize = (8,8))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd4da8150>




![png](output_13_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image = blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7f2bd4d84f10>




![png](output_14_1.png)





Feature Matching Notebook
```python
# In this notebook, we used the downloadable package, Open Computer Vision, to detect similar features on various images through the performance of specific
# commands regarding the images. Within this notebook, an image of various features was chosen, as well as an extremely specific image. By 
# running various comands throughout the notebook, we were able to view exactly where the specific image was similar to the image containing
# multiple features crossed in similarity.
```


```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12, 10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
mini_wheats = cv2.imread('miniw.jpg', 0)
display(mini_wheats)
```


![png](output_3_0.png)



```python
cereals = cv2.imread('cereal.jpg', 0)
display(cereals)
```


![png](output_4_0.png)



```python
orb = cv2.ORB_create()

kp1, des1 = orb.detectAndCompute(mini_wheats, mask = None)
kp2, des2 = orb.detectAndCompute(cereals, mask = None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
mini_wheats_matches = cv2.drawMatches(mini_wheats, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(mini_wheats_matches)
```


![png](output_9_0.png)



```python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(mini_wheats, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k = 2)
```


```python
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('Length of total matches:', len(matches))
print('Length of good matches:', len(good))
```

    Length of total matches: 410
    Length of good matches: 25



```python
sift_matches = cv2.drawMatchesKnn(mini_wheats, kp1, cereals, kp2, good, None, flags = 2)
display(sift_matches)
```


![png](output_15_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(mini_wheats, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_params = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(mini_wheats, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```


![png](output_19_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(mini_wheats, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
indx_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_param = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance < 0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0, 255, 0),
                   singlePointColor = (255, 0, 0),
                   matchesMask = matchesMask,
                   flags = 0)
```


```python
flann_matches = cv2.drawMatchesKnn(mini_wheats, kp1, cereals, kp2, matches, None, **draw_params)
                                 
display(flann_matches)
```


![png](output_25_0.png)






Object Detection Notebook
```python
# In this notebook, we used the downloadable package, Open Computer Vision, to detect similar objecta on two images through the performance of specific
# commands regarding the images. Within this notebook, an image of one object was chosen, as well as an image with numerous of that same specific object.
# By running various comands throughout the notebook, we were able to view exactly where the specific object image was similar to the image containing
# numerous objects crossed in similarity.
```


```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
```


```python
%matplotlib inline
```


```python
full = cv2.imread('flowertraining.jpg')
```


```python
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7f89567a93d0>




![png](output_7_1.png)



```python
test = cv2.imread('flowertesting.jpg')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f89566cfa50>




![png](output_10_1.png)



```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (259, 194, 3)
    Training image shape: (183, 276, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, test, method)
    
    mini_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
        
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255, 0, 0), 10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title("Detection of template")
    
    plt.suptitle(m)
    
    plt.show()
    print('/n')
    print('/n')
```


![png](output_13_0.png)


    /n
    /n



![png](output_13_2.png)


    /n
    /n



![png](output_13_4.png)


    /n
    /n



![png](output_13_6.png)


    /n
    /n



![png](output_13_8.png)


    /n
    /n



![png](output_13_10.png)


    /n
    /n



```python
# Heatmap colors did not display as yours did; used to multiple pictures to attempt this!
```


```python

```
