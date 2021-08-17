# HOPE
**HOPE** (**ho**mologous 3’ extension mediated <u>p</u>rime <u>e</u>ditor) utilizes paired prime editing guide RNAs (pegRNAs) encoding the same edits in both sense and antisense DNA strands to achieve high editing efficiency in  human HEK293T cells as well as mismatch repair (MMR)-deficient HCT116 cells

<img src="https://i.loli.net/2021/08/12/oeTBOfb1kCLDlGw.png" alt="image-20210812193249364" style="zoom:40%;" />



## 1.Calculation of editing efficiency

The desired editing efficiency of HOPE (substitution/insertion/deletion) were calculated by CRISPResso2 (https://github.com/pinellolab/CRISPResso2), which is designed to enable rapid and intuitive interpretation of genome editing experiments, comprehensively and powerfully.

Herein, we use CRISPResso2 to calculate the efficiency for substitution and insertion/deletion with “PE (Prime-edited) mode” and “HDR mode”, respectively.

The analysis of substitution events needs following parameters: “—prime editing pegRNA spacer seq”, “—prime editing pegRNA extension seq”, “—prime editing pegRNA scaffold seq”, “—prime editing nicking guide seq”. Specifically, for insertions and deletions, the HDR amplicon sequence was used as HDR reference with “-e” parameter, and “—prime editing pegRNA spacer seq”, “—prime editing nicking guide seq” were also needed. And “—discard indel reads” parameter was used for the calculation of the substitution, insertion, and deletion events, while “—ignore substitutions” was used for calculation of undesired indel type. The precise editing frequencies were calculated based on the “CRISPResso quantification of editing frequency” files that were obtained from CRISPResso2. The corresponding formulas are shown below:
$$
\text { Substitution ratio }=\frac{\text { [\# Unmodified reads (Prime - edited amplicon)] }}{[\text {\# Reads aligned all amplicons }]} \times 100 \%
$$

$$
\text { Insertion or Deletion ratio }=\frac{\text { [\# Unmodified reads (HDR amplicon) }]}{[\text {\# Reads aligned all amplicons }]} \times 100 \%
$$

$$
\text { Indel ratio }=\frac{\text { [\# Modified reads (corresponding targeted amplicon)] }}{[\text {\# Reads aligned corresponding targeted amplicon] }} \times 100 \%
$$



## 2.Indel analysis of HOPE

In order to evaluate the undesired indel after editing by HOPE, the indeL level of each editing result needs to be evaluated. We carry out the indel analysis with “HOPE_indel_analysis.py” (Python3).

```python
python HOPE_indel_analysis.py -ifa Test_ref.fa -iat Test_allele_filt.txt \
					-event_indel Test_event_indel.txt \
					-idx_indel Test_idx_indel.txt \
					-ins Test_ins.txt \
					-del Test_del.txt
```

The script needs two inputs. The "-ifa" parameter is the input of the corresponding amplicon sequence, and the "-iat" is the input filtered from the "Alleles_frequency_table.txt" file analyzed by CRISPResso2. For example, in PE mode, select the table after filting "Primed-edited" and "Modified" as input, for example:

```bash
Aligned_Sequence	Reference_Sequence	Reference_Name	Read_Status	n_deleted	n_inserted	n_mutated	#Reads	%Reads
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGA----------------GGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	16	0	0	69	0.0312864158010003	idx_ 1
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCAC------------------GCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	18	0	0	68	0.0308329894850437	idx_ 2
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATG-GGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	0	1	0	67	0.0303795631690872	idx_ 3
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCAC-GCAAACGGTTGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	1	0	1	67	0.0303795631690872	idx_ 4
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCA-----ACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	5	0	0	65	0.0294727105371742	idx_ 5
```

The results of INDEL can be obtained through analysis, including the positions of insertion, length of insertion sequence, reads statistics and results of insertion sequence.

```bash
Ins_Ref_idx	Length_ins(bp)	c_Reads	Ins_seq
139	1	70	G
33	2	1	TC
181	2	1	AG
154	1	1	A
15	1	1	G
...
```

Furthermore, you can also obtain the start and end positions of the deletion sequence, the deletion length, and the reads statistics.

```bash
Del_Ref_Start	Del_Ref_End	Length_del(bp)	c_Reads
103	120	18	75
103	103	1	73
102	106	5	72
103	113	11	64
117	127	11	64
...
```

The results of cumulative statistical reads of insertion and deletion events at each location of amplicon reference can be obtained.

```bash
Idx	Event_type	c_Reads
...
100	del-count	80
100	ins-count	0
101	del-count	80
101	ins-count	0
102	del-count	152
102	ins-count	0
103	del-count	429
103	ins-count	0
...
```

In addition, you can also obtain statistics on reads that begin to occur at certain locations of the amplicon reference for insertion and deletion.

```bash
idx	ins_count	del_count
...
100	0	70
101	0	1
102	0	72
103	0	277
104	0	0
105	0	51
106	0	0
...
```



## 3. Off-target analysis

The off-target efficiencies for each site could be calculated as descripted above. In order to evaluate the efficiency of off-target, we need to calculate the corresponding background efficiency according to the sites. Here, we choose to remove origin and terminal 15bp of low sequencing quality of the amplicon reference, and remove 25bp before and after around the editing site. Left region used as the background for BG efficiency calculation. Then the editing efficiency of the background is calculated after removing more than 10% of the edits if exsit (which is considered to be SNP) in BG region.

The R script here shows a example for calculating the off-target BG ratio. The input “.bmat” file contain the information of substitution, insertion, deletion and ambigious statistal reads information for each index of amplicon reference, which could be found in Detect-seq pipeline (http://detect-seq.com/).

