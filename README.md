# HOPE
**HOPE** (`ho`mologous 3’ extension mediated `p`rime `e`ditor) utilizes paired prime editing guide RNAs (pegRNAs) encoding the same edits in both sense and antisense DNA strands to achieve high editing efficiency in  human HEK293T cells as well as mismatch repair (MMR)-deficient HCT116 cells

<img src="https://i.loli.net/2021/08/12/oeTBOfb1kCLDlGw.png" alt="image-20210812193249364" style="zoom:40%;" />



## 1.Calculation of editing efficiency

The desired editing efficiencies of `HOPE` (substitution/insertion/deletion) are calculated by CRISPResso2 (https://github.com/pinellolab/CRISPResso2), a software designed to enable rapid and intuitive interpretation of genome editing experiments.

Herein, we use CRISPResso2 to calculate the efficiencies of desired substitution and insertion/deletion with `“PE (Prime-edited) mode”` and `“HDR mode”`, respectively.

The analysis of substitution events needs following parameters: “—prime editing pegRNA spacer seq”, “—prime editing pegRNA extension seq”, “—prime editing pegRNA scaffold seq”, “—prime editing nicking guide seq”. Specifically, for insertions and deletions, the HDR amplicon sequence was used as HDR reference with “-e” parameter, and “—prime editing pegRNA spacer seq”, “—prime editing nicking guide seq” were also needed. And “—discard indel reads” parameter was used for the calculation of the substitution, insertion, and deletion events, while “—ignore substitutions” was used for calculation of undesired indel type. The precise editing frequencies were calculated based on the “CRISPResso quantification of editing frequency” files that were obtained from CRISPResso2. The corresponding formulas are shown below:
$$
\text { Substitution ratio }=\frac{\text { [\# Unmodified reads (Prime - edited amplicon)] }}{[\text {\# Reads aligned all amplicons }]} \times \text{100 \%}
$$

$$
\text { Insertion or Deletion ratio }=\frac{\text { [\# Unmodified reads (HDR amplicon) }]}{[\text {\# Reads aligned all amplicons }]} \times \text{100 \%}
$$

$$
\text { Indel ratio }=\frac{\text { [\# Modified reads (corresponding targeted amplicon)] }}{[\text {\# Reads aligned corresponding targeted amplicon] }} \times \text{100 \%}
$$



## 2.Indel analysis of HOPE

To evaluate the detailed undesired indel types of HOPE, we carry out the indel analysis with “HOPE_indel_analysis.py” (Python3).

For help info, please run `python HOPE_indel_analysis.py -h`:

```bash
$python HOPE_indel_analysis.py -h

usage: HOPE_indel_analysis.py [-h] -ifa INPUT_FASTA -iat INPUT_ALLELE_TABLE -event_indel OUTPUT_EVENT_INDEL -idx_indel OUTPUT_IDX_INDEL -ins
                              OUTPUT_INS_INFO -del OUTPUT_DEL_INFO

Get all the related insertions and deletions information based on the CRISPResso2 allele table result file

optional arguments:
  -h, --help            show this help message and exit
  -ifa INPUT_FASTA, --input_fasta INPUT_FASTA
                        The input reference fasta file
  -iat INPUT_ALLELE_TABLE, --input_allele_table INPUT_ALLELE_TABLE
                        The input CRISPResso2 allele table file
  -event_indel OUTPUT_EVENT_INDEL, --output_event_indel OUTPUT_EVENT_INDEL
                        Output each reference index exist event position plus total indel information
  -idx_indel OUTPUT_IDX_INDEL, --output_idx_indel OUTPUT_IDX_INDEL
                        Output the only the idx indel info, not plus each idx
  -ins OUTPUT_INS_INFO, --output_ins_info OUTPUT_INS_INFO
                        Output each insertion event position, length and count
  -del OUTPUT_DEL_INFO, --output_del_info OUTPUT_DEL_INFO
                        Output each deletion event position, length and count
```

Test example command for calculate indel analysis results.

```bash
python HOPE_indel_analysis.py -ifa Test_ref.fa -iat Test_allele_filt.txt \
	-event_indel Test_event_indel.txt \
  	-idx_indel Test_idx_indel.txt \
  	-ins Test_ins.txt \
  	-del Test_del.txt
```

The script needs two inputs. The ``"-ifa"`` parameter is the input of the corresponding amplicon sequence, and the ``"-iat"`` is the input filtered from the "Alleles_frequency_table.txt" file analyzed by CRISPResso2. For example, filter the column” Reference_Name“ as ”Prime-edited“ and column ”Read_Status“ as ”MODIFIED“.

```bash
Aligned_Sequence	Reference_Sequence	Reference_Name	Read_Status	n_deleted	n_inserted	n_mutated	#Reads	%Reads
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGA----------------GGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	16	0	0	69	0.0312864158010003
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCAC------------------GCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	18	0	0	68	0.0308329894850437
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATG-GGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	0	1	0	67	0.0303795631690872
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCAC-GCAAACGGTTGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	1	0	1	67	0.0303795631690872
CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCA-----ACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	CAGGTGAAGGTGTGGTTCCAGAACCGGAGGACAAAGTACAAACGGCAGAAGCTGGAGGAGGAAGGGCCTGAGTCCGAGCAGAAGAACAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACTAGGGTGGGCAAC	Prime-edited	MODIFIED	5	0	0	65	0.0294727105371742
```

The details of INDEL can be obtained through the pipeline which is mentioned above, including the positions of insertion, length of insertion sequence, reads statistics and insertion sequence.

```bash
Ins_Ref_idx	Length_ins(bp)	c_Reads	Ins_seq
139	1	70	G
33	2	1	TC
181	2	1	AG
154	1	1	A
15	1	1	G
...
```

Furthermore, you can also obtain the start and end positions of the deletion sequence, deletion sequence length, and the reads statistics.

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

Cas-OFFinder (https://github.com/snugel/cas-offinder) is used to predict the off-target sites. The substitution or indel ratio at the predicted off-target editing sites, or a ± 25 bp region surrounding the off-target sites could be calculated as descripted above. Average substitution/indel ratios observed out of potential off-target editing regions (± 25 bp surrounding the predicted off-target sites) are defined as `“BackGround”` (excluding the ± 15 bp from the two side of the amplicon sequence). Note that SNPs are not taken into the calculation (SNP serve as mutation ratio more than 10%).

The R script shows an example for calculating the BG ratio. The input “.bmat” file contain the information of substitution, insertion, deletion and ambiguous statical reads for each location index of amplicon reference, which could be found in Detect-seq pipeline (http://detect-seq.com/).

The `“.bmat”` file looks like:

```bash
chr_name	chr_index	ref_base	A	G	C	T	del_count	insert_count	ambiguous_count	deletion	insertion	ambiguous	mut_num
F-AS-MIS3-10	1	G	0	342389	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	2	G	0	342388	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	3	A	342392	0	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	4	G	0	342392	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	5	A	342392	0	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	6	G	0	342391	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	7	A	342392	0	0	0	0	0	0	.	.	.	0
F-AS-MIS3-10	8	T	0	0	0	342386	0	0	0	.	.	.	0
F-AS-MIS3-10	9	G	0	342392	0	0	1	8	0	AC	A,A,A,A,A,A,A,A	.	0
```

