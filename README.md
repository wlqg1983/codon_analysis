#### script: PR2.py ####

'''
Codon Third Base Content Analyzer

Description: Calculate the proportion of third base (A/U/G/C) in specified codons,
             used for codon usage bias analysis.

Input format:
  Gene_Name1
  AAA 10 AAC 5 AAG 8 AAU 3
  ACA 7 ACC 12 ACG 4 ACU 6
  (codon frequency pairs, space-separated)
  
  Gene_Name2
  GCA 15 GCC 8 GCG 5 GCU 10
  ...

Output:
  {input_name}_with_labels.txt  : Intermediate file with gene labels
  {input_name}_out.txt          : Result table (Gene A U G C)

Target codons (32 codons):
  GCA GCC GCG GCU  GGA GGC GGG GGU
  CUA CUC CUG CUU  CCA CCC CCG CCU
  CGA CGC CGG CGU  UCA UCC UCG UCU
  ACA ACC ACG ACU  GUA GUC GUG GUU

Calculation:
  For each gene, extract frequencies of target codons,
  calculate proportion of A/U/G/C at the 3rd position.

Usage:
  python script.py <input_files> [input_files...]

Input methods:
  1. Single file    : python script.py data.txt
  2. Wildcards      : python script.py *.txt
  3. Multiple files : python script.py file1.txt file2.txt

Output files location:
  Same directory as input file

Dependency:
  No external dependencies (uses only Python standard library)

Example output (out.txt):
  Gene        A       U       G       C
  rps12       0.3245  0.2876  0.1987  0.1892
  rpl16       0.2987  0.3124  0.1876  0.2013
'''


########################################################################################

#### script: CDS-final.py ####

'''
FASTA Gene Sequence Filtering and Deduplication Tool

Description: Filter complete CDS sequences (length>300bp + ATG start + stop codon), 
             keep the longest sequence per gene

Usage:
  python script.py <input> [input...]

Input methods:
  1. Single file    : python script.py input.fasta
  2. Wildcards      : python script.py *.fasta *.fas
  3. Directory      : python script.py /path/to/dir/

Output:
  {original_filename}_filtered_{unique_gene_count}.fasta

Filtering criteria:
  - Sequence length > 300 bp
  - Starts with ATG (start codon)
  - Ends with TAA/TAG/TGA (stop codon)

Deduplication rules:
  - Use the first word of FASTA header as gene name
  - Keep the longest sequence for each gene

Dependency:
  pip install biopython
'''

########################################################################################
#### citation ####

XU Rui, GUAN Xinyue, LIU Jing, YIN Xiujuan, WANG Liqiang. (2026) Comprehensive Analysis of Chloroplast Genomes from 15 Mangrove Species (Rhizophoraceae): Codon Usage, Phylogenetic Relationships, and Evolutionary Insights. 
