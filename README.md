# HpGP

## Step 1: Generate a report file (require blastn installed)

./BlastGFF_all.pl 11_165.gff 11_165.list.txt polished_assembly_11_165_2.fasta 26695.fasta > 11_165.report.txt

## Step 2: Generate potential numbers and real numbers for forward and reverse sequences.

Running command: Rscript run.R polished_assembly_NCI_57.fasta NCI_57.list.txt NCI_57.report.txt

