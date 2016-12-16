# HpGP

## Step 1: Generate a report file (require blastn installed)

./BlastGFF_all.pl 11_165.gff 11_165.list.txt polished_assembly_11_165_2.fasta 26695.fasta > 11_165.report.txt

## Step 2: Generate potential numbers and real numbers for forward and reverse sequences.

Rscript run.all.R Rabkin_4023_assembly.fasta 4023.list.txt 4023.report.txt > 4023.summary.txt
