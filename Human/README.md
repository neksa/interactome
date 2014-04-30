bash-3.2$ grep ">" proteins.fasta | wc -l
   20259
bash-3.2$ grep ">" proteins_with_isoforms.fasta | wc -l
   40144

/usr/local/bin/deltablast -db pdbaa -out HumanDeltaBLAST_hits.xml -outfmt 5 -query /Users/agoncear/data/Uniprot/Proteome/9606/proteins.fasta -inclusion_ethresh 0.05 -domain_inclusion_ethresh 0.05

