## Published lncRNA Data directory
These data were gathered from [PLncDB](https://www.tobaccodb.org/plncdb/) and represent high confidence predicted or experimentally validated lncRNAs from various sources (NCBI, RNAcentral, PLncDB, etcetera). Files include both `.fa` sequences and `.gff3` annotations from each individual database. Sequences will be concatonated together into a full database, then condensed using `CD-HIT` to remove similar sequences from multiple sources. 
```bash
cat *.fa >> species_PLncDB_lncRNA.fa
cat *.gff3 >> species_PLncDB_lncRNA.gff3
grep -c ">" species_PLncDB_lncRNA.fa # number of truth sequences in raw dataset
```