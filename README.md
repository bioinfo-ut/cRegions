# cRegions - using scripts in command-line interface

* Make sure you have Python (support for both 2+ and 3+), PERL and R installed

### Convert protein multiple sequence alignment into corresponding codon-based nucleotide alignment using pal2nal.
pal2nal author is Mikita Suyama, more information at http://www.bork.embl.de/pal2nal

Usage: <br>
```
perl  pal2nal.pl  protein.aln  nucleotide.fasta [options]
```

<code>protein.aln<code> - A protein multiple sequence alignment (in FASTA format)<br>
<code>nucleotide.fasta</code> -  Protein-coding sequences (CDS) of the respective proteins (in FASTA format). Also, either mRNA or the full genome can be used instead of CDS. The coding sequence must not contain introns.

Options:<br>
<b>-output (clustal|paml|fasta|codon), default = clustal</b>

<b>-codontable (1(default)|2|3|4|5|6|9|10|11|12|13|14|15|16|21|22|23)</b><br>
    NCBI GenBank codon table<br>
    1  Universal code<br>
    2  Vertebrate mitochondrial code<br>
    3  Yeast mitochondrial code<br>
    4  Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code<br>
    5  Invertebrate mitochondrial<br>
    6  Ciliate, Dasycladacean and Hexamita nuclear code<br>
    9  Echinoderm and Flatworm mitochondrial code<br>
    10  Euplotid nuclear code<br>
    11  Bacterial, archaeal and plant plastid code<br>
    12  Alternative yeast nuclear code<br>
    13  Ascidian mitochondrial code<br>
    14  Alternative flatworm mitochondrial code<br>
    15  Blepharisma nuclear code<br>
    16  Chlorophycean mitochondrial code<br>
    21  Trematode mitochondrial code<br>
    22  Scenedesmus obliquus mitochondrial code<br>
    23  Thraustochytrium mitochondrial code<br>

  (More info see pal2nal.v14/README)<br>
Output:<br>
<code>pal2nal.fasta</code> - A single file with corresponding codon alignment.
<br>



### Henikoff position-based weights.
Henikoff position-based sequence weights are calculated based on the codon alignment.

Usage: <br>
```
python henikoff_weights.py -i pal2nal.fasta -o weights.txt
```

Options:<br>
<b>-i, --input -</b> a codon alignment in FASTA format, (pal2nal output). Can also be RNA sequence<br>
<b>-o, --output -</b> output file, default is weights.txt

Output:<br>
<code>weights.txt</code> - A single file with weights for each sequence.
<br>
Citation: Henikoff S., Henikoff JG. 1994. Position-based sequence weights. Journal of Molecular Biology 243:574–578. DOI: 10.1016/0022-2836(94)90032-9
<br>



### Calculate observed  frequences and predicted proportions for each nucleotide.
text

Usage: <br>
```
python msaPositionReader.py -i pal2nal.fasta -g 1 -w weights.txt -o /results 
```

Options:<br>
<b>'-i, --input -</b> a codon alignment in FASTA format<br>
<b>'-g, --code -</b> genetic code identifier (same as in pal2nal.pl)<br>
<b>'-o, --output -</b> dir for output files<br>
<b>'-w, --weights -</b> (OPTIONAL) a weights file in tab separated format<br>
<b>'-t, --table -</b> (OPTIONAL) a codon table with codon usage proportions. Examples in Codon Usage Database https://www.kazusa.or.jp/codon/ (Codon Usage Table with Amino Acids)<br>



Output:<br>
<code>weighted_predicted.tsv</code> - 
<code>weighted_predicted_uniform.tsv</code> - 
<code>predicted.tsv</code> - 
<code>predicted_uniform.tsv</code> - 
<code>observed.tsv</code> - 
<code>proportion_observed.tsv</code> - 
<code>codonUsageScript.txt</code> - 
<code>codon_usage_table.txt</code> - 
<code>codon_usage_bias.tsv</code> - 
<code>codon_usage_over_all_sequences.tsv</code> - 
<code>codon_usage_over_all_sequences.tsv</code> - 
<br>



## Executing cRegions scripts with an example on non-structural polyprotein of Alphaviruses
```
perl pal2nal.pl  ALPHA_NON_STRUCTURAL_MAFFT.fasta  ALPHA_NON_STRUCTURAL_genome.fasta  -output fasta -codontable 1
python henikoff_weights.py -i pal2nal.fasta -o weights.txt
```


TODO
2) <code>msaPositionReader.py<code><br>


3) <code>msaStatistics.R </code>