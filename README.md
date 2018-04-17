# cRegions - using scripts in command-line interface

* Make sure you have Python (support for both 2+ and 3+), PERL and R installed

### Convert protein multiple sequence alignment into corresponding codon-based nucleotide alignment using pal2nal.
pal2nal author is Mikita Suyama, more information at http://www.bork.embl.de/pal2nal

Usage: <br>
```
perl  pal2nal.pl  protein.aln  nucleotide.fasta [options]
```

<code>protein.aln</code> - A protein multiple sequence alignment (in FASTA format)<br>
<code>nucleotide.fasta</code> -  Protein-coding sequences (CDS) of the respective proteins (in FASTA format). Also, either mRNA or the full genome can be used instead of CDS. The coding sequence must not contain introns. <br>

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
<br>
Citation: Henikoff S., Henikoff JG. 1994. Position-based sequence weights. Journal of Molecular Biology 243:574–578. DOI: 10.1016/0022-2836(94)90032-9
<br>
<br>



### Calculate observed and predicted nucleotide frequences or proportions for each position in the codon alignment.

Next step is to calculate observed frequencies for each nucleotide (A;C;G;T) in each position in the codon alignment. Also, predicted nucleotide values are calculated based on the codon usage bias of genes analysed or assuming uniform codon usage of synonymus codons.

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
<code>observed.tsv</code> - Observed nucleotide frequencies for each position in the codon alignment.<br>
<code>proportion_observed.tsv</code> - Observed nucleotide proportions for each position in the codon alignment.<br>
<code>predicted.tsv</code> - Predicted nucleotide proportions for each position in the codon alignment using codon usage bias of these genes.<br>
<code>weighted_predicted.tsv</code> - Predicted nucleotide proportions for each position in the codon alignment using codon usage bias of these genes. Results are weighted using values in the weights.txt file.<br>
<code>predicted_uniform.tsv</code> - Predicted nucleotide proportions for each position in the codon alignment assuming uniform codon usage.<br>
<code>weighted_predicted_uniform.tsv</code> - Predicted nucleotide proportions for each position in the codon alignment assuming uniform codon usage. Results are weighted using values in the weights.txt file.<br>
<code>codon_usage_table.txt</code> - Codon usage bias in a table format, can be used as an input in msaPositionReader.py script ot the cRegions webpage.<br>
<code>codon_usage_bias.tsv</code> -  Codon usage bias tsv file. Information about the number of sequences and the length of the MSA is provided from which the codon usage was calculated.<br>
<code>codonUsageScript.txt</code> - A file used as an input for highcharts plots.<br>
<br>
<br>


### Calculate three different metrics based on observed and predicted values

cRegions algorithm uses three different metrics to compare observed and predicted values.<br>
1)  Chi-Square Goodness of Fit Test (chisq.test) for each column in the codon alignment. This test tries to fit a statistical model (predicted nucleotide proportions) to the observed data, estimating how “close” observed values are to expected values. Bonferroni correction is used to show the threshold with significance level p = 0.05. We used Bonferroni correction as it is most strict to avoid false positive hits.
2) Root-mean-square deviation (RMSD). <br>
3) Maximum difference (MAXDIF), which selects only a single nucleotide for each column that has the highest absolute difference between predicted and observed values. <br><br>


Usage: <br>
```
Rscript msaStatistics.R [mode] [allowed gap] [skip gap] [window size] [codon position] [input/output dir] observed.tsv proportion_observed.tsv weighted_predicted_uniform.tsv weighted_predicted.tsv codon_usage_bias.tsv
```
Options:<br>
<b>mode - </b> 1 - Single positions; 2 - sliding window mode; 3 both (web tool uses mode 3)<br>
<b>allowed gap - </b> how many gaps (percentage) are allowed in a single column in the codon alignment, default 20 <br>
<b>skip gap - </b> when a position in the sliding window mode in the codon aligment has more than 'skip gap' (percentage), then this position is excluded from the window and the next one is included , default 90 ( applicable when window size > 1 ) <br>
<b>window size - </b> windows size in the sliding window mode, default 1<br>
<b>codon position - </b> position in three-nucleotide codon which is used in the arithmetic mean calculation, default 3<br>
<b>input/output dir - </b> input files (observed.tsv, predicted.tsv etc. ) must be located at this directory and also results are written into this directory <br>
<code>weighted_predicted_uniform.tsv</code> - predicted_uniform.tsv can be also used instead of weighted_predicted_uniform.tsv<br>
<code>weighted_predicted.tsv</code> - predicted.tsv can be also used instead of weighted_predicted.tsv <br>


Output:<br>
<code>raw_values.tsv</code> - Metric values (RMSD and MAXDIF) and p-value for Chi-Square Goodness of Fit Test for each position in the codon alignment (also for uniform codon usage). NA means that metric could not be calculated (e.g., more gaps than allowed, less than two nucleotides in one position). <br>
<code>raw_values_sliding_window.tsv</code> - Metric values (RMSD and MAXDIF) and p-value for Chi-Square Goodness of Fit Test for each position in the codon alignment in the sliding window mode. <br>

* msaStatistics.R generates multiple other files for webtool.

<br>
<br>

### Executing cRegions scripts on the non-structural polyprotein of Alphaviruses as an example

```
perl pal2nal.pl  ALPHA_NON_STRUCTURAL_MAFFT.fasta  ALPHA_NON_STRUCTURAL_genome.fasta  -output fasta -codontable 1
python henikoff_weights.py -i pal2nal.fasta -o weights.txt
python msaPositionReader.py -i pal2nal.fasta -g 1 -w weights.txt -o /results
Rscript msaStatistics.R 3 20 90 1 3 /results observed.tsv proportion_observed.tsv weighted_predicted_uniform.tsv weighted_predicted.tsv codon_usage_bias.tsv
```