# cRegions

## Using cRegions scripts in command-line interface
* Make sure you have Python (support for both 2+ and 3+), PERL and R installed

### Conversion of protein sequence alignments into the corresponding codon-based DNA alignments using pal2nal (author Mikita Suyama, more information at http://www.bork.embl.de/pal2nal)

Usage: 
perl  pal2nal.pl  protein.aln  nucleotide.fasta [options] 

protein.aln - A protein multiple sequence alignment (in FASTA format)
nucleotide.fasta -  Protein-coding sequences (CDS) of the respective proteins (in FASTA format). Also, either mRNA or the full genome can be used instead of CDS. The coding sequence must not contain introns.

Options:
-output (clustal|paml|fasta|codon), default = clustal

-codontable (1(default)|2|3|4|5|6|9|10|11|12|13|14|15|16|21|22|23)<br>
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


(More info see pal2nal.v14/README):


### Henikoff position-based weights
Henikoff position-based sequence weights are calculated based on the codon alignment

Usage: 
python henikoff_weights.py -i pal2nal_output.fasta -o weights.txt

-i, --input - a codon alignment in FASTA format, (pal2nal output). Can also be RNA sequence
-o, --output - output file, default is weights.txt



Citation: Henikoff S., Henikoff JG. 1994. Position-based sequence weights. Journal of Molecular Biology 243:574–578. DOI: 10.1016/0022-2836(94)90032-9


##Executing cRegions scripts with an example on non-structural polyprotein of Alphaviruses
perl pal2nal.pl  ALPHA_FULL_MAFFT.fasta  ALPHA_FULL_genome.fasta  -output fasta -codontable 1
python henikoff_weights.py -i pal2nal.fasta -o weights.txt


TODO
2) msaPositionReader.py
3) msaStatistics.R 