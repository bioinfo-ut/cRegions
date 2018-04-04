# cRegions

##Excecuting cRegions scripts in command-line interface

### 1) Conversion of protein sequence alignments into the corresponding codon-based DNA alignments using pal2nal (author Mikita Suyama, more information http://www.bork.embl.de/pal2nal)

Usage:  pal2nal.pl  protein.aln  nucleotide.fasta [options] 

1) A protein multiple sequence alignment (in FASTA format)
2) Protein-coding sequences (CDS) of the respective proteins (in FASTA format). Also, either mRNA or the full genome can be used instead of CDS. The coding sequence must not contain introns.

Options:

-output (clustal|paml|fasta|codon), default = clustal

-codontable (1(default)|2|3|4|5|6|9|10|11|12|13|14|15|16|21|22|23)
    NCBI GenBank codon table
    1  Universal code
    2  Vertebrate mitochondrial code
    3  Yeast mitochondrial code
    4  Mold, Protozoan, and Coelenterate Mitochondrial code
    and Mycoplasma/Spiroplasma code
    5  Invertebrate mitochondrial
    6  Ciliate, Dasycladacean and Hexamita nuclear code
    9  Echinoderm and Flatworm mitochondrial code
    10  Euplotid nuclear code
    11  Bacterial, archaeal and plant plastid code
    12  Alternative yeast nuclear code
    13  Ascidian mitochondrial code
    14  Alternative flatworm mitochondrial code
    15  Blepharisma nuclear code
    16  Chlorophycean mitochondrial code
    21  Trematode mitochondrial code
    22  Scenedesmus obliquus mitochondrial code
    23  Thraustochytrium mitochondrial code


(More info see pal2nal.v14/README):


Command example:
pal2nal.pl  ALPHA_FULL_MAFFT.fasta  ALPHA_FULL_genome.fasta  -output fasta -codontable 1


1) henikoff_weights.py
2) msaPositionReader.py
3) msaStatistics.R 