# -*- coding: cp1257 -*-
#Last edit: 19.03.2018
#input:
	#2) a codon alignment
	#3) optional: weights file and codon table
#output
	#codon_usage_bias.tsv
	#observed.tsv
	#proportion_observed.tsv
	#predicted.tsv
	#weighted_predicted.tsv
	#predicted_uniform.tsv
	#weighted_predicted_uniform.tsv
	

'''
codon_usage_bias.tsv
info: 	[number of sequences]	[alignment length]
ttt	PHE	0.4101
ttc	PHE	0.5899
tta	LEU	0.1001

observed.tsv
pos	A	C	G	T	gap	totalACGT
1	49	0	0	0	0	49
2	0	0	0	49	0	49
3	0	0	49	0	0	49
4	0	0	49	0	0	49
5	3	25	21	0	0	49

proportion_observed.tsv
pos	 A 	 C 	 G 	 T	gap
1	1.0	0.0	0.0	0.0	0.0
2	0.0	0.0	0.0	1.0	0.0
3	0.0	0.0	1.0	0.0	0.0

predicted_uniform.tsv
pos	 A 	 C 	 G 	 T
6	0.25	0.25	0.25	0.25
7	0.2667	0.0	0.6	0.1333
8	0.6	0.3333	0.0667	0.0

predicted.tsv
pos	 A 	 C 	 G 	 T
11	0.7143	0.2365	0.0492	0.0
12	0.4243	0.0777	0.4301	0.0679
13	0.0	0.0	1.0	0.0

'''


import sys,getopt
import os
import re
import collections
import copy

#TODO use dictionary dict = {'ttt': ['PHE','F',0]} etc.
#The Standard Code (transl_table=1) 

CODON_TABLE_1={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/3)],'tag':['STOP','*',(1.0/3)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/3)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#The Vertebrate Mitochondrial Code (transl_table=2)
CODON_TABLE_2={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/4)],'tag':['STOP','*',(1.0/4)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',0.5],'tgg':['TRP','W',0.5],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/4)],'cgc':['ARG','R',(1.0/4)],'cga':['ARG','R',(1.0/4)],'cgg':['ARG','R',(1.0/4)],
				 'att':['ILE','I',(1.0/2)],'atc':['ILE','I',(1.0/2)],'ata':['MET','M',(1.0/2)],'atg':['MET','M',(1.0/2)],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['STOP','*',(1.0/4)],'agg':['STOP','*',(1.0/4)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#The Yeast Mitochondrial Code (transl_table=3)
CODON_TABLE_3={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',0.5],'ttg':['LEU','L',0.5],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['THR','T',(1.0/8)],'ctc':['THR','T',(1.0/8)],'cta':['THR','T',(1.0/8)],'ctg':['THR','T',(1.0/8)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/2)],'atc':['ILE','I',(1.0/2)],'ata':['MET','M',(1.0/2)],'atg':['MET','M',(1.0/2)],
				 'act':['THR','T',(1.0/8)],'acc':['THR','T',(1.0/8)],'aca':['THR','T',(1.0/8)],'acg':['THR','T',(1.0/8)],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#The Mold,Protozoan,and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
CODON_TABLE_4={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#5. The Invertebrate Mitochondrial Code (transl_table=5)
CODON_TABLE_5={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',0.25],'cgc':['ARG','R',0.25],'cga':['ARG','R',0.25],'cgg':['ARG','R',0.25],
				 'att':['ILE','I',(1.0/2)],'atc':['ILE','I',(1.0/2)],'ata':['MET','M',(1.0/2)],'atg':['MET','M',(1.0/2)],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/4)],'agc':['SER','S',(1.0/4)],'aga':['SER','S',(1.0/4)],'agg':['SER','S',(1.0/4)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#6. The Ciliate,Dasycladacean and Hexamita Nuclear Code (transl_table=6)
CODON_TABLE_6={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['GLN','Q',(1.0/4)],'tag':['GLN','Q',(1.0/4)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',1.0],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',(1.0/4)],'cag':['GLN','Q',(1.0/4)],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#9. The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
CODON_TABLE_9={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',0.5],'tag':['STOP','*',0.5],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',0.5],'tgg':['TRP','W',0.5],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',0.25],'cgc':['ARG','R',0.25],'cga':['ARG','R',0.25],'cgg':['ARG','R',0.25],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',(1.0/3)],'aac':['ASN','N',(1.0/3)],'aaa':['ASN','N',(1.0/3)],'aag':['LYS','K',1.0],
				 'agt':['SER','S',(1.0/4)],'agc':['SER','S',(1.0/4)],'aga':['SER','S',(1.0/4)],'agg':['SER','S',(1.0/4)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#10. The Euplotid Nuclear Code (transl_table=10)
CODON_TABLE_10={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',(1.0/3)],'tgc':['CYS','C',(1.0/3)],'tga':['CYS','C',(1.0/3)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}
#11. The Bacterial,Archaeal and Plant Plastid Code (transl_table=11)
#No changes compared to standard
CODON_TABLE_11={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/3)],'tag':['STOP','*',(1.0/3)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/3)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}
				 
#12. The Alternative Yeast Nuclear Code (transl_table=12)
CODON_TABLE_12={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/5)],'ttg':['LEU','L',(1.0/5)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/3)],'tag':['STOP','*',(1.0/3)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/3)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/5)],'ctc':['LEU','L',(1.0/5)],'cta':['LEU','L',(1.0/5)],'ctg':['SER','S',(1.0/3)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/3)],'agc':['SER','S',(1.0/3)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#The Ascidian Mitochondrial Code (transl_table=13)
CODON_TABLE_13={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/4)],'cgc':['ARG','R',(1.0/4)],'cga':['ARG','R',(1.0/4)],'cgg':['ARG','R',(1.0/4)],
				 'att':['ILE','I',(1.0/2)],'atc':['ILE','I',(1.0/2)],'ata':['MET','M',0.5],'atg':['MET','M',0.5],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['GLY','G',(1.0/6)],'agg':['GLY','G',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',(1.0/6)],'ggc':['GLY','G',(1.0/6)],'gga':['GLY','G',(1.0/6)],'ggg':['GLY','G',(1.0/6)],
				 '---':['GAP','-',0.0]}

#14. The Alternative Flatworm Mitochondrial Code (transl_table=14)
CODON_TABLE_14={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',(1.0/3)],'tac':['TYR','Y',(1.0/3)],'taa':['TYR','Y',(1.0/3)],'tag':['STOP','*',1.0],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/4)],'cgc':['ARG','R',(1.0/4)],'cga':['ARG','R',(1.0/4)],'cgg':['ARG','R',(1.0/4)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',(1.0/3)],'aac':['ASN','N',(1.0/3)],'aaa':['ASN','N',(1.0/3)],'aag':['LYS','K',1.0],
				 'agt':['SER','S',(1.0/4)],'agc':['SER','S',(1.0/4)],'aga':['SER','S',(1.0/4)],'agg':['SER','S',(1.0/4)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

# Chlorophycean Mitochondrial Code (transl_table=16)
CODON_TABLE_16={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/7)],'ttg':['LEU','L',(1.0/7)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['LEU','L',(1.0/7)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/2)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/7)],'ctc':['LEU','L',(1.0/7)],'cta':['LEU','L',(1.0/7)],'ctg':['LEU','L',(1.0/7)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#Trematode Mitochondrial Code (transl_table=21)
CODON_TABLE_21={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/4)],'cgc':['ARG','R',(1.0/4)],'cga':['ARG','R',(1.0/4)],'cgg':['ARG','R',(1.0/4)],
				 'att':['ILE','I',(1.0/2)],'atc':['ILE','I',(1.0/2)],'ata':['MET','M',(1.0/2)],'atg':['MET','M',(1.0/2)],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',(1.0/3)],'aac':['ASN','N',(1.0/3)],'aaa':['ASN','N',(1.0/3)],'aag':['LYS','K',1.0],
				 'agt':['SER','S',(1.0/4)],'agc':['SER','S',(1.0/4)],'aga':['SER','S',(1.0/4)],'agg':['SER','S',(1.0/4)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}


#The Standard Code (transl_table=1)
CODON_TABLE_22={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/7)],'ttg':['LEU','L',(1.0/7)],
				 'tct':['SER_TCN','S',(1.0/3)],'tcc':['SER_TCN','S',(1.0/3)],'tca':['STOP','*',(1.0/3)],'tcg':['SER_TCN','S',(1.0/3)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/3)],'tag':['LEU','L',(1.0/7)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/3)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/7)],'ctc':['LEU','L',(1.0/7)],'cta':['LEU','L',(1.0/7)],'ctg':['LEU','L',(1.0/7)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#23. Thraustochytrium Mitochondrial Code (transl_table=23) - like trans_table 11
CODON_TABLE_23={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['STOP','*',(1.0/4)],'ttg':['LEU','L',(1.0/5)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/4)],'tag':['STOP','*',(1.0/4)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['STOP','*',(1.0/4)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/5)],'ctc':['LEU','L',(1.0/5)],'cta':['LEU','L',(1.0/5)],'ctg':['LEU','L',(1.0/5)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}

#24. Pterobranchia Mitochondrial Code (transl_table=24)
CODON_TABLE_24={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['TRP','W',(1.0/2)],'tgg':['TRP','W',(1.0/2)],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/4)],'cgc':['ARG','R',(1.0/4)],'cga':['ARG','R',(1.0/4)],'cgg':['ARG','R',(1.0/4)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',(1.0/3)],'aag':['LYS','K',(1.0/3)],
				 'agt':['SER','S',(1.0/3)],'agc':['SER','S',(1.0/3)],'aga':['SER','S',(1.0/3)],'agg':['LYS','K',(1.0/3)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',0.25],'ggc':['GLY','G',0.25],'gga':['GLY','G',0.25],'ggg':['GLY','G',0.25],
				 '---':['GAP','-',0.0]}


#Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
CODON_TABLE_25={'ttt':['PHE','F',0.5],'ttc':['PHE','F',0.5],'tta':['LEU','L',(1.0/6)],'ttg':['LEU','L',(1.0/6)],
				 'tct':['SER_TCN','S',(1.0/4)],'tcc':['SER_TCN','S',(1.0/4)],'tca':['SER_TCN','S',(1.0/4)],'tcg':['SER_TCN','S',(1.0/4)],
				 'tat':['TYR','Y',0.5],'tac':['TYR','Y',0.5],'taa':['STOP','*',(1.0/2)],'tag':['STOP','*',(1.0/2)],
				 'tgt':['CYS','C',0.5],'tgc':['CYS','C',0.5],'tga':['GLY','G',(1.0/5)],'tgg':['TRP','W',1.0],
				 'ctt':['LEU','L',(1.0/6)],'ctc':['LEU','L',(1.0/6)],'cta':['LEU','L',(1.0/6)],'ctg':['LEU','L',(1.0/6)],
				 'cct':['PRO','P',0.25],'ccc':['PRO','P',0.25],'cca':['PRO','P',0.25],'ccg':['PRO','P',0.25],
				 'cat':['HIS','H',0.5],'cac':['HIS','H',0.5],'caa':['GLN','Q',0.5],'cag':['GLN','Q',0.5],
				 'cgt':['ARG','R',(1.0/6)],'cgc':['ARG','R',(1.0/6)],'cga':['ARG','R',(1.0/6)],'cgg':['ARG','R',(1.0/6)],
				 'att':['ILE','I',(1.0/3)],'atc':['ILE','I',(1.0/3)],'ata':['ILE','I',(1.0/3)],'atg':['MET','M',1.0],
				 'act':['THR','T',0.25],'acc':['THR','T',0.25],'aca':['THR','T',0.25],'acg':['THR','T',0.25],
				 'aat':['ASN','N',0.5],'aac':['ASN','N',0.5],'aaa':['LYS','K',0.5],'aag':['LYS','K',0.5],
				 'agt':['SER','S',(1.0/2)],'agc':['SER','S',(1.0/2)],'aga':['ARG','R',(1.0/6)],'agg':['ARG','R',(1.0/6)],
				 'gtt':['VAL','V',0.25],'gtc':['VAL','V',0.25],'gta':['VAL','V',0.25],'gtg':['VAL','V',0.25],
				 'gct':['ALA','A',0.25],'gcc':['ALA','A',0.25],'gca':['ALA','A',0.25],'gcg':['ALA','A',0.25],
				 'gat':['ASP','D',0.5],'gac':['ASP','D',0.5],'gaa':['GLU','E',0.5],'gag':['GLU','E',0.5],
				 'ggt':['GLY','G',(1.0/5)],'ggc':['GLY','G',(1.0/5)],'gga':['GLY','G',(1.0/5)],'ggg':['GLY','G',(1.0/5)],
				 '---':['GAP','-',0.0]}

				 
				 
				 
def readFileToArray(filePath):
	try:
		file = open(filePath, "r")
	except:
		print ("ERROR - Reading " + str(filePath))

	array = []
	while True :
	
		line = file.readline()

		if line == "" :
			break

	   
		array.append(line.split(chr(10))[0])

	file.close()
	return array
	

def fileOpenAndGetSeq(filePath):
	
	seq = ""
	listSeq = []
	
	
	try:
		file = open(filePath, "r")
	except:
		print ("ERROR - Reading " + str(filePath))
		
	entry = []
	while True :
		line = file.readline()

		if line == "" :
			entry = []
			entry.append(name)
			entry.append(seq)
			listSeq.append(entry)
			break
		
		#chr(13)=CR in linux server
		#chr(10)=LF in linux server \n
		line = line.split(chr(10))[0]
		
		if ">" not in line:
			seq = seq+line
			

		if ">" in line:
			if seq != "":
				entry = []
				entry.append(name)
				
				seq = seq.lower();
				seq = seq.replace('u','t')
				
				entry.append(seq)
				listSeq.append(entry)
				#empty seq
				seq = ""
				#new name
				name = line.split(">")[1]
			else:
				name = line.split(">")[1]
	   
		   
	# Faili sulgemine
	file.close()
	return listSeq

def dictClearValues(d):
	if sys.version_info[0] < 3:
		arr = d.iteritems()
	else:
		arr = d.items()
		
	for key, value in arr:
		value[-1] = 0

def fillCodonCounts(codonAlignment,codonCounts):

	
	for a in range(0, len(codonAlignment)):

		sequence = codonAlignment[a][1]
		seqLength = len(sequence) #sequence length of the first sequence or total alignment
		codon = ""
		
		for i in range(0,seqLength):
			nucleotide_removeCR = (sequence[i]).split(chr(13))[0]
			nucleotide = nucleotide_removeCR.split(chr(10))[0]
			nucleotide = nucleotide.lower();
			nucleotide = nucleotide.replace('u','t')
			codon = codon + nucleotide
			if (len(codon) == 3):
				codonCount = codonCounts.get(codon)
				if (codonCount is not None):
					codonCount[-1] += 1 #'ttt':['PHE', 0] => 'ttt': ['PHE', 1]
				codon = ""
				

def calculateCodonUsageBias(codonCounts):

	totalAminoAcidFreq = {}
	
	if sys.version_info[0] < 3:
		arr = codonCounts.iteritems()
	else:
		arr = codonCounts.items()

	for codon, value in arr:
		aminoAcid = value[0]
		count = value[-1]
		result = totalAminoAcidFreq.get(aminoAcid)
		if (result is None):
			totalAminoAcidFreq[aminoAcid] = count
		else:
			totalAminoAcidFreq[aminoAcid] += count
			
	codonUsageBias = {}

	if sys.version_info[0] < 3:
		arr = codonCounts.iteritems()
	else:
		arr = codonCounts.items()
	
	for codon, value in arr:
		val = value[:] # copy value in order to keep old one
		aminoAcid = val[0]
		count = val[-1]
		total = totalAminoAcidFreq.get(aminoAcid)
		prop = 0
		if (total != 0):
			prop = count * 1.0 / total
			val.append (prop)
		codonUsageBias[codon] = val
	
	return codonUsageBias

def codonUsageBiasToSortedList(codonUsageBias):
	CODONLIST = ['ttt','ttc','tta','ttg',
			'tct','tcc','tca','tcg',
			'tat','tac','taa','tag',
			'tgt','tgc','tga','tgg',
			'ctt','ctc','cta','ctg',
			'cct','ccc','cca','ccg',
			'cat','cac','caa','cag',
			'cgt','cgc','cga','cgg',
			'att','atc','ata','atg',
			'act','acc','aca','acg',
			'aat','aac','aaa','aag',
			'agt','agc','aga','agg',
			'gtt','gtc','gta','gtg',
			'gct','gcc','gca','gcg',
			'gat','gac','gaa','gag',
			'ggt','ggc','gga','ggg']
	codonUsageList = []
	for a in range(0,len(CODONLIST)):
		codon = CODONLIST[a]
		result = codonUsageBias.get(codon)
		if (result is not None):
			arr = [codon, result[0], result[2], result[-1]]
			codonUsageList.append(arr)
		else:
			print ('ERROR - could not find codon')
	return codonUsageList
	
	

def createDataFileCodonUsage(filePath, codonUsageBiasList, tag):

	try:
		dataFile = open(filePath + tag +"codon_usage_bias.tsv", "w")
		#File for final user
		dataFileUser = open(filePath + tag +"codon_usage_over_all_sequences.tsv", "w")
		#datafile for js script
		dataFileScript = open(filePath + tag + "codonUsageScript.txt", "w")
	except:
		print ("ERROR - Creating " + str(filePath))
		
	dataFileScript.write("[")
	AminoAcidsInArray = []
	dataFile.write("info: "+ "\t" + str(len(codonAlignment)) +"\t" +  msaLength +"\t" +"\n")
	dataFileUser.write("codon"+ "\t" + "AA" +"\t" +"\t" + "freq" +"\t" + "prop" +"\n")
	for i in range(0, len(codonUsageBiasList)):
		codon = codonUsageBiasList[i][0]
		aminoacid = codonUsageBiasList[i][1]
		count = codonUsageBiasList[i][2]
		prop = round(codonUsageBiasList[i][-1],4)
		
		dataFile.write(codon + "\t")
		dataFile.write(aminoacid + "\t")
		dataFile.write(str(count) + "\t")
		dataFile.write(str(prop) + "\n")

		#File for final user
		if (not aminoacid == "GAP"):
			dataFileUser.write(codon + "\t")
			dataFileUser.write(aminoacid + "\t")
			dataFileUser.write(str(count) + "\t")
			dataFileUser.write(str(prop) + "\n")

		arrayStr = "data: ["
		if ((not aminoacid in AminoAcidsInArray) and aminoacid != "GAP"):
			AminoAcidsInArray.append(aminoacid) #add new aa into array

		#name is codon for highcharts
		dataFileScript.write("\n{name:'"+ codon +"',")
		dataFileScript.write("\n maxPointWidth: 20, \n")
		for j in range(0, 21): #0->20 => len= 21
			if(j < len (AminoAcidsInArray) and aminoacid == AminoAcidsInArray[j] ):
				if (j < 20):
					arrayStr = arrayStr +str(prop*100)+", "
				else:
					arrayStr = arrayStr +str(prop*100)
			else:
				if (j < 20):
					arrayStr = arrayStr + "'" + "" + "',"
				else:
					arrayStr = arrayStr + "'" + ""+ "'"
		dataFileScript.write(arrayStr+"]},")
	dataFileScript.write("]")

	dataFileScript.close()
	dataFile.close()
	dataFileUser.close()
	
def createCodonUsageFrequenceTable(filePath, codonUsageBiasList):

	try:
		dataFileUser = open(filePath + "codon_usage_table.txt", "w")
	except:
		print ("ERROR - Creating " + str(filePath))

	for i in range(0, len(codonUsageBiasList)):

		#File for final user
		if (i<4 or (i>=16 and i<16+4) or (i>=2*16 and i<2*16+4) or (i>=3*16 and i<3*16+4)):
			for j in range(0, 13,4):
				dataFileUser.write(codonUsageBiasList[i+j][0] + " ")
				dataFileUser.write(codonUsageBiasList[i+j][1] + " ")
				dataFileUser.write(str(round(codonUsageBiasList[i+j][-1],2)) + "  ")
			dataFileUser.write("\n")

		if((i+1)%16==0):
			dataFileUser.write("\n")
		


	dataFileUser.close()

def findObsPosition(filePath, fileName, codonAlignment):
	
	allSeq = len(codonAlignment) #Number of sequences in the MSA input
	position = [] #Position there are A,C,G,T and last is gap count
	
	msaLength = len(codonAlignment[0][1])
	for a in range(0,msaLength,1):
		position.append([0,0,0,0,0])

	for i in range (0,allSeq,1):
		sequence = codonAlignment[i][1]
				
		for j in range(0,len(sequence),1):
			#one sequence
			nucleotide_removeCR = (sequence[j]).split(chr(13))[0]
			nucleotide = nucleotide_removeCR.split(chr(10))[0]
			nucleotide = nucleotide.lower()
			nucleotide = nucleotide.replace('u','t') #if RNA sequence
			
			#Compare nucleotides, and add to positions [0,0,0,0,0]
			if nucleotide == "a":
				position[j][0] = position[j][0] + 1
			elif nucleotide == "c" :
				position[j][1] = position[j][1] + 1
			elif nucleotide == "g" :
				position[j][2] = position[j][2] + 1
			elif nucleotide == "t" :
				position[j][3] = position[j][3] + 1
			elif nucleotide == "-" :
				position[j][4] = position[j][4] + 1

	try:
		fileResult = open(filePath + "proportion_" + fileName, "wt")
		fileRawCounts = open(filePath + fileName, "wt")
	except:
		print ("ERROR- opening result file")

	fileResult.write("pos"+"\t"+"A"+"\t"+"C"+"\t"+"G"+"\t"+"T"+"\t"+ "gap"+ "\n")
	fileRawCounts.write("pos"+"\t"+"A"+"\t"+"C"+"\t"+"G"+"\t"+ "T" + "\t"+ "gap" +"\t"+"totalACGT"+ "\n")

	for i in range (0,len(position),1):
		sumACGT = position[i][0]+position[i][1]+position[i][2]+position[i][3]
		gap = position[i][4]

		fileResult.write(str(i+1)+"\t")
		fileRawCounts.write(str(i+1)+"\t")

		controlSum = 0
		for j in range(0,4):
			fileRawCounts.write(str(position[i][j]) + "\t")
			if (sumACGT==0):
				fileResult.write(str("0.0") + "\t")
			else:
				fileResult.write(str(round(position[i][j]*1.0/sumACGT,4)) + "\t")
				controlSum = controlSum + position[i][j]*1.0/sumACGT


		fileResult.write(str(round(position[i][4]*1.0/(gap+sumACGT),4)) + "\n")
		fileRawCounts.write(str(gap) + "\t" + str(sumACGT) + "\n")

		
	fileResult.close()



def createAminoAcidDict(codonUsageBias):

	dict = {} #ARG:{'a':0.0,'c':0.0,'g':0.0,'t':0.0}

	if sys.version_info[0] < 3:
		arr = codonUsageBias.iteritems()
	else:
		arr = codonUsageBias.items()
		
	for codon, value in arr:
		val = value[:] 
		aminoAcid = val[0]
		prop = val[-1]
		aa = dict.get(aminoAcid)
		
		#if (codon[:2] == 'tc' and aminoAcid == 'SER'):
			#aminoAcid == "SER_TCN"
		
		if (aa is None):
			dict[aminoAcid] = { 0:{'a':0.0,'c':0.0,'g':0.0,'t':0.0},
								1:{'a':0.0,'c':0.0,'g':0.0,'t':0.0},
								2:{'a':0.0,'c':0.0,'g':0.0,'t':0.0}}
			dict[aminoAcid][0][codon[0]] = prop
			dict[aminoAcid][1][codon[1]] = prop
			dict[aminoAcid][2][codon[2]] = prop
		else:
			dict[aminoAcid][0][codon[0]] += prop
			dict[aminoAcid][1][codon[1]] += prop
			dict[aminoAcid][2][codon[2]] += prop
			

	return dict
	


def calculatePredictedValues(codonAlignment,codonUsageBias,aaDictionary,weights):

	#print (aaDictionary) #   'MET': {0: {'a': 1.0, 'c': 0.0, 't': 0.0, 'g': 0.0}, 1: {'a': 0.0, 'c': 0.0, 't': 1.0, 'g': 0.0}, 2: {'a': 0.0, 'c': 0.0, 't': 0.0, 'g': 1.0}}
	#print(codonUsageBias)  'ctt': ['LEU', 'L', 0.16666666666666666], 'atg': ['MET', 'M', 1.0], 'aca': ['THR', ....
	#print(codonAlignment) ['HPV65REF', 'atggcagat---------aaag....]
	
	result = {} #  a:0.023,c:0.2...
	msaLength = len(codonAlignment[0][1]) #Length of one sequence
	totalSeq = len(codonAlignment)

	for i in range(0,msaLength,3):
		result[i] = {'a':0.0,'c':0.0,'g':0.0,'t':0.0}
		result[i+1] = {'a':0.0,'c':0.0,'g':0.0,'t':0.0}
		result[i+2] = {'a':0.0,'c':0.0,'g':0.0,'t':0.0}

		for j in range(0,len(codonAlignment)):
			name = codonAlignment[j][0]
			seq = codonAlignment[j][1]
			p1 = seq[i]
			p2 = seq[i + 1]
			p3 = seq[i + 2]
			codon = p1 + p2 + p3
			
			aminoAcid = codonUsageBias.get(codon)
			
			if (codon != '---' and aminoAcid is not None):
				aminoAcid = aminoAcid[0]
				propForCodonEachPos = aaDictionary.get(aminoAcid)
				p_prop = []
				p_prop.append(propForCodonEachPos[0])
				p_prop.append(propForCodonEachPos[1])
				p_prop.append(propForCodonEachPos[2])
				
				weight = 1
				if(weights):
					weight = weights.get(name)
					if(weight is not None):
						weight = float(weight)
				
				for k in range(0,3):
					result[i+k]['a'] += p_prop[k]['a'] * weight
					result[i+k]['c'] += p_prop[k]['c'] * weight
					result[i+k]['g'] += p_prop[k]['g'] * weight
					result[i+k]['t'] += p_prop[k]['t'] * weight

	#normalize
	if sys.version_info[0] < 3:
		arr = result.iteritems()
	else:
		arr = result.items()
		
	for position, value in arr:
		total = sum(value.values(),0.0)
		if sys.version_info[0] < 3:
			values = value.iteritems()
		else:
			values = value.items()
		for nuc, prop in values:
			if (total > 0):
				value[nuc] = round(prop/total,4)
			else:
				value[nuc] = 0
			

	return result
def createDataFile(filePath, predicted):

	try:
		dataFile = open(filePath, "w")
	except:
		print ("ERROR - Creating " + str(filePath))
		
	
	dataFile.write("pos"+"\t"+" A "+"\t"+" C "+"\t"+" G "+"\t"+" T "+ "\n")
	
	
	if sys.version_info[0] < 3:
		arr = predicted.iteritems()
	else:
		arr = predicted.items()

	counter = 1
	for position, value in arr:
		dataFile.write(str(position+1))
		sum = str(value['a']+value['c']+value['g']+value['t'])
		dataFile.write("\t" + str(value['a']) + "\t" + str(value['c']) + "\t" + str(value['g']) + "\t" + str(value['t']))
		if (counter != len(predicted) ):
			dataFile.write("\n")
		counter +=1

	dataFile.close()

def userCodonTableTOcodonCounts(userCodonTable,codonCounts):
			
	codonUsageBias = {}

	if sys.version_info[0] < 3:
		arr = codonCounts.iteritems()
	else:
		arr = codonCounts.items()
		
	for codon, value in arr:
		if (codon == '---'):
			continue
		val = value[:] # copy value in order to keep old one

		prop = 0
		for b in range(0,len(userCodonTable)):
			line = userCodonTable[b]
			line = str(line).lower()
			line = line.replace('u','t') #replace u with t
			line = line.replace('*','STOP') #replace * to stop
			splitedLine = line.split(codon)
			if (len(splitedLine)==2):
				try:
					val[-1] = float(splitedLine[-1].split(" ")[2])
				except:
					break

				codonUsageBias[codon] = val
				break
				
	
	codonSum = {}
	if sys.version_info[0] < 3:
		arr = codonUsageBias.iteritems()
	else:
		arr = codonUsageBias.items()
		
	for codon, value in arr:
		aminoAcid = value[0]
		prop = value[-1]
		res = codonSum.get(aminoAcid)
		if(res is None):
			codonSum[aminoAcid] = prop
		else:
			codonSum[aminoAcid] += prop
	

	if sys.version_info[0] < 3:
		arr = codonUsageBias.iteritems()
	else:
		arr = codonUsageBias.items()
		
	for codon, value in arr:
		aminoAcid = value[0]
		prop = value[-1]
		if(codonSum[aminoAcid] > 0):
			codonUsageBias[codon][-1] = prop/codonSum[aminoAcid]

	

	return codonUsageBias

	
def loadWeightsToDictionary(filePath):
	try:
		file = open(filePath, "r")
	except:
		print ("ERROR - Reading " + str(filePath))

	dict = {}
	while True :
		line = file.readline()
		if line == "" :
			break
		lineList = line.split("\n")
		lineList = lineList[0].split("\t")
		name = lineList[0]
		weight = lineList[1]
		is_in_dict = dict.get(name)
		if(is_in_dict is None):
			dict[name] = weight
	file.close()
	return dict
	

def usage():
	print('python positionReaderTerminal.py -c <codon alignment in fasta format> -w <weights.txt> -g 1 -o /results')
	print('-c, --codonalignmen \t a odon alignment in fasta format')
	print('-w, --weights \t (OPTIONAL) a weights file in tab separated format')
	print('-t, --table \t (OPTIONAL) a codon table with codon usage proportions')
	print('-g, --code \t genetic code identifier (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)')
	print('-o, --output \t dir for output files')
	print('-v, --verbose \t print process steps')
	
def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hp:c:w:t:g:o:v", ["help","codonalignment=","weights=","table=","code=","output=","verpose"])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err)  # will print something like "option -a not recognized"
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt == "-v":
			global verbose
			verbose = True
		elif opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-c", "--codonalignment"):
			global codonAlignmentFilePath
			codonAlignmentFilePath = arg
		elif opt in ("-w", "--weights"):
			global weightsFilePath
			weightsFilePath = arg
		elif opt in ("-t", "--table"):
			global userCodonTablePath
			userCodonTablePath = arg
		elif opt in ("-g", "--code"):
			global codonTable
			codonTable = arg
		elif opt in ("-o", "--output"):
			global outputFolder
			outputFolder = arg
		else:
			assert False, "unhandled option"

if __name__ == '__main__':

	msaLength = "default 0"
	input = None
	output = None
	codonAlignmentFilePath = None
	codonTable = None
	outputFolder = None
	userCodonTablePath = None
	weightsFilePath = None
	weights = False
	verbose = False
	main()


	isUserCodonTable = False
	if(not userCodonTablePath == None):
		isUserCodonTable = True

	if(weightsFilePath):
		try:
			weights = loadWeightsToDictionary(weightsFilePath)
		except:
			weights = False
			print("ERROR - Loading weights")

	seqCount = 0
	codonAlignment = []
	codonAlignment = fileOpenAndGetSeq(codonAlignmentFilePath)
	msaLength = str(len(codonAlignment[0][1]))
	
	if (codonTable == "2"):
		codonCounts = copy.deepcopy(CODON_TABLE_2)
		codonUsageUniform = CODON_TABLE_2
	elif (codonTable == "3"):
		codonCounts = copy.deepcopy(CODON_TABLE_3)
		codonUsageUniform = CODON_TABLE_3
	elif (codonTable == "4"):
		codonCounts = copy.deepcopy(CODON_TABLE_4)
		codonUsageUniform = CODON_TABLE_4
	elif (codonTable == "5"):
		codonCounts = copy.deepcopy(CODON_TABLE_5)
		codonUsageUniform = CODON_TABLE_5
	elif (codonTable == "6"):
		codonCounts = copy.deepcopy(CODON_TABLE_6)
		codonUsageUniform = CODON_TABLE_6
	elif (codonTable == "9"):
		codonCounts = copy.deepcopy(CODON_TABLE_9)
		codonUsageUniform = CODON_TABLE_9
	elif (codonTable == "10"):
		codonCounts = copy.deepcopy(CODON_TABLE_10)
		codonUsageUniform = CODON_TABLE_10
	elif (codonTable == "11"):
		codonCounts = copy.deepcopy(CODON_TABLE_11)
		codonUsageUniform = CODON_TABLE_11
	elif (codonTable == "12"):
		codonCounts = copy.deepcopy(CODON_TABLE_12)
		codonUsageUniform = CODON_TABLE_12
	elif (codonTable == "13"):
		codonCounts = copy.deepcopy(CODON_TABLE_13)
		codonUsageUniform = CODON_TABLE_13
	elif (codonTable == "14"):
		codonCounts = copy.deepcopy(CODON_TABLE_14)
		codonUsageUniform = CODON_TABLE_14
	elif (codonTable == "16"):
		codonCounts = copy.deepcopy(CODON_TABLE_16)
		codonUsageUniform = CODON_TABLE_16
	elif (codonTable == "21"):
		codonCounts = copy.deepcopy(CODON_TABLE_21)
		codonUsageUniform = CODON_TABLE_21
	elif (codonTable == "22"):
		codonCounts = copy.deepcopy(CODON_TABLE_22)
		codonUsageUniform = CODON_TABLE_22
	elif (codonTable == "23"):
		codonCounts = copy.deepcopy(CODON_TABLE_23)
		codonUsageUniform = CODON_TABLE_23
	elif (codonTable == "24"):
		codonCounts = copy.deepcopy(CODON_TABLE_24)
		codonUsageUniform = CODON_TABLE_24
	elif (codonTable == "25"):
		codonCounts = copy.deepcopy(CODON_TABLE_25)
		codonUsageUniform = CODON_TABLE_25
	else:
		codonCounts = copy.deepcopy(CODON_TABLE_1)
		codonUsageUniform = CODON_TABLE_1

		
	
	#Calculate codon usage
	dictClearValues(codonCounts)
	fillCodonCounts(codonAlignment,codonCounts)
	codonUsageProportion = calculateCodonUsageBias(codonCounts)

	#Create codonUsage file and file for highcharts
	codonUsageList = codonUsageBiasToSortedList(codonUsageProportion)
	createDataFileCodonUsage(outputFolder,codonUsageList,"")
	createCodonUsageFrequenceTable(outputFolder,codonUsageList)
	

	if (isUserCodonTable):
		try:
			userCodonTable = readFileToArray(userCodonTablePath)
			codonUsageProportion = userCodonTableTOcodonCounts(userCodonTable,codonCounts) 
			userCodonUsageList = codonUsageBiasToSortedList(codonUsageProportion)
			createDataFileCodonUsage(outputFolder,userCodonUsageList,"user_")
			
		except:
			print("ERROR - Codon Usage Table has wrong format")

	#Observed values
	findObsPosition(outputFolder,"observed.tsv", codonAlignment)
	
	#Predicted values Uniform
	aminoAcidDictUniform = createAminoAcidDict(codonUsageUniform)
	predicted_uniform = calculatePredictedValues(codonAlignment,codonUsageUniform,aminoAcidDictUniform,False)
	createDataFile(outputFolder + "predicted_uniform.tsv", predicted_uniform)
	
	if ( weights ):
		predicted_uniform_weighted = calculatePredictedValues(codonAlignment,codonUsageUniform,aminoAcidDictUniform,weights)
		createDataFile(outputFolder + "weighted_predicted_uniform.tsv", predicted_uniform_weighted)
	
	#Predicted values non-Uniform
	aminoAcidDict = createAminoAcidDict(codonUsageProportion)
	predicted = calculatePredictedValues(codonAlignment,codonUsageProportion,aminoAcidDict,False)
	createDataFile(outputFolder + "predicted.tsv", predicted)
	
	if ( weights ):
		predicted_weighted = calculatePredictedValues(codonAlignment,codonUsageProportion,aminoAcidDict,weights)
		createDataFile(outputFolder + "weighted_predicted.tsv", predicted_weighted)
