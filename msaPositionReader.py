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
	

import sys,getopt
import os
import re
import collections
import copy

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
	name = ""
	listSeq = []
	entry = []
	
	
	try:
		file = open(filePath, "r")
	except:
		print ("ERROR - Reading " + str(filePath))
		
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
				
				seq = seq.lower()
				seq = seq.replace('u','t')
			
				entry.append(seq)
				listSeq.append(entry)
				seq = "" #clear seq
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

def fillCodonCounts(codonAlignment,codonCounts,weights):

	
	for a in range(0, len(codonAlignment)):

		sequence = codonAlignment[a][1]
		name = codonAlignment[a][0]
		seqLength = len(sequence) #sequence length of the first sequence or total alignment
		codon = ""
		weight = 1
		if (weights):
			weightResult = weights.get(name)
			if(weightResult is not None):
					weight = float(weightResult)
		
		for i in range(0,seqLength):
			nucleotide_removeCR = (sequence[i]).split(chr(13))[0]
			nucleotide = nucleotide_removeCR.split(chr(10))[0]
			nucleotide = nucleotide.lower()
			nucleotide = nucleotide.replace('u','t')
			codon = codon + nucleotide
			if (len(codon) == 3):
				codonCount = codonCounts.get(codon)
				if (codonCount is not None):
						codonCount[-1] += weight #'ttt':['PHE', 0] => 'ttt': ['PHE', 1]
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
		dataFileScript = open(filePath + tag + "codonUsageScript.txt", "w")
	except:
		print ("ERROR - Creating " + str(filePath))
		
	dataFileScript.write("[")
	AminoAcidsInArray = []
	dataFile.write("number_of_seq:"+str(len(codonAlignment)) +"\t msa_length:" + msaLength +"\t" +"\n")
	if (tag == "adjusted_"):
		dataFile.write("codon"+ "\t" + "AA" +"\t" +"\t" + "total_weight" +"\t" + "prop" +"\n")
	else:
		dataFile.write("codon"+ "\t" + "AA" +"\t" +"\t" + "freq" +"\t" + "prop" +"\n")
	for i in range(0, len(codonUsageBiasList)):
		codon = codonUsageBiasList[i][0]
		aminoacid = codonUsageBiasList[i][1]
		count = codonUsageBiasList[i][2]
		prop = round(codonUsageBiasList[i][-1],4)
		
		dataFile.write(codon + "\t")
		dataFile.write(aminoacid + "\t")
		dataFile.write(str(count) + "\t")
		dataFile.write(str(prop) + "\n")

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
	
def createCodonUsageFrequenceTable(filePath, codonUsageBiasList, tag):

	try:
		dataFileUser = open(filePath + tag + "codon_usage_table.txt", "w")
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

	dict = {}

	if sys.version_info[0] < 3:
		arr = codonUsageBias.iteritems()
	else:
		arr = codonUsageBias.items()
		
	for codon, value in arr:
		val = value[:] 
		aminoAcid = val[0]
		prop = val[-1]
		aa = dict.get(aminoAcid)
		
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

def normalizeSelectedWeights (dictionary, totalProp):
	#input = {'a': 5.5, 'c': 14.0, 't': 0.1, 'g': 0.125} 
	if sys.version_info[0] < 3:
		arr = dictionary.iteritems()
	else:
		arr = dictionary.items()

	total = 0
	for nuc, prop in arr:
		if ( prop > 1):
			total += prop - 1


	if sys.version_info[0] < 3:
		arr = dictionary.iteritems()
	else:
		arr = dictionary.items()

	for nuc, prop in arr:
		if (total > 0 and prop > 1):
			dictionary[nuc] = round((prop-1)/total,4) * totalProp

	return dictionary

def normalizeWeights (dictionary):
	#input = {'a': 5.5, 'c': 14.0, 't': 2.5, 'g': 0.0} 
	total = sum(dictionary.values(),0.0)

	if sys.version_info[0] < 3:
		arr = dictionary.iteritems()
	else:
		arr = dictionary.items()

	for nuc, prop in arr:
		if (total > 0):
			dictionary[nuc] = round(prop/total,4)
		else:
			dictionary[nuc] = 0
	return dictionary

'''
prediction = {'a': 0.5, 'c': 0, 't': 0.5, 'g': 0.0} Based on all positions
observed = {'a': 0.28, 'c': 0, 't': 0.72, 'g': 0.0} Due to similar sequences
weighted / adjusted predicted should be something = ~ {'a': 0.3, 'c': 0, 't': 0.7, 'g': 0.0}
'''
def calculatePredictedValues(codonAlignment,codonUsageBias,aaDictionary,weights):
	
	predictions = {} #  a:0.023,c:0.2...
	nucWeights = {} #  a:0.023,c:0.2...
	totalSeqInMsa = len(codonAlignment) #Length of one sequence
	defaultWeight = 1.0/totalSeqInMsa
	msaLength = len(codonAlignment[0][1]) #Length of one sequence
	count = -1

	for i in range(0,msaLength,3):
		for k in range(0,3):
			predictions[i+k] = {'a':0.0,'c':0.0,'g':0.0,'t':0.0}
			nucWeights[i+k]   = {'a':{'count':0, 'weight':0.0},'c':{'count':0, 'weight':0.0},'g':{'count':0, 'weight':0.0},'t':{'count':0, 'weight':0.0}} 

		'''
		# ---------- FOR DEBUGGING TODO: Remove-----------
		count +=1
		showPosition = 21
		# ---------- FOR DEBUGGING TODO: Remove-----------
		'''

		for j in range(0,len(codonAlignment)):
			name = codonAlignment[j][0]
			seq = codonAlignment[j][1]
			codon = seq[i] + seq[i + 1] + seq[i + 2]

			codonResult = codonUsageBias.get(codon) #  gta: ['VAL', 'V', 600, 0.35377358490566035] 
			if (codon != '---' and codonResult is not None):
				aminoAcid = codonResult[0]
				propForCodonEachPos = aaDictionary.get(aminoAcid) # VAL': {0: {'a': 0.0, 'c': 0.0, 't': 0.0, 'g': 1.0}, 1: {'a': 0.0, 'c': 0.0, 't': 1.0, 'g': 0.0}, 2: {'a': 0.38502358490566035, 'c': 0.047759433962264154, 't': 0.35377358490566035, 'g': 0.21344339622641509}}
				p_prop = [] # nucleotide probabilites for current amino acid for each position in the codon
				p_prop.append(propForCodonEachPos[0]) # {'a': 0.0, 'c': 0.0, 't': 0.0, 'g': 1.0}
				p_prop.append(propForCodonEachPos[1]) # {'a': 0.0, 'c': 0.0, 't': 1.0, 'g': 0.0}
				p_prop.append(propForCodonEachPos[2]) # {'a': 0.385, 'c': 0.0478, 't': 0.354, 'g': 0.213}

				'''
				# ---------- FOR DEBUGGING TODO: Remove-----------
				if (weights and count == showPosition):
					print("<br>")
					print (name)
					print (aminoAcid)
					print (codon)
					print("Predicted for 3rd pos :  ")
					print (p_prop[-1])
					print("<br>")
				# ---------- FOR DEBUGGING TODO: Remove-----------
				'''

				for k in range(0,3): # 0 1 2
					weight = defaultWeight
					if(weights):
						weightResult = weights.get(name)
						if(weightResult is not None):
							weight = float(weightResult)
						else:
							print('ERROR - No weight found')

					currentNucleotide = codon[k]
					nucWeights[i+k][currentNucleotide]['count'] += 1
					nucWeights[i+k][currentNucleotide]['weight'] += weight
					currentNucleotide = codon[k]
					for x in ['a','c','g','t']:
						predictions[i+k][x] += p_prop[k][x]
		
		#Normalize predictions
		for k in range(0,3):
			predictions[i+k] = normalizeWeights(predictions[i+k])
		
		
		'''
		# ---------- FOR DEBUGGING TODO: Remove-----------
		if (weights and count == showPosition):
			
			print("<br>predictions:<br>")
			print(predictions[i])
			print("<br>")
			print(predictions[i+1])
			print("<br>")
			print(predictions[i+2])
			print("<br>")
			print("<br>nucWeights:<br>")
			print(nucWeights[i])
			print("<br>")
			print(nucWeights[i+1])
			print("<br>")
			print(nucWeights[i+2])
			print("<br>")
		# ---------- FOR DEBUGGING TODO: Remove-----------
		'''


		if(weights):
			for k in range(0,3):
				difNucleotideCounter = 0
				difNucleotideTotal = 0
				totalProp = 0
				for x in ['a','c','g','t']:
					nucCount = nucWeights[i+k][x]['count']
					if (nucCount > 0):
						difNucleotideCounter += 1
						difNucleotideTotal += nucCount
						'''
						#Total proportion that should be divided based on weights. (It is important in cases when a nucleotide is not observed but has prob > 0)
						prediction  {'a': 0.385, 'c': 0.0478, 't': 0.354, 'g': 0.213}
						observed  {'a': 6, 'c': 0, 't': 5, 'g': 3}
						average weights  {'a': 0.3, 'c': 0.1, 't': 0.3, 'g': 0.3}
						totalProp = 0.385 + 0.354 + 0.213 = 0.952
						adjusted values  {'a': 0.385 x 1/0.3 = 1.28 , 't': 0.354 x 1/0.3 = 1.18, 'g': 0.213 x 1/0.3 = 0.71}
						total adjusted = 1.28 + 1.18 + 0.71 = 3.17
						weighted predictions {'a': 1.28 / 3.17 x 0.952 = 0.384 , 'c': 0.0478, 't': 1.18 / 3.17 x 0.952 ...}
						'''

						totalProp += predictions[i+k][x] 

				for x in ['a','c','g','t']:
					nucCount = nucWeights[i+k][x]['count']

					# Weighting enabled only if there is at least two different nucleotides in a column
					if (nucCount > 0 and difNucleotideCounter > 1):
						avgWeight = nucWeights[i+k][x]['weight'] / nucCount
						
						
						'''
						# ---------- FOR DEBUGGING TODO: Remove-----------
						if (weights and count == showPosition):
							print (x + " avgWeight : " + str(avgWeight) + "<br>")
							print (str(predictions[i+k][x] * 1/avgWeight) + "<br>")
						# ---------- FOR DEBUGGING TODO: Remove-----------
						'''


						predictions[i+k][x] = predictions[i+k][x] * 1/avgWeight + 1 # +1 for detecting nucleotides that where weighted in normalizeSelecteWeights

	
				predictions[i+k] = normalizeSelectedWeights(predictions[i+k],totalProp)
		
		
		'''
		# ---------- FOR DEBUGGING TODO: Remove-----------
		if (weights and count == showPosition):

			print("<br>adjusted predictions:<br>")
			print(predictions[i])
			print("<br>")
			print(predictions[i+1])
			print("<br>")
			print(predictions[i+2])
			print("<br>")
		# ---------- FOR DEBUGGING TODO: Remove-----------
		'''


	return predictions

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
	print('python positionReaderTerminal.py -i <codon alignment in fasta format> -w <weights.txt> -g 1 -o /results')
	print('-i, --input \t a codon alignment in FASTA format')
	print('-w, --weights \t (OPTIONAL) a weights file in tab separated format')
	print('-t, --table \t (OPTIONAL) a codon table with codon usage proportions')
	print('-g, --code \t genetic code identifier (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)')
	print('-o, --output \t dir for output files')
	print('-v, --verbose \t print process steps')
	
def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hp:i:w:t:g:o:v", ["help","input=","weights=","table=","code=","output=","verpose"])
	except getopt.GetoptError as err:
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
		elif opt in ("-i", "--input"):
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
	fillCodonCounts(codonAlignment,codonCounts,False)
	codonUsageProportion = calculateCodonUsageBias(codonCounts)
	

	#Create codonUsage file and file for highcharts
	codonUsageList = codonUsageBiasToSortedList(codonUsageProportion)
	createDataFileCodonUsage(outputFolder,codonUsageList,"")
	createCodonUsageFrequenceTable(outputFolder,codonUsageList, "")


	#Calclulate adjusted codon usage
	if ( weights ):
		dictClearValues(codonCounts)
		fillCodonCounts(codonAlignment,codonCounts,weights)
		codonUsageProportionAdjusted = calculateCodonUsageBias(codonCounts)

		#Create codonUsage file and file for highcharts
		codonUsageListAdjusted = codonUsageBiasToSortedList(codonUsageProportionAdjusted)
		createDataFileCodonUsage(outputFolder,codonUsageListAdjusted,"adjusted_")
		createCodonUsageFrequenceTable(outputFolder,codonUsageListAdjusted, "adjusted_")
		
	

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
		predicted_weighted = calculatePredictedValues(codonAlignment,codonUsageProportionAdjusted,aminoAcidDict,weights)
		createDataFile(outputFolder + "weighted_predicted.tsv", predicted_weighted)
