#!/usr/bin/env python2

import optparse
from collections import defaultdict

################################# Command line options

desc='Converte a Transdecoder longest_orfs.GFF3 file with its BLAST hit file into a BED12 file of putative multiple CDS regions in a mRNA annotation\n Cite: Cicconardi et al 2023 Tribe-wide genomics reveals the evolutionary dynamics of genome size, content, and selection across the adaptive radiation of Heliconiini butterflies. In prep.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 25-09-2016 - Author: FCicconardi')

parser.add_option('-b', '--Blast-outfmt6', dest='blast', help='BLAST hits (outfmt6) of the longest_orfs file. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-g', '--GFF3-longest_orfs', dest='gff', help='longest_orfs.gff3 file. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['blast','gff']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

blastDCT=defaultdict(list)

tbl=open(opts.blast, 'r').readlines()

for row in tbl:
	qID=row.split()[0]
	rID=row.split()[1]
	score=float(row.split()[11])
	Eval=float(row.split()[10])
	blastDCT[qID].append((rID,score,Eval))


gff=open(opts.gff, 'r').readlines()

for line in gff:
	if len(line) > 1:
		feat=line.split()[2]
		if feat == 'CDS':
			info=line.strip().split()[-1].split(';')
			for el in info:
				if 'Parent' in el:
					cds=el.split('=')[1]
					if cds in blastDCT:
						transc=line.strip().split()[0]
						start=int(line.strip().split()[3])
						end=int(line.strip().split()[4])
						hit=blastDCT[cds][0][0]
						Eval=blastDCT[cds][0][2]
						print '%s\t%s\t%s\t%s\t%s\t+' % (transc,start,end,hit,Eval)
						
		
			
