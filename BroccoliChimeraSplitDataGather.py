#!/usr/bin/env python2

import optparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

################################# Command line options

desc='Filter out fragments from OGs\n Cite: Cicconardi et al 2023 Tribe-wide genomics reveals the evolutionary dynamics of genome size, content, and selection across the adaptive radiation of Heliconiini butterflies. In prep.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 04-05-2021 - Author: FCicconardi')
parser.add_option('-c', '--Broccoli-ChimeraTbl', dest='tbl', help='chimeric_proteins.txt from Broccoli', action='store', metavar='<FILE>')
parser.add_option('-g', '--Broccoli-table_OGs', dest='ogs', help='table_OGs_protein_names.txt from Broccoli', action='store', metavar='<FILE>')
parser.add_option('-f', '--Fasta', dest='fasta', help='Proteomic fasta file for all the species', action='store', metavar='<ARG>')

(opts, args) = parser.parse_args()

mandatories = ['tbl','ogs','fasta']
for m in mandatories:
	if not opts.__dict__[m]:
		print "\nWARNING!  output file is not specified\n"
		parser.print_help()
		exit(-1)

############################## Reading files and parameters

AssInfo={'Dple':'FragChromosomes','Bany':'FragChromosomes','Jcoe':'FragChromosomes','Mcin':'Chromosomes','Smor':'Chromosomes','Ptel':'FragChromosomes','Dpha':'FragChromosomes','Diul':'Chromosomes','Pdid':'FragChromosomes','Avfl':'FragChromosomes','Avcr':'Scaffolds','Avpe':'FragmentedScaffolds','Djun':'FragChromosomes','Eali':'FragmentedScaffolds','Etal':'Scaffolds','Elyb':'FragmentedScaffolds','Eisa':'Chromosomes','Elam':'FragmentedScaffolds','Evib':'FragmentedScaffolds','Hhor':'FragmentedScaffolds','Hcly':'FragmentedScaffolds','Htel':'Scaffolds','Hhec':'Scaffolds','Hher':'Scaffolds','Herd':'Chromosomes','Hpet':'FragmentedScaffolds','Heet':'FragmentedScaffolds','Hlat':'Scaffolds','Hhim':'Chromosomes','Hric':'FragmentedScaffolds','Hper':'FragmentedScaffolds','Hcha':'FragChromosomes','Hert':'FragmentedScaffolds','Hdem':'Scaffolds','Hleu':'FragmentedScaffolds','Hsar':'Chromosomes','Hant':'FragmentedScaffolds','Hele':'FragmentedScaffolds','Hcon':'FragmentedScaffolds','Hsap':'FragmentedScaffolds','Hhew':'FragmentedScaffolds','Haoe':'FragChromosomes','Hheb':'FragmentedScaffolds','Hhie':'FragmentedScaffolds','Hdor':'FragChromosomes','Hxan':'FragmentedScaffolds','Hege':'FragmentedScaffolds','Hbur':'Scaffolds','Hwal':'FragmentedScaffolds','Hnat':'FragChromosomes','Hmel':'Chromosomes','Hcyd':'Chromosomes','Hpac':'FragmentedScaffolds','Htim':'FragmentedScaffolds','Hheu':'FragmentedScaffolds','Hbes':'Scaffolds','Hism':'FragmentedScaffolds','Hnum':'Scaffolds','Heth':'FragmentedScaffolds','Hhel':'Scaffolds','Hatt':'FragmentedScaffolds','Helv':'Scaffolds','Hpar':'Scaffolds'}



RefSP=['Eisa','Jcoe','Dpha','Diul','Djun','Hmel','Hsar','Pdid','Hcyd','Hcha']


print 'Read Proteome'
ProtDCT=defaultdict(list)
fasta_file=open(opts.fasta, 'r')
for record in SeqIO.parse(fasta_file, 'fasta'):
	ProtDCT[record.id].append(str(record.seq))


OGdict=defaultdict(list)

indcol=open(opts.ogs, 'r').readline().split()[0]

OGS = pd.read_csv(opts.ogs, sep='\t', index_col=indcol, low_memory=False)

for OG, row in OGS.iterrows():
	j=-1
	for i in row:
		j+=1
		if 'fasta' in row.index[j]:
			sp=row.index[j].replace('.fasta','')
		else: sp=row.index[j]
		if sp in RefSP:
			OGdict[OG].append(i) 

ChimTbl=open(opts.tbl, 'r').readlines()

OutManCheck=open('PutativeChimeric.ToManualCheck.AA.fasta', 'w')


print 'Read', opts.tbl
for el in ChimTbl:
	if el.startswith('#'): continue
	else:
		if int(el.strip().split('\t')[2]) >= 2:
		#if len(el.strip().split('\t')) > 3:
			print el.strip()
			OGs=el.strip().split('\t')[3]
			fID=el.split()[1]
			sp=fID.split('.')[0]
			if AssInfo[sp] != 'Chromosomes' or AssInfo[sp] != 'FragChromosomes':

				Chimfasta_out=open('Chimeric.'+fID+'.AA.fasta', 'w')
				ChFasta=ProtDCT[fID][0]
				print >> Chimfasta_out, '>%s\n%s' % (fID,ChFasta)
				
				OGfasta_out=open('OGs.'+fID+'.AA.fasta', 'w')
			
				for og in OGs.split():
					if og in OGdict:
						GOT=0
						for rf in RefSP:
							if sp != rf:
								if GOT == 0:
									for i in range(0,len(OGdict[og])):
										RefID=OGdict[og][i]
										if RefID.startswith(rf):
											if RefID in ProtDCT:
												SeqAA=ProtDCT[RefID][0]
												print '\t',og,RefID
												print >> OGfasta_out, '>%s_%s\n%s' % (og,RefID,SeqAA)
												GOT=1
												break
		elif int(el.strip().split('\t')[2]) < 2:
			fID=el.split()[1]
			print >> OutManCheck, '>%s\n%s' % (fID,ProtDCT[fID][0])
			






