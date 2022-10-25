#!/usr/bin/env python2
import optparse
import re
from collections import defaultdict

################################# Command line options

desc='Correct the GFF3 generated from GUSHR into a BED12 file\n Cite: Cicconardi et al 2023 Tribe-wide genomics reveals the evolutionary dynamics of genome size, content, and selection across the adaptive radiation of Heliconiini butterflies. In prep.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 25-09-2016 - Author: FCicconardi')

parser.add_option('-g', '--GFF3', dest='gff', help='GFF file from GUSHR. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['gff']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

file=open(opts.gff, 'r')

gtf=file.readline()
gtf_trs=defaultdict(list)
gtf_scf=defaultdict(list)

while gtf:
	if gtf.startswith('#'): gtf=file.readline()
	else:
		el=gtf.strip().split('\t')
		scf=el[0]
		start=int(el[3])
		end=int(el[4])
		strand=el[6]
		if start > end:
			start=int(el[4])
			end=int(el[3])
		if el[2] == 'transcript' or el[2] == 'mRNA':
			info=el[-1].split(';')
			#print info
			for i in range(0,len(info)):
				if 'ID=' in info[i]:
					#print info[i]
					if ':' in info[i]:
						transid=info[i].split(':')[1]
					elif '=' in info[i]:
						transid=info[i].split('=')[1]
			score='100'
			key='%s\t%s\t%s\t%s\t%s' % (scf,start,end,transid,score)
			gtf_scf[scf].append(key)
			gtf_trs[key].append(gtf.strip())
			gtf=file.readline()
		elif el[2] == 'gene':
			gtf=file.readline()
		elif el[2] == 'CDS' or el[2].endswith('utr'):
			gtf_trs[key].append(gtf.strip())
			gtf=file.readline()
		else: gtf=file.readline()


for scf in sorted(gtf_scf.keys()):
	for i in range(0,len(gtf_scf[scf])):
		tr=gtf_scf[scf][i]
		ExonsCoords=[]
		AllExStarts=[]
		AllExEnds=[]
		EXONS=defaultdict(list)
		for j in range(0,len(gtf_trs[tr])):
			el=gtf_trs[tr][j].split('\t')
			tmp=tr.split('\t')
			id=tmp[3].replace('_P','-P')
			score=tmp[4]
			tr_start=int(tmp[1])
			scf=el[0]
			start=int(el[3])
			end=int(el[4])
			strand=el[6].strip()
			if el[2] == 'mRNA':
				CDSstart=''
				CDSend=''
				transcript_info='%s\t%s\t%s\t%s\t%s\t%s' % (scf, start-1, end, id, score, strand)
			elif el[2] == 'CDS' or el[2].endswith('utr'):
				EXONS[start].append(gtf_trs[tr][j])
				if CDSstart == '':
					if el[2] == 'CDS': CDSstart=start-1
				elif CDSend == '':
					if strand == '+':
						if el[2] == 'three_prime_utr': CDSend=start-1
					elif strand == '-':
						if el[2] == 'five_prime_utr': CDSend=start-1
		if CDSend == '': CDSend=end
		STARTS=sorted(EXONS.keys())
		if len(EXONS.keys()) == 1:
			ExonsCoords.append((start,end))
			AllExStarts.append(start)
			AllExEnds.append(end)
		elif len(EXONS.keys()) == 2:
			pos=STARTS[0]
			tStart=int(EXONS[STARTS[0]][0].split()[3])
			tEnd=int(EXONS[STARTS[0]][0].split()[4])
			for p in range(1,len(STARTS)):
				pos=STARTS[p]
				pFeat=EXONS[STARTS[p-1]][0].split()[2]
				pStart=int(EXONS[STARTS[p-1]][0].split()[3])
				pEnd=int(EXONS[STARTS[p-1]][0].split()[4])
				tFeat=EXONS[STARTS[p]][0].split()[2]
				tStart=int(EXONS[STARTS[p]][0].split()[3])
				tEnd=int(EXONS[STARTS[p]][0].split()[4])
				if tFeat == pFeat:
					if p == 1:
						ExonsCoords.append((pStart,pEnd))
						AllExStarts.append(pStart)
						AllExEnds.append(pEnd)
					ExonsCoords.append((tStart,tEnd))
					AllExStarts.append(tStart)
					AllExEnds.append(tEnd)
				else:
					if tStart-pEnd == 1:
						ExonsCoords.append((pStart,tEnd))
						AllExStarts.append(pStart)
						AllExEnds.append(tEnd)
        
		elif len(EXONS.keys()) > 2:
			pos=STARTS[0]
			tStart=int(EXONS[STARTS[0]][0].split()[3])
			tEnd=int(EXONS[STARTS[0]][0].split()[4])

			for p in range(1,len(STARTS)-1):
				PrevCoord=''
				pos=STARTS[p]
				pFeat=EXONS[STARTS[p-1]][0].split()[2]
				pStart=int(EXONS[STARTS[p-1]][0].split()[3])
				pEnd=int(EXONS[STARTS[p-1]][0].split()[4])
				tFeat=EXONS[STARTS[p]][0].split()[2]
				tStart=int(EXONS[STARTS[p]][0].split()[3])
				tEnd=int(EXONS[STARTS[p]][0].split()[4])
				nFeat=EXONS[STARTS[p+1]][0].split()[2]
				nStart=int(EXONS[STARTS[p+1]][0].split()[3])
				nEnd=int(EXONS[STARTS[p+1]][0].split()[4])
				if tFeat == pFeat and tFeat == nFeat:
					if p == 1:
						ExonsCoords.append((pStart,pEnd))
						AllExStarts.append(pStart)
						AllExEnds.append(pEnd)
					if (tStart,tEnd) in ExonsCoords: continue
					else:
						ExonsCoords.append((tStart,tEnd))
						AllExStarts.append(tStart)
						AllExEnds.append(tEnd)
				else:
					if tFeat == pFeat and tFeat != nFeat:
						if nStart-tEnd == 1:
							if p == 1:
								ExonsCoords.append((pStart,pEnd))
								ExonsCoords.append((tStart,nEnd))
								AllExStarts.append(pStart)
								AllExStarts.append(tStart)
								AllExEnds.append(nEnd)
							else:
								ExonsCoords.append((tStart,nEnd))
								AllExEnds.append(nEnd)
						elif nStart-tEnd > 1:
							ExonsCoords.append((tStart,tEnd))
							AllExStarts.append(tStart)
							AllExEnds.append(tEnd)

					elif tFeat != pFeat and tFeat != nFeat:
						if tStart-pEnd == 1 and nStart-tEnd == 1:
							if p == 1:
								ExonsCoords.append((pStart,nEnd))
								AllExStarts.append(pStart)
								AllExEnds.append(nEnd)
							elif p > 1:
								PrevCoord=ExonsCoords[-1]
								if pStart == PrevCoord[0]:
									ExonsCoords.remove(PrevCoord)
									ExonsCoords.append((pStart,nEnd))
									AllExStarts.append(pStart)
									AllExEnds.append(nEnd)

					elif tFeat != pFeat and tFeat == nFeat:
						if tStart-pEnd == 1:
							if p == 1:
								ExonsCoords.append((pStart,tEnd))
								AllExStarts.append(pStart)
								AllExEnds.append(tEnd)
							elif p > 1:
								PrevCoord=ExonsCoords[-1]
								if tEnd != PrevCoord[1]:
									ExonsCoords.append((pStart,tEnd))
									AllExStarts.append(pStart)
									AllExEnds.append(tEnd)
						elif tStart-pEnd > 1:
							ExonsCoords.append((tStart,tEnd))
							AllExStarts.append(tStart)
							AllExEnds.append(tEnd)


			if nFeat == tFeat:
				if (nStart,nEnd) not in ExonsCoords:
						ExonsCoords.append((nStart,nEnd))
						AllExStarts.append(nStart)
						AllExEnds.append(nEnd)

			elif nFeat != tFeat and len(EXONS.keys()) == 3:
				PrevCoord=ExonsCoords[-1]
				if nEnd != PrevCoord[1]:
					ExonsCoords.append((tStart,nEnd))
					AllExStarts.append(tStart)
					AllExEnds.append(nEnd)
					
			pos=STARTS[p+1]
			#print p+1,pos,nEnd-nStart+1,EXONS[pos][0]
		lengths=[]
		starts=[]
		TR_start=min(AllExStarts)
		TR_end=max(AllExEnds)
		for Coord in ExonsCoords:
			lengths.append(str(int(Coord[1])-int(Coord[0])+1))
			starts.append(str(int(Coord[0])-TR_start))
		print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t255,0,0\t%s\t%s\t%s' % (scf,TR_start-1,TR_end,id,score,strand,CDSstart,CDSend,len(lengths),','.join(lengths),','.join(starts))

