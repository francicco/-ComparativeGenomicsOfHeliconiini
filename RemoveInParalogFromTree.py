#!/usr/bin/env python2

import optparse
from Bio import SeqIO
from ete2 import Tree
from collections import defaultdict

################################# Command line options

desc='Remove paralogs from the Ortholog Group (OG)\n Cite: Cicconardi et al 2023 Tribe-wide genomics reveals the evolutionary dynamics of genome size, content, and selection across the adaptive radiation of Heliconiini butterflies. In prep.'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 06-01-2022 - Author: FCicconardi')
parser.add_option('-t', '--in-tree', dest='tree', help='Tree file name.', action='store', type='string', metavar='<FILE>')
parser.add_option('-f', '--fasta', dest='fasta', help='Fasta file name.', action='store', type='string', metavar='<FILE>')

(opts, args) = parser.parse_args()

mandatories = ['tree','fasta']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING!  output file is not specified\n"
                parser.print_help()
                exit(-1)


############################## Reading files and parameters

fasta_dct=defaultdict(list)

for record in SeqIO.parse(opts.fasta, 'fasta'):
	fasta_dct[record.id].append(record.seq)

SpList=defaultdict(list)

tree=open(opts.tree, 'r').readline()

Tree=Tree(tree.strip(), format=1)
c=0
NodeSupport=defaultdict(list)
for node in Tree.traverse("postorder"):
	if node.is_leaf():
		sp=node.name.split('.')[0]
		SpList[sp].append(node.name)
	else:
		if len(node.name) > 1:
			if float(node.name) < 0.9:
				removed_node = node.delete()
			else:
				c+=1
				Nname='Node'+str(c)
				NodeSupport[Nname].append(node.name)
				node.name=Nname

EueidesTaxa=[]
HeliconiusTaxa=[]
print
for sp in SpList:
	MinDist=defaultdict(list)
	AncestorNode=[]
	if len(SpList[sp]) > 1:
		for leaf in SpList[sp]:
			#print leaf
			AncestorNode.append(Tree.search_nodes(name=leaf)[0].up.name)
		if len(set(AncestorNode)) == 1:
			for leaf in SpList[sp]:
				MinDist[Tree.search_nodes(name=leaf)[0].dist].append(leaf)
			bestPara=Tree.search_nodes(name=MinDist[min(MinDist.keys())][0])[0].name
			if bestPara.split('_')[1].startswith('H'):
				HeliconiusTaxa.append(SpList[sp][0])
			elif bestPara.split('_')[1].startswith('E'):
				EueidesTaxa.append(SpList[sp][0])
		else:
			for leaf in SpList[sp]:
				print leaf
				if len(Tree.search_nodes(name=leaf)[0].up.name) == 0: print 'Root'
				else:
					print Tree.search_nodes(name=leaf)[0].up.name, NodeSupport[Tree.search_nodes(name=leaf)[0].up.name][0]
				print Tree.search_nodes(name=leaf)[0].up
				print
				if leaf.split('_')[1].startswith('H'):
					HeliconiusTaxa.append(SpList[sp][0])
				elif leaf.split('_')[1].startswith('E'):
					EueidesTaxa.append(SpList[sp][0])
	else:
		if SpList[sp][0].split('_')[1].startswith('H'):
			HeliconiusTaxa.append(SpList[sp][0])
		elif SpList[sp][0].split('_')[1].startswith('E'):
			EueidesTaxa.append(SpList[sp][0])

print

SpList=[]
for sp in EueidesTaxa:
	SpList.append(sp.split('.')[0])

if len(SpList) == len(set(SpList)):
	out=open(opts.fasta.replace('NT.fasta','Eueides.NT.fasta'), 'w')
	for sp in EueidesTaxa:
		if len(str(fasta_dct[sp][0])) > 1:
			print >> out, '>%s\n%s' % (sp,str(fasta_dct[sp][0]).replace('-',''))

SpList=[]
for sp in HeliconiusTaxa:
	SpList.append(sp.split('.')[0])

if len(SpList) == len(set(SpList)):
	out=open(opts.fasta.replace('NT.fasta','Heliconius.NT.fasta'), 'w')
	for sp in HeliconiusTaxa:
		if len(str(fasta_dct[sp][0])) > 1:
			print >> out, '>%s\n%s' % (sp,str(fasta_dct[sp][0]).replace('-',''))


for node in Tree.traverse("postorder"):
	if node.name in NodeSupport:
		node.name=NodeSupport[node.name][0]

treeout=open(opts.tree.replace('NT.fasta.treefile','NT.Pruned.treefile'), 'w')
print >> treeout, Tree.write(format=1)
