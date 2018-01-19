#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
try:
	import pip
except ImportError, e:
	print "Module pip not found!"
	print "Please install pip manually to proceed!"
	sys.exit(1)
def install(package):
    pip.main(['install', package])

for mod in ['pip','scipy','string','math','re','csv',
            'sys','os','commands','datetime','operator',
            'getopt','pickle','shutil','glob','types',
            'math','copy','pyExcelerator','xlrd','xlwt','xlutils','types']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)
		#pass # module doesn't exist, deal with it.


import numpy as np
import unicodedata
import textwrap as tw

import itertools as IT


from Bio import pairwise2
from Bio import SeqIO

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna, generic_protein

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML




class NestedDict(dict):
	def __getitem__(self, key):
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

def numerize(s):
	try:
		if s=='NAN':
			return s
		float(s)
		if float(s).is_integer():
			return int(float(s))
		elif float(s)==0:
			return float(s)
		else:
			return float(s)

	except ValueError:
		return s


def filename(ifile):
	if ifile.split('.')[0] == '':
		ipat = ''
		iname = ''
		itype = ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep = "\\"
		elif "/" in ifile.split('.')[0]:
			sep = "/"
		else:
			ipat = ''
			iname = ifile.split('.')[0]
			itype = ifile.split('.')[1]
			return ipat, iname, itype
		allpath = ifile.split('.')[0]
		iname = allpath.split(sep)[-1]
		ipath = allpath.split(sep)[:-1]
		ipat = '/'.join(ipath)
		itype = ifile.split('.')[1]
	return ipat, iname, itype


def readcsv(ifile,delim=','):
	f=open(ifile,'r')
	sheet=[]
	rdr=csv.reader(f, delimiter=delim)
	data=[ln for ln in rdr]
	f.close()
	return data

def maptable(ifile):
	print ifile
	info=NestedDict()
	ilist=[]
	ipt, inm, itp = filename(ifile)

	data=readcsv(ifile)
	headers=data[0]
	#Automatically find variables
	headin={ hd : headers.index(hd) for hd in headers}
	nec=['File','Plate','Type']

	#filein=headers.index('File')
	# platein=headers.index('Plate')
	# strainin=headers.index('Strain')
	# typein=headers.index('Type')
	#print metin, ecoin,plate,well
	for ln in data[1:]:
		index=ln[headin['Plate']]+'-'+ln[headin['Well']]
		for hd in headin.keys():
			info[index][hd]=numerize(ln[headin[hd]].strip().encode('ascii','ignore'))
		#print info[fl]['Replicate']

	#print genes
	return info




def index_genbank_features(gb_record, feature_type, qualifier) :
	answer = dict()
	for (index, feature) in enumerate(gb_record.features) :
		if feature.type==feature_type :
			if qualifier in feature.qualifiers :
				#There should only be one locus_tag per feature, but there
				#are usually several db_xref entries
				for value in feature.qualifiers[qualifier] :
					if value in answer :
						print "WARNING - Duplicate key %s for %s features %i and %i" % (value, feature_type, answer[value], index)
					else:
						answer[value] = index
	return answer




def findamplicons(hspsmap,strand,ampmin,ampmax):
	ordkeys=sorted(hspsmap.keys())
	amplicons=NestedDict()
	if len(hspsmap.keys())>1:
		for ind1, pos1 in enumerate(ordkeys):
			qs1=hspsmap[pos1].query_start
			qe1=hspsmap[pos1].query_end
			ts1=hspsmap[pos1].sbjct_start
			te1=hspsmap[pos1].sbjct_end
			for ind2, pos2 in enumerate(ordkeys):
				qs2 = hspsmap[pos2].query_start
				qe2 = hspsmap[pos2].query_end
				ts2 = hspsmap[pos2].sbjct_start
				te2 = hspsmap[pos2].sbjct_end
				if ind1!=ind2 and qs2>qe1:
					index='{}:{}'.format(ind1,ind2)
					if strand=='-1':
						dist=ts1-te2
					else:
						dist=te2-ts1
					if dist>ampmin and dist<ampmax:
						amplicons[dist]['Index']=index
						amplicons[dist]['Primer1']=hspsmap[pos1]
						amplicons[dist]['Primer2']=hspsmap[pos2]
	else:
		pos=ordkeys[0]
		qs = hspsmap[pos].query_start
		qe = hspsmap[pos].query_end
		ts = hspsmap[pos].sbjct_start
		te = hspsmap[pos].sbjct_end
		if qe>30:
			index='X:0'
			amplicons[0]['Primer2'] = hspsmap[pos]
		else:
			index = '0:X'
			amplicons[0]['Primer1'] = hspsmap[pos]
		amplicons[0]['Index'] = index

	return amplicons




def printamplicons(amplicons,strand):
	#header=[]
	table=[]
	for ampk in sorted(amplicons.keys()):
		#print amplicons[ampk]['Index']
		#prim1_id,prim2_id=amplicons[min(amplicons.keys())]['Index'].split(':')
		pr1_ind,pr2_ind=amplicons[ampk]['Index'].split(':')
		if pr1_ind!='X':
			pr1=amplicons[ampk]['Primer1']
			qs1 = pr1.query_start
			qe1 = pr1.query_end
			ts1 = pr1.sbjct_start
			te1 = pr1.sbjct_end
			pr1q=pr1.query
			pr1s=pr1.sbjct
			pr1m=pr1.match
			e1=pr1.expect
		else:
			qs1 = 0
			qe1 = 0
			ts1 = 0
			te1 = 0
			pr1q=''
			pr1s=''
			pr1m=''
			e1=''
			
		if pr2_ind!='X':
			pr2=amplicons[ampk]['Primer2']
			qs2 = pr2.query_start
			qe2 = pr2.query_end
			ts2 = pr2.sbjct_start
			te2 = pr2.sbjct_end
			pr2q=pr2.query
			pr2s=pr2.sbjct
			pr2m=pr2.match
			e2=pr2.expect

		else:
			qs2 = 0
			qe2 = 0
			ts2 = 0
			te2 = 0
			pr2q=''
			pr2s=''
			pr2m=''
			e2=''


		if strand == '-1':
			query_header = "{}:{}({})".format(qe2, qs2,e2) + ' ' * ((qe2-qs2) +15) + "{}:{}({})".format(
				qe1, qs1,e1)
			subject_footer = "{}:{}".format(te2, ts2) + ' ' * 20 + "{}:{}".format(
				te1, ts1)
			query = pr2q + '.' * 20 + pr1q
			sbjct = pr2s + '.' * 20 + pr1s
			match = pr2m + '<' * 20 + pr1m
		else:
			query_header = "{}:{}(e={})".format(qs1, qe1,e1) + ' ' * ((qe1-qs1)+15) + "{}:{}(e={})".format(
				qs2, qe2,e2)
			subject_footer = "{}:{}".format(ts1, te1) + ' ' * 20 + "{}:{}".format(
				ts2, te2)
			query = pr1q + '.' * 20 + pr2q
			sbjct = pr1q + '.' * 20 + pr2s
			match = pr1m + '>' * 20 + pr2m


		print '\n\tAmplicon length: {}'.format(ampk)
		print('\t'+query_header)
		print('\t'+query)
		print('\t'+match)
		print('\t'+sbjct)
		print('\t'+subject_footer)
		print '\tMapping: {}..{}\n'.format(ts1, te2)
		ln=[0 if ampk==0 else len(amplicons.keys()),ampk,'{}..{}'.format(ts1, te2), \
		    e1,qs1, qe1,ts1, te1, e2,qs2, qe2, ts2, te2]
		table.append(ln)
	return table



def csvwriter(table,oname,sep=','):
	f=open(oname,"wb")
	ofile=csv.writer(f, delimiter=sep) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
	for row in table:
		#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
		ofile.writerow(row)
	f.close()


def getindices(s):
    return [i for i, c in enumerate(s) if c.isupper()]


def orderfeatures(fragment,strand,alen):
	ftheader=['Type','Distance','Start','End','Compound','Group']
	feats = [feat for feat in fragment.features if feat.strand == strand and feat.type in ['Operon', 'gene']]
	fmap={}
	features={}
	isgene=False
	isoperon=False
	overlap=False
	flen = fragment.__len__()

	#Order by distance
	fmap={fid:ftpos(ft,flen)[0]-alen for fid,ft in enumerate(feats)}
	order = sorted(fmap, key=fmap.get)
	ofeats = [feats[i] for i in order]


	fmap2={fid:ftpos(ft,flen)[0]-alen for fid,ft in enumerate(ofeats)}

	grouping=groupfeatures(fmap2,2)

	firstgene=False
	firstoperon=False
	overlapping={}
	firsts={}

	for fid,ft in enumerate(ofeats):
		frstart,frend=ftpos(ft,flen)
		dist = frstart-alen
		Compound = True if type(ft.location).__name__=="CompoundLocation" else False

		if ft.type=='gene':
			isgene=True
		if ft.type=='Operon':
			isoperon=True
		if dist<0:
			overlapping[fid]=ft
		if ft.type=='gene' and not firstgene:
			firstgene=True
			firsts[fid]=ft
		if ft.type=='Operon' and not firstoperon:
			firstoperon=True
			firsts[fid]=ft

		features[fid] = {hd:val for hd,val in IT.izip(ftheader,[ft.type,dist,frstart,frend,Compound,grouping[fid]]) }

	return ofeats, features, isgene, isoperon, overlapping, firsts


def ftpos(ft,flen):
	fstart = ft.location._start if type(ft.location).__name__ == "FeatureLocation" else ft.location.nofuzzy_start
	fend = ft.location._end if type(ft.location).__name__ == "FeatureLocation" else ft.location.nofuzzy_end

	frstart = fstart if strand == 1 else flen - fend
	frend = fend if strand == 1 else flen - fstart
	return frstart,frend


def groupfeatures(fmap,thres):
	grouping={}
	group=1
	for i in fmap.keys()[1:]:
		if fmap[i]-fmap[i-1]<thres:
			grouping[i]= group
			grouping[i-1] = group
			group+=1
	left=[k for k in fmap.keys() if k not in grouping.keys()]
	for i in left:
		grouping[i] = group
		group += 1

	return grouping











os.chdir("/Users/Povilas/Dropbox/Projects/2015-Metformin/Annotations/Ecoli/UAL")

genome_file="../Genome/Ecoli_K12_MG1655_U00096.3.gb"

info_file="Joined.csv"
genome=SeqIO.read(genome_file, "genbank")
data=maptable(info_file)



# for id,f in enumerate(genome.features[1:]):
# 	print id,f.type, f.location.start, f.location.end, f.location.strand, ','.join(f.qualifiers['gene']) if f.type == 'gene' else ''


# Collect primers
# Create input file of suitable primers to search ... suitably numbered.
# Note the key is constructed to be a suitable BLAST primer sequence to search.
#



#genome_file = "../Genome/Ecoli_K12_MG1655_U00096.3.gb"#.format(gen)
#genome = SeqIO.read(genome_file, "genbank")



gene_name_index=index_genbank_features(genome,'gene','gene')
CDS_name_index=index_genbank_features(genome,'CDS','gene')
gene_synonym_index=index_genbank_features(genome,'gene','gene_synonym')
CDS_synonym_index=index_genbank_features(genome,'CDS','gene_synonym')



len(gene_name_index.keys())
len(CDS_name_index.keys())
len(gene_synonym_index.keys())
len(CDS_synonym_index.keys())




#genome.features[CDS_name_index['pyrA']]





blast_leader="UAL_primers"

blast_input_file = blast_leader+ ".blast_input"

primer_list = []

primer_id = 0
idlist=[]
namelist=[]

names_gene=[]
names_CDS=[]
names_genesyn=[]
names_CDSsyn=[]
names_Manual=[]

names_UAL=[]

strandannot={}

fh = open(blast_input_file, 'w')

for key in data.keys():
	print key
	name=data[key]['Gene_Name']
	prim1 = Seq(data[key]['Primer.1'], generic_dna)
	prim2 = Seq(data[key]['Primer.2'], generic_dna)
	#id='{}:{}'.format(key,name)
	if prim1!='' or prim2!='':
		if data[key]['Strand_manual']!='NA':
			strand=data[key]['Strand_manual']
			data[key]['Strand_annotation'] = 'Manual'
			names_Manual.append(key)
		elif name in gene_name_index.keys():
			strand=genome.features[gene_name_index[name]].strand
			data[key]['Strand_annotation']='gene'
			names_gene.append(key)
		elif name in CDS_name_index.keys():
			strand = genome.features[CDS_name_index[name]].strand
			data[key]['Strand_annotation']='CDS'
			names_CDS.append(key)
		elif name in gene_synonym_index.keys():
			strand=genome.features[gene_synonym_index[name]].strand
			data[key]['Strand_annotation'] = 'genesyn'
			names_genesyn.append(id)
		elif name in CDS_synonym_index.keys():
			strand = genome.features[CDS_synonym_index[name]].strand
			data[key]['Strand_annotation'] = 'CDSsyn'
			names_CDSsyn.append(key)
		else:
			strand=data[key]['Strand']
			data[key]['Strand_annotation']='UAL'
			names_UAL.append(key)
		#Set primer orientation based on strand
		if strand==-1:
			query = prim2.reverse_complement() + 'N' * 20 + prim1.reverse_complement()
		else:
			query = prim1 + 'N' * 20 + prim2.reverse_complement()

		seq_record = SeqRecord(query, id='{}:{}:{}'.format(key,name,strand), description='')
		SeqIO.write(seq_record, fh, "fasta")
		primer_list.append(seq_record)
		namelist.append(name)
		idlist.append(key)
		primer_id += 1
fh.close()



for key in names_Manual:#names_UAL:
	print key, data[key]['Gene_Name'],data[key]['Strand'],data[key]['Strand_manual'],data[key]['Strand_annotation']



len(names_Manual)

len(namelist)

anlist=[nm for nm in namelist if nm in gene_name_index.keys()]

unanlist=[nm for nm in namelist if not nm in gene_name_index.keys()]

len(anlist)
len(unanlist)

unanlist


len(data.keys())

len(namelist)

len(idlist)


#Strand annotated

len(names_gene)+len(names_genesyn)
(len(names_gene)+len(names_genesyn))*100/1860

#Strand annotated by main gene name

1594*100/1860


#Unique amplicons
1811*100/1860

len(names_Manual)*100/1860

len(names_gene)
len(names_genesyn)
len(names_Manual)
len(names_UAL)




gen=3
ECgenome='Ecoli_K12_MG1655_U00096.{}'.format(gen)
# Make Blast DB
leader_filename='../Genome/'+ECgenome

blast_database_file = leader_filename + ".fasta"

#os.system("makeblastdb -in %s -dbtype nucl -title %s_BLAST_DB -out %s_BLAST_DB" % (blast_database_file, leader_filename, leader_filename))

#
# Do BLAST via BioPython with reading in XML format (although file is larger than necessary in this case ... does make it more portable).
#

blast_output_file = blast_leader +'_'+ ECgenome+ ".blast_output"

blastn_cline = NcbiblastnCommandline(query=blast_input_file, db=leader_filename + "_BLAST_DB", task='blastn-short',
                                     evalue='0.1', outfmt=5, out=blast_output_file)

print(blastn_cline)
# Should not really produce standard output or error.
stdout, stderr = blastn_cline()

result_handle = open(blast_output_file)
blast_records = NCBIXML.parse(result_handle)

recordmap={ record.query.split(':')[0] : record for record in list(blast_records)}


E_VALUE_THRESH = 0.05
ampmin=50#?
ampmax=1500


header=['Plate','Well','Gene_Name','Strand','Amplicons','Amplicon_length','Mapping', \
        'Primer1_e','Primer1_query_start','Primer1_query_end','Primer1_target_start','Primer1_target_end', \
        'Primer2_e','Primer2_query_start', 'Primer2_query_end', 'Primer2_target_start', 'Primer2_target_end']
table=[]
table.append(header)

for key in recordmap.keys(): #['AZ09-B11']: #
	record=recordmap[key]
	print record.query
	index,name,strand=record.query.split(':')
	plate,well=index.split('-')
	for alignment in record.alignments:
		hspsmap={ str(hsp.sbjct_start) : hsp for hsp in alignment.hsps if hsp.expect < E_VALUE_THRESH }
		ordkeys=sorted(hspsmap.keys())

		for ind, hkey in enumerate(ordkeys):
			hsp=hspsmap[hkey]
			print('****Alignment:{}****'.format(ind))
			#print('sequence:', alignment.title)
			#print('length:', alignment.length)
			print 'e value: {}, Strand: {}'.format(hsp.expect,hsp.frame)
			print("{}:{}".format(hsp.query_start,hsp.query_end))
			print(hsp.query)
			print(hsp.match)
			print(hsp.sbjct)
			print("{}:{}".format(hsp.sbjct_start,hsp.sbjct_end))
		print 'Possible alignments: {}'.format(len(ordkeys))
		if len(ordkeys)>0:
			amplicons=findamplicons(hspsmap,strand,ampmin,ampmax)
			# minkey=sorted(amplicons.keys())[0]
			#
			# prim1_hsp = amplicons[minkey]['Primer1'] #hspsmap[ordkeys[int(prim1_id)]]
			# prim2_hsp = amplicons[minkey]['Primer2'] #hspsmap[ordkeys[int(prim2_id)]]
			print 'Possible amplicons: {}'.format(len(amplicons.keys()))
			amptable=printamplicons(amplicons,strand)
			for ln in amptable:
				table.append([plate,well,name,strand]+ln)

csvwriter(table,'Amplicons_{}.csv'.format(ECgenome))




#Map amplicons to genome
#Now position coordinates of amplicons are oriented on forward strand

genome = SeqIO.read(genome_file, "genbank")
ti={hd:ind for ind, hd in enumerate(header)}
amplicons_features='{}_amplicons.gb'.format(ECgenome)
index_prev=''
amplicon=0
for rw in table[1:]:
	plate=rw[ti['Plate']]
	well=rw[ti['Well']]
	gene_name=rw[ti['Gene_Name']]
	index='{}-{}:{}'.format(plate,well,gene_name)
	key='{}-{}'.format(plate,well)
	amplicons=rw[ti['Amplicons']]
	#Only perfect matches
	if amplicons==1:
		if data[key]['Strand_annotation']=='UAL':
			col='255 0 0'
		else:
			col='152 251 152'

		if index!=index_prev:
			index_prev=index
			amplicon=1
		else:
			amplicon = amplicon+1
		indexa='{}:Amp-{}'.format(index,amplicon)
		if int(rw[ti['Amplicon_length']])>0:
			start=int(rw[ti['Primer1_target_start']])
			end=int(rw[ti['Primer2_target_end']])
			if end>start:
				strd=1
				truestart=start
				trueend=end
			else:
				strd=-1
				truestart=end
				trueend=start
			#print '{} Start:{}, End:{}'.format(indexa,start,end)
			seq_feature=SeqFeature(FeatureLocation(truestart,trueend, strand=strd), type="Amplicon", \
			                       id=indexa,\
			                       qualifiers={'Plate':plate,'Well':well,'Gene':gene_name,\
			                                   'Amplicon':amplicon,\
			                                   'Amplicons':amplicons,\
			                                   'Strand':strd,'colour':col,\
			                                   #'note':'{} {}'.format(indexa,data[key]['Description']),\
			                                   'Strand_annotation':data[key]['Strand_annotation']})
			genome.features.append(seq_feature)



faf = open(amplicons_features, 'w')
SeqIO.write(genome, faf, "gb")
faf.close()




#Map TF binding sites

TF_file="../Transcription_Factors/TF_binding_clean.csv"


TF_data=readcsv(TF_file,delim=',')

genomeamp = SeqIO.read(amplicons_features, "genbank")
amplicons_TFBS='{}_amplicons_TFBS.gb'.format(ECgenome)


TFindex={itm:ind  for ind, itm in enumerate(TF_data[0])}

for ln in TF_data[1:]:
	TFBS_start = int(ln[TFindex['TFBS_start']])
	TFBS_end = int(ln[TFindex['TFBS_end']])
	TFBSid = ln[TFindex['TFBSid']]
	if TFBS_start!=0 and TFBS_end!=0:

		TFid=ln[TFindex['TFid']]
		TF_name=ln[TFindex['TF_name']]
		strd=1 if ln[TFindex['Strand']]=='forward' else -1

		print TFid, TF_name, strd
		#Only perfect matches
		col='255 0 0'
		seq_feature=SeqFeature(FeatureLocation(TFBS_start,TFBS_end, strand=strd), type="TFBS", \
		                       id=TFBSid,\
		                       qualifiers={'TF_name':TF_name, \
		                                   'ID': TFid, \
		                                   'TF_id': TFid, \
		                                   'TU': ln[TFindex['TU']], \
		                                   'Promoter_name': ln[TFindex['Promoter_name']],\
		                                   'Effect': ln[TFindex['Effect']], \
		                                   'Evidence': ln[TFindex['Evidence']],
		                                   'Confidence': ln[TFindex['Confidence']], \
		                                   'Strand':strd,\
		                                   'colour':col})
		genomeamp.features.append(seq_feature)

faf = open(amplicons_TFBS, 'w')
SeqIO.write(genomeamp, faf, "gb")
faf.close()




#Map promoter regions


Prom_file="../Transcription_Factors/Promoters.csv"
Prom_data=readcsv(Prom_file,delim=',')

genomeprom = SeqIO.read(amplicons_TFBS, "genbank")
amplicons_Proms='{}_amplicons_TFBS_Prom.gb'.format(ECgenome)


Promindex={itm:ind  for ind, itm in enumerate(Prom_data[0])}

for ln in Prom_data[1:]:
	Promid = ln[Promindex['Promid']]
	Promoter_name = ln[Promindex['Promoter_name']]
	TSS = int(ln[Promindex['TSS']])
	Prom_sequence=ln[Promindex['Promoter_sequence']]
	Sigma = ln[Promindex['Sigma_factor']]
	strd = 1 if ln[Promindex['Strand']] == 'forward' else -1

	Promlen=len(Prom_sequence)
	TSSid=getindices(Prom_sequence)[0]+1 if Promlen >0 else 0


	print Promid, Promoter_name, TSS, TSSid, Promlen

	Prom_start = TSS-TSSid if strd==1 else TSS+TSSid
	Prom_end = TSS+(Promlen-TSSid) if strd==1 else TSS-(Promlen-TSSid)

	if TSS > 0 and Promlen > 0:
		#Only perfect matches
		col='0 0 255'
		seq_feature=SeqFeature(FeatureLocation(Prom_start,Prom_end, strand=strd), type="Promoter", \
		                       id=Promid,\
		                       qualifiers={'Promoter_name':Promoter_name, \
		                                   'ID': Promid,\
		                                   'Promoter_id': Promid, \
		                                   'TSS': TSS, \
		                                   'Evidence': ln[Promindex['Evidence']],
		                                   'Confidence': ln[Promindex['Confidence']], \
		                                   'Strand':strd,\
		                                   'colour':col})
		genomeprom.features.append(seq_feature)

faf = open(amplicons_Proms, 'w')
SeqIO.write(genomeprom, faf, "gb")
faf.close()





#Map operons

Oper_file="../Transcription_Factors/Operons.csv"
Oper_data=readcsv(Oper_file,delim=',')

genomeoperon = SeqIO.read(amplicons_Proms, "genbank")
amplicons_Operons='{}_amplicons_TFBS_Prom_Operon.gb'.format(ECgenome)

#Header index
Operindex={itm:ind  for ind, itm in enumerate(Oper_data[0])}

col='255 153 0'

for ln in Oper_data[1:]:
	Operon = ln[Operindex['Operon']]
	#In pyhton indexing starts at 0
	Oper_start = int(ln[Operindex['Start']])-1
	Oper_end = int(ln[Operindex['End']])-1
	strd = 1 if ln[Operindex['Direction']] == 'forward' else -1
	Oper_genes=ln[Operindex['Genes']]


	print Operon, Oper_start, Oper_end, strd, Oper_genes

	# if strd==1:
	# 	truestart=Oper_start
	# 	trueend=Oper_end
	# else:
	# 	truestart=Oper_end
	# 	trueend=Oper_start
	# if TSS > 0 and Promlen > 0:
	# 	#Only perfect matches

	seq_feature=SeqFeature(FeatureLocation(Oper_start,Oper_end, strand=strd), type="Operon", \
	                       id=Operon,\
	                       qualifiers={'Operon_name':Operon, \
	                                   'Genes': Oper_genes,\
	                                   'Gene_no': int(ln[Operindex['Genes_no']]),\
	                                   'Evidence': ln[Operindex['Evidence']],
	                                   'Confidence': ln[Operindex['Confidence']], \
	                                   'Strand':strd,\
	                                   'colour':col})
	genomeoperon.features.append(seq_feature)


faf = open(amplicons_Operons, 'w')
SeqIO.write(genomeoperon, faf, "gb")
faf.close()

#When GB files is saved, features become list?



#Something wrong with features that are read....
#genomeoperon=SeqIO.read("Ecoli_K12_MG1655_U00096.3_amplicons_TFBS_Prom_Operon.gb", "genbank")






#Find nearest features

expansion=5000

amplicons_feat=[ind for ind, itm in enumerate(genomeoperon.features) if itm.type=='Amplicon']


amplicons_map={ str(itm.qualifiers['Gene'][0]) : itm for itm in genomeoperon.features if itm.type=='Amplicon'}


operons_map={ str(itm.qualifiers['Operon_name']) : itm for itm in genomeoperon.features if itm.type=='Operon'}

#operons_map2={ str(itm.qualifiers['Operon_name'][0]) : itm for itm in genomeoperon2.features if itm.type=='Operon'}



amplicons_map['ribE']

ampl=amplicons_map['yagE']


opr=operons_map['yagEF']


opr.location._end-opr.location._start



astart=ampl.location._start
aend=ampl.location._end

aend-astart

[ ft for ft in genomeoperon[astart:aend+expansion].features if ft.strand==ampl.strand]



header=['UAL_promoter','Plate','Well','Strand','Type','Distance','Group','Name','Genes']
mappingdata=[]
mappingdata.append(header)


for ampfid in amplicons_feat: #[9799]:#
	amplicon=genomeoperon.features[ampfid]
	strand=amplicon.strand
	agene=amplicon.qualifiers['Gene'][0]
	aplate=amplicon.qualifiers['Plate'][0]
	awell=amplicon.qualifiers['Well'][0]

	astart=amplicon.location._start
	aend=amplicon.location._end
	alen = aend - astart

	extend=expansion
	while 1:
		fragment = genomeoperon[astart:aend+extend] if strand == 1 else genomeoperon[astart-extend:aend]
		feats, features, isgene, isoperon, overlapping, firsts=orderfeatures(fragment,strand,alen)
		if len(feats)>2 and isgene and isoperon:
			break
		extend=extend+2000
		if extend>100000:
			break


	flen=fragment.__len__()
	ftypes=[ft.type for ft in feats]

	print '\t', strand, alen, amplicon.qualifiers['Gene'],amplicon.qualifiers['Plate'], amplicon.qualifiers['Well'], len(feats), extend, ftypes
	#print features

	if len(overlapping)>0:
		selfeats=overlapping
	else:
		selfeats=firsts

	for fid,ft in selfeats.iteritems():
		finfo=features[fid]

		if ft.type=='Operon':
			ftname=ft.qualifiers['Operon_name']
			ftgenes=ft.qualifiers['Genes']
		else:
			ftname=ft.qualifiers['gene'][0]
			ftgenes=ftname

		#if finfo['Distance']<0:
		row=[agene,aplate,awell, str(strand), ft.type,finfo['Distance'],finfo['Group'], ftname, ftgenes]
		mappingdata.append(row)




csvwriter(mappingdata, 'Gene_Operon_mapping_{}.csv'.format(ECgenome))




fmap={fid:ftpos(ft,flen)[0]-alen for fid,ft in enumerate(feats)}



#How to deal with compound locations
ft=feats[2]


feats[1].location

type(feats[2].location).__name__=="CompoundLocation"




#For genes on reverse strand positions are not reversed

#Find downstream operon
#Search downstream till the first operon
#Check distance between amplicon end and operon start


amplicons_map['yncA'].location._start
amplicons_map['yncA'].location._end



operons_map['mnaT-ydcZ'].location._start
operons_map['mnaT-ydcZ'].location._end








len(genome.features)


len(genomeprom.features)


genome[184013:184309].features



genome[184300:184114].features



genome[182220:182490].features()


#Test whether amplicons overlap with intended genes









#
#
# for key in data.keys():
# 	print key
# 	name=data[key]['Gene_Name']
# 	prim1 = Seq(data[key]['Primer.1'], generic_dna)
# 	prim2 = Seq(data[key]['Primer.2'], generic_dna)
# 	strand=data[key]['Strand']
# 	if prim1!='' or prim2!='':
# 		if strand=='-1':
# 			query = prim2.reverse_complement() + 'N' * 20 + prim1.reverse_complement()
# 		else:
# 			query = prim1 + 'N' * 20 + prim2.reverse_complement()
# 		seq_record = SeqRecord(query, id='{}:{}:{}'.format(key,name,strand), description='')
# 		SeqIO.write(seq_record, fh, "fasta")
# 		primer_list.append(seq_record)
# 		idlist.append(key)
# 		primer_id += 1
#
# fh.close()
#






#Rubbish
#
# [[hspsmap[key].sbjct_start,hspsmap[key].sbjct_end,hspsmap[key].query_start,hspsmap[key].query_end] for key in sorted(hspsmap.keys())]
#
# for hsp in alignment.hsps:
# 	print hsp.sbjct_start
#
#
# # for seq_record in SeqIO.parse(genome_file, "genbank"):
# #     print(seq_record.id)
# #     print(repr(seq_record.seq))
# #     print(len(seq_record))
#
#
#
#
#
# genome.features[0].qualifiers
#
#
#
#
# gf=genome.features[10]
#
#
#
#
# genome.seq[7397:8423].reverse_complement()
#
#
# start = gf.location.nofuzzy_start
#
# end = gf.location.nofuzzy_end
#
#
# genome.seq[7397:8423].reverse_complement().translate(table=11,to_stop=True)
#
#
# #Indexing function example
# index_genbank_features(genome,"CDS","locus_tag")
#
# len(genome)
#
# genome.format('fasta')
#
# selection=genome[1000:2000]
#
#
# print selection.features[1]
#
#
# SeqFeature(FeatureLocation(9, 108, strand=-1), type="promoter")
#
#
#
#
#
#
#
#
# prim1 = Seq("CCGCTCGAGTCAGATCCCGTGGATTAACA", generic_dna)
# prim2 = Seq("CGGGATCCACGCGCAGTACGGAGTTC", generic_dna)
#
#
# query=prim1+'N'*20+prim2.reverse_complement()
#
#
# str(query)
#
# matrix = matlist.blosum62
# gap_open = -10
# gap_extend = -0.5
#
# alns = pairwise2.align.globalds(genome.seq,query, matrix, gap_open, gap_extend)
#
# top_aln = alns[0]
# top_aln
# aln_1, aln_2, score, begin, end = top_aln
#
#
#
#

