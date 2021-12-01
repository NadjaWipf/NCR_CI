from collections import defaultdict
import math
import sys

COMPLEMENT = {'A':'T','T':'A','G':'C','C':'G'}

TRA = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
       "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def readVCFnaive(fileName):
	S = defaultdict(dict)
	with open(fileName,'r') as IN:
		for l in IN:
			if l.startswith('#'):
				continue
			sl = l.split()
			chrom = sl[0]
			pos = sl[1]
			ref = sl[3]
			alt = sl[4].split(',')
			S[chrom][int(pos)] = (ref, alt )
	return S


def readFasta( fileName ):
	S={}
	with open(fileName,'r') as IN:
		seqName = ""
		seq = ""

		for l in IN:
			if l.startswith('>'):
				if not seqName == "":
					S[seqName] = seq
				seqName = l[1:].split()[0]
				seq=''

			else:
				seq += l.strip()
	return S

def readFastaFor1tr( fileName , trName):
	
	seq = ""
	with open(fileName,'r') as IN:
		seqName = ""
		store = False
		for l in IN:
			if l.startswith('>'):
				if store :
					break

				seqName = l[1:].split()[0]
				if seqName == trName:
					store=True

			elif store:
				seq += l.strip()
	return seq


transcriptFastaFile="Anopheles-coluzzii-Ngousso_TRANSCRIPTS_AcolN1.0.fa"
gtfFile="Anopheles-coluzzii-Ngousso_BASEFEATURES_AcolN1.0.gtf"
transcript=sys.argv[1] # 'ACON004707-RD'

transcriptSTR = '"' + transcript + '";'

SNPFile = sys.argv[2]

OUTfile = sys.argv[3]

## reading GTF
chromosom = None
exonList = []
geneReverseStrand = False

with open(gtfFile) as IN :

	for l in IN: 
		if l.startswith('#'):
			continue
		sl = l.strip().split()
		
		if sl[2] != "exon" :
			continue

		if sl[11] == transcriptSTR:
			
			chromosom = sl[0]

			exonStop = int(sl[4])
			exonStart = int(sl[3])

			geneReverseStrand = sl[6] == '-'
			if not geneReverseStrand:
				exonList.append( (exonStart,exonStop) )
			else:
				exonList.insert(0, (exonStart,exonStop) )


# we sort exons by their order of appearance on the transcript 
exonList = sorted(exonList , reverse= geneReverseStrand)

a,b = zip(exonList[0] , exonList[-1])
trSTART = min(a)
trSTOP = max(b)

## reading VCF

SNPs = readVCFnaive(SNPFile)
#print(SNPs)
DNApositions=[]

#filtering SNPs which are in the boundaries
for pos in SNPs[chromosom]:
	if pos <= trSTOP and pos >= trSTART:
		DNApositions.append( pos )

print("read",len(DNApositions),"SNP positions in gene boundaries")
## DNA pos to RNA pos

# we sort SNPs by their order of appearance on the transcript 
DNApositions = sorted(DNApositions , reverse= geneReverseStrand)

DNAtoRNApositions={}

curPos = 1
for exon in exonList:
	#print(curPos,exon)

	if geneReverseStrand : # reverse strand

		while DNApositions[0] > exon[0] : # the current SNP is before the (reversed) end of this exon

			if DNApositions[0] > exon[1]: # the current SNP is before the (reversed) beginning of this exon
				#print( DNApositions[0], "intronic" )
				DNAtoRNApositions[DNApositions[0]] = None
			else:
				RNApos = curPos + exon[1] - DNApositions[0]

				#print( DNApositions[0] , RNApos , math.ceil(RNApos/3) )
				DNAtoRNApositions[DNApositions[0]] = RNApos
			
			DNApositions.pop(0)

			if len( DNApositions )==0:
					break
	else : # forward strand
		while DNApositions[0] < exon[1] : # the current SNP is before the end of this exon

			if DNApositions[0] < exon[0]: # the current SNP is before the (reversed) beginning of this exon
				#print( DNApositions[0], "intronic" )
				DNAtoRNApositions[DNApositions[0]] = None
			else :
				RNApos = curPos + DNApositions[0] - exon[0]

				#print( DNApositions[0] , RNApos , math.ceil(RNApos/3) )
				DNAtoRNApositions[DNApositions[0]] = RNApos

			DNApositions.pop(0)

			if len( DNApositions )==0:
					break

	curPos += (exon[1] - exon[0])+1 # add exon length to current position
	if len( DNApositions )==0:
		break

print("assigned", len(DNAtoRNApositions) , "positions along the transcript")
## reading fasta

transcriptSeq = readFastaFor1tr( transcriptFastaFile , transcript)
#print( transcript )
#print( transcriptSeq[:10] )


## finding codon and if synonymous or not
DNApostoAltNames = {}

for DNApos in DNAtoRNApositions.keys() :

	RNApos = DNAtoRNApositions[DNApos]
	if RNApos is None :
		#print("intronic position")
		altNames = ['intronic']*len(SNPs[chromosom][DNApos][1])
		#print(chromosom,DNApos, ','.join(altNames) , SNPs[chromosom][DNApos][0] , ','.join( SNPs[chromosom][DNApos][1] )  , sep='\t')
		DNApostoAltNames[DNApos] = altNames
		continue
	codon = math.ceil(RNApos/3)

	codonPosition = (RNApos-1)%3 # 0, 1 or 2

	# NB : additional offset of -1 because of 1-based indexing
	codonSeq = transcriptSeq[ codon*3 -3 : codon*3  ]
	
	#print(transcript , DNApos , ': codon' , codon )
	#print( '\tref:', codonSeq , TRA[codonSeq]) 
	## checking that vcf and fasta reference match (could catch messing up position conversions)
	refFromVCF = SNPs[chromosom][DNApos][0]
	refFromFasta = codonSeq[codonPosition]
	if geneReverseStrand:
		refFromFasta = COMPLEMENT[refFromFasta]

	if not refFromFasta == refFromVCF:
		print(transcript , DNApos , ': codon' , codon ,"!!diff fasta/vcf : ",refFromFasta,refFromVCF )

	altNames = []
	for alt in SNPs[chromosom][DNApos][1]:
		if alt == '*':
			#print('\t',alt,':','indel' , 'indel')
			name='indel'
		else:
			x = list(codonSeq)
			x[ codonPosition ] = alt
			codonAlt = ''.join(x)
			#print('\t',alt,':',codonAlt , TRA[codonAlt])
			if TRA[codonAlt] == TRA[codonSeq]:
				name = "synonymous"
			else:
				name = TRA[codonSeq] + str(codon) + TRA[codonAlt]
		altNames.append(name)

	#print(chromosom,DNApos, ','.join(altNames) , refFromVCF , ','.join( SNPs[chromosom][DNApos][1] )  , sep='\t')
	DNApostoAltNames[DNApos] = altNames

print('assigned name for',len(DNApostoAltNames),'positions')

S = defaultdict(dict)
with open(SNPFile,'r') as IN , open(OUTfile,'w') as OUT:
	for l in IN:
		if not l.startswith('#'):
			
			sl = l.split('\t')
			pos = int(sl[1])
			sl[2] = ','.join(DNApostoAltNames[pos] )
			l = '\t'.join(sl)
		print(l , end='',file=OUT)



