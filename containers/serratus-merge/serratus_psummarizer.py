#!/usr/bin/python2

import sys
import os

SRA = os.getenv("SUMZER_SRA", "SRAxxxxxx")
MAXALNS = int(os.getenv("SUMZER_MAXALNS", "10000000"))
MAXX = int(os.getenv("SUMZER_MAXX", "100"))
THROWX = (os.getenv("SUMZER_THROWX", "YES") != "NO")
COMMENT = os.getenv("SUMZER_COMMENT", "")

CVG_BINS = 25
MINCVGPCT = 50

'''
#     ReadLabel               RefLabel  QLo     QHi      QL      TLo     THi      TL   PctId     Evalue
#             0                      1    2       3       4        5       6       7       8          9
ERR2756788.4545 Coronaviridae.repl1a.1  152     298     299     1230    1278    3990    57.1    1.1e-11
ERR2756788.4927 Adenoviridae.Q9WF13.1   238     158     301     97      123     485     63.0    3.0e-06
ERR2756788.5782 Coronaviridae.nucp.2    68      283     297     91      158     444     47.2    2.8e-12
ERR2756788.6146 Poxviridae.Q6VZU9.1     132     16      301     76      114     325     53.8    6.5e-09
ERR2756788.6146 Poxviridae.Q6VZU9.1     152     301     301     50      99      325     46.0    1.9e-08
ERR2756788.6232 Coronaviridae.repl1a.1  156     1       302     2980    3031    4117    73.1    7.4e-21
'''

# Include bitscore, gaps and mismatches
'''
#        ReadLabel                   RefLabel  QLo     QHi       QL     TLo     THi      TL    PctId     Evalue BtSc  MiM     Gap
#                0                          1    2       3        4       5       6       7        8          9   10   11      12
ERR2756788.31345/1      Coronaviridae.R1AB.15   151     2       151     209     258     349     86.0    7.2e-22 96.3    7       0
ERR2756788.31345/2      Coronaviridae.R1AB.15   9       149     151     181     227     349     80.9    3.7e-18 84.0    9       0
ERR2756788.37873/1      Coronaviridae.R1AB.9    1       129     151     542     584     584     76.7    2.1e-13 68.2    10      0
ERR2756788.43195/1      Coronaviridae.R1AB.15   151     2       151     67      116     349     70.0    7.0e-17 79.7    15      0
ERR2756788.46979/1      Adenoviridae.CAPSH.1    129     1       151     1       43      498     48.8    1.6e-08 52.0    22      0
ERR2756788.47790/2      Coronaviridae.R1AB.15   151     2       151     59      108     349     64.0    2.5e-14 71.2    18      0
ERR2756788.49252/1      Coronaviridae.NCAP.2    148     11      151     75      120     337     45.7    4.8e-10 57.0    25      0
ERR2756788.59541/1      Coronaviridae.SPIKE.3   149     3       151     167     215     651     75.5    8.2e-18 82.8    12      0
ERR2756788.59541/2      Coronaviridae.SPIKE.3   89      148     149     162     181     651     85.0    2.3e-04 38.1    3       0
'''

def CmpKey__(i):
	global d, Keys

	ki = Keys[i]
	ni = d[ki]
	return ni

def GetOrder(Dict):
	global d, Keys

	d = Dict
	Keys = list(d.keys())
	N = len(Keys)
	Order = list(range(0, N))
	Order.sort(key=CmpKey__)
	Order.reverse()
	return Order

fRep = open(sys.argv[1], "w")
fEcho = None
if len(sys.argv) > 2:
	fEcho = open(sys.argv[2], "w")

AlnCount = 0
ExceptCount = 0

FamilyToAlnCount = {}
FamilyToPctIds = {}
FamilyToGenes = {}
FamilyToTotalAlnLength = {}

FamilyGeneToLength = {}
FamilyGeneToAlnCount = {}
FamilyGeneToPctIds = {}
FamilyGeneToBins = {}
FamilyGeneToTotalAlnLength = {}

SumReadLength = 0

def IncCount(D, k, n=1):
	try:
		D[k] += n
	except:
		D[k] = n

def AppendList(D, k, x):
	try:
		D[k].append(x)
	except:
		D[k] = [ x ]

def GetMedian(D, k):
	try:
		v = D[k]
	except:
		return -1
	v.sort()
	n = len(v)
	Median = v[n/2]
	return Median

def IncCvg(FamilyGene, TLo, THi, TL):
	try:
		v = FamilyGeneToBins[FamilyGene]
	except:
		FamilyGeneToBins[FamilyGene] = [ 0 ]*CVG_BINS

	Mid = (TLo + THi)/2
	Bin = ((Mid-1)*CVG_BINS)//TL
	if Bin < 0:
		Bin = 0
	if Bin >= CVG_BINS:
		Bin = CVG_BINS - 1
	
	FamilyGeneToBins[FamilyGene][Bin] += 1

def MakeCartoon(FamilyGene):
	try:
		w = FamilyGeneToBins[FamilyGene]
	except:
		return "?"

	Max = max(w)
	if Max < 4:
		Max = 4
	s = ""
	for i in w:
		if i == 0:
			s += '_'
		elif i <= Max/4:
			s += '.'
		elif i <= Max/2:
			s += 'o'
		else:
			s += 'O'
	return s 

def GetCvgPct(FamilyGene):
	try:
		w = FamilyGeneToBins[FamilyGene]
	except:
		return 0

	N = len(w)
	n = 0
	for i in w:
		if i > 2:
			n += 1
	Pct = (n*100)//N
	return Pct

def AddFamilyGene(Family, Gene):
	global FamilyToGenes
	try:
		v = FamilyToGenes[Family]
	except:
		FamilyToGenes[Family] = [ Gene ]
		return

	if not Gene in FamilyToGenes[Family]:
		FamilyToGenes[Family].append(Gene)

def GetGeneCount(Family):
	try:
		Genes = FamilyToGenes[Family]
		return len(Genes)
	except:
		return 0

def GetCoveredGenePcts(Family):
	try:
		Genes = FamilyToGenes[Family]
	except:
		return 0

	Pcts = []
	for Gene in Genes:
		FamilyGene = Family + "." + Gene
		CvgPct = GetCvgPct(FamilyGene)
		if CvgPct >= MINCVGPCT:
			Pcts.append(CvgPct)
	return Pcts

def GetCoveredGenesCount(Family):
	Pcts = GetCoveredGenePcts(Family)
	n = len(Pcts)
	return n

def GetAlnCount(Family):
	try:
		AlnCount = FamilyToAlnCount[Family]
	except:
		AlnCount = 0
	return AlnCount

def GetTotalAlnLength(Family):
	try:
		n = FamilyToTotalAlnLength[Family]
	except:
		n = 0
	return n

def GetScore(Family):
	NA = GetAlnCount(Family)
	Pcts = GetCoveredGenePcts(Family)
	Pcts.sort()
	N = len(Pcts)
	TopPct = 0
	SecondPct = 0
	if N > 0:
		TopPct = Pcts[N-1]
	if N > 1:
		SecondPct = Pcts[N-2]
	Score = (TopPct + SecondPct)/4
#	print >> sys.stderr, Family, Pcts, Score, NA
	if NA > 5000:
		Score += 50
	elif NA > 2000:
		Score += 40
	elif NA > 1000:
		Score += 20
	elif NA > 100:
		Score += 10
	elif NA > 1:
		Score += 5

	return Score

def Rep(s):
	global SRA
	fRep.write(SRA + ";" + s + "\n")

def Report():
	global AlnCount

	s = "comment:%s;" % (COMMENT)
	Rep(s)

	MeanReadLength = 0
	if AlnCount > 0:
		MeanReadLength = SumReadLength/AlnCount
	s = "totals:alns=%d;readlength=%d;" % (AlnCount, MeanReadLength)
	if ExceptCount > 0:
		s += "excepts=%d;" % ExceptCount
	if AlnCount == MAXALNS:
		s += "truncated=yes;"
	else:
		s += "truncated=no;"
	Rep(s)

	Order = GetOrder(FamilyToAlnCount)
	Families = list(FamilyToAlnCount.keys())
	n = len(Families)
	for k in Order:
		Family = Families[k]
		AlnCount = FamilyToAlnCount[Family]
		TotalAlnLength = FamilyToTotalAlnLength[Family]
		AvgAlnlength = 0.0
		if AlnCount > 0:
			AvgAlnlength = float(TotalAlnLength)/AlnCount
		MedPctId = GetMedian(FamilyToPctIds, Family)
		GeneCount = GetGeneCount(Family)
		CoveredGeneCount = GetCoveredGenesCount(Family)
		Score = GetScore(Family)

		s = "family:%s;" % Family
		s += "score=%.0f;" % Score
		s += "alns=%d;"% AlnCount
		s += "avgcols=%.0f;" % AvgAlnlength
		s += "pctid=%.1f;" % MedPctId
		s += "genes=%d/%d;" % (CoveredGeneCount, GeneCount)
		Rep(s)

	Order = GetOrder(FamilyGeneToAlnCount)
	FGs = list(FamilyGeneToAlnCount.keys())
	n = len(FGs)
	for k in Order:
		FamilyGene = FGs[k]
		AlnCount = FamilyGeneToAlnCount[FamilyGene]
		MedPctId = GetMedian(FamilyGeneToPctIds, FamilyGene)
		Length = FamilyGeneToLength[FamilyGene]
		CvgPct = GetCvgPct(FamilyGene)
		Cartoon = MakeCartoon(FamilyGene)
		TotalAlnLength = FamilyGeneToTotalAlnLength[FamilyGene]
		AvgAlnlength = 0.0
		if AlnCount > 0:
			AvgAlnlength = float(TotalAlnLength)/AlnCount
		s = "gene:%s;" % FamilyGene
		s += "alns=%d;"% AlnCount
		s += "avgcols=%.0f;" % AvgAlnlength
		s += "medpctid=%.1f;" % MedPctId
		s += "len=%d;" % Length
		s += "cvgpct=%d;" % CvgPct
		s += "cvg=%s;" % Cartoon
		Rep(s)

def DoAln_():
	global SumReadLength

	Fields = Line.split()
	if len(Fields) < 10:
		print >> sys.stderr, len(Fields), Fields, Line
		assert False

	RefLabel = Fields[1]
	QLo = int(Fields[2])
	QHi = int(Fields[3])
	QL = int(Fields[4])
	TLo = int(Fields[5])
	THi = int(Fields[6])
	TL = int(Fields[7])
	PctId = float(Fields[8])
	# Evalue = float(Fields[9])
	# BtSc = float(Fields[10])
	# MiM = int(Fields[11])
	# Gap = int(Fields[12])

	Fields = RefLabel.split('.')
	assert len(Fields) == 3
	Family = Fields[0]
	Gene = Fields[1]

	FamilyGene = Family + "." + Gene

	SumReadLength += QL

	if THi > TLo:
		AlnLength = THi - TLo + 1
	else:
		AlnLength = TLo - THi + 1

	FamilyGeneToLength[FamilyGene] = TL
	AddFamilyGene(Family, Gene)
	IncCount(FamilyToAlnCount, Family)
	IncCount(FamilyGeneToAlnCount, FamilyGene)
	IncCount(FamilyToTotalAlnLength, Family, AlnLength)
	IncCount(FamilyGeneToTotalAlnLength, FamilyGene, AlnLength)
	AppendList(FamilyGeneToPctIds, FamilyGene, PctId)
	AppendList(FamilyToPctIds, Family, PctId)
	IncCvg(FamilyGene, TLo, THi, TL)

def DoAln():
	global ExceptCount

	if THROWX:
		DoAln_()
		return
	try:
		DoAln_()
	except:
		ExceptCount += 1
		return

for Line in sys.stdin:
	if fEcho != None:
		fEcho.write(Line)
	DoAln()
	AlnCount += 1
	if AlnCount >= MAXALNS or ExceptCount >= MAXX:
		Report()
		sys.exit(0)
Report()
