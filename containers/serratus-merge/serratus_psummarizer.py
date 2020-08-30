#!/usr/bin/python2

import sys
import os

SRA = os.getenv("SUMZER_SRA", None)
MAXALNS = int(os.getenv("SUMZER_MAXALNS", "10000000"))
MAXX = int(os.getenv("SUMZER_MAXX", "100"))
THROWX = (os.getenv("SUMZER_THROWX", "YES") != "NO")
SUMZER_COMMENT = os.getenv("SUMZER_COMMENT", None)

CVG_BINS = 25
MINCVGPCT = 50

ReportFileName = sys.argv[1]
EchoFileName = None
if len(sys.argv) > 2:
	EchoFileName = sys.argv[2]

fRep = open(ReportFileName, "w")
fEcho = None
if len(sys.argv) > 2:
	fEcho = open(EchoFileName, "w")
ExceptCount = 0

'''
#     ReadLabel               RefLabel  QLo     QHi      QL      TLo     THi      TL   PctId     Evalue
#             0                      1    2       3       4        5       6       7       8          9
ERR2756788.4545 Coronaviridae.repl1a.1  152     298     299     1230    1278    3990    57.1    1.1e-11
ERR2756788.4927 Adenoviridae.Q9WF13.1   238     158     301     97      123     485     63.0    3.0e-06
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

def CountToSymbol(i):
	#		    012345678901  2^11=2048
	Symbols = "_.:uwaomUWAOM^"
	n = 1
	for c in Symbols:
		if i < n:
			return c
		n *= 2
	return "^"

SumReadLength = 0
AlnCount = 0

IdToAlnCount = {}
IdToSumPctId = {}
IdToSumBasesPctId = {}
IdToSumBases = {}
IdToCoverageVec = {}

Genes = []
Families = []

AccToLength = {}

def IncDict(D, k, n=1):
	try:
		D[k] += n
	except:
		D[k] = n

def IncCvg(Id, TLo, THi, TL):
	global IdToCoverageVec

	try:
		v = IdToCoverageVec[Id]
	except:
		IdToCoverageVec[Id] = [ 0 ]*CVG_BINS

	Mid = (TLo + THi)//2
	Bin = ((Mid-1)*CVG_BINS)//TL
	if Bin < 0:
		Bin = 0
	if Bin >= CVG_BINS:
		Bin = CVG_BINS - 1
	IdToCoverageVec[Id][Bin] += 1

def MakeCartoon(Id):
	try:
		w = IdToCoverageVec[Id]
	except:
		return "~"

	s = ""
	for i in w:
		s += CountToSymbol(i)
	return s 

def GetDict(Dict, Key, Default = None):
	try:
		return Dict[Key]
	except:
		return Default

def GetCvgBins(Id, MinCount):
	try:
		v = IdToCoverageVec[Id]
	except:
		return 0

	N = len(v)
	if N == 0:
		return 0
	n = 0
	for x in v:
		if x >= MinCount:
			n += 1
	return n

def GetPctId(Id):
	SumBases = GetDict(IdToSumBases, Id)
	SumBasesPctId = GetDict(IdToSumBasesPctId, Id)
	PctId = 0
	if SumBases > 0:
		PctId = float(SumBasesPctId)/SumBases
	if PctId < 50:
		PctId = 50
	if PctId > 100:
		PctId = 100
	return PctId

def GetIdentityWeight(PctId):
	FractId = PctId/100.0
	if FractId < 0.5:
		FractId = 0.5
	if FractId > 1.0:
		FractId = 1.0
	Weight = 1.0/FractId
	return Weight

def GetScore(Id):
	Cvg1 = GetCvgBins(Id, 1)
	Cvg8 = GetCvgBins(Id, 8)
	RawScore = Cvg8*4 + (Cvg1 - Cvg8)
	PctId = GetPctId(Id)
	Weight = GetIdentityWeight(PctId)
	Score = RawScore*Weight
	if Score > 100:
		Score = 100
	return Score

def GetGene(Label):
	Fields = Label.split('.')
	if len(Fields) == 1:
		return Label
	Gene = Fields[0] + "." + Fields[1]
	return Gene

def GetFamily(Label):
	Family = Label.split('.')[0]
	return Family

def DoAln_():
	global SumReadLength

	Fields = Line.split()
	try:
		RefLabel = Fields[1]
		QL = int(Fields[4])
		TLo = int(Fields[5])
		THi = int(Fields[6])
		TL = int(Fields[7])
		PctId = float(Fields[8])
	except:
		return

	Acc = RefLabel
	Family = GetFamily(RefLabel)
	Gene = GetGene(RefLabel)
	if Gene not in Genes:
		Genes.append(Gene)
	if Family not in Families:
		Families.append(Family)

#	n = RefLabel.find('.')
#	assert n > 0
#	Family = RefLabel[:n]

	SumReadLength += QL

	if THi > TLo:
		AlnLength = THi - TLo + 1
	else:
		AlnLength = TLo - THi + 1

	AccToLength[Acc] = TL
#	AddFamilyAcc(Family, Acc)

	for Id in [ Acc, Gene, Family ]:
		IncCvg(Id, TLo, THi, TL)
		IncDict(IdToAlnCount, Id)
		IncDict(IdToSumBases, Id, AlnLength)
		IncDict(IdToSumPctId, Id, PctId)
		IncDict(IdToSumBasesPctId, Id, PctId*AlnLength)

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

def Rep(s):
	global SRA
	if SRA != None:
		fRep.write("sra=" + SRA + ";")
	fRep.write(s + "\n")

def GetLine(Id):
	global Genes

	Cartoon = MakeCartoon(Id)
	PctId = GetPctId(Id)
	Score = GetScore(Id)
	Alns = GetDict(IdToAlnCount, Id)
	SumBases = GetDict(IdToSumBases, Id)
	IsGene = (Id in Genes)
	IsFamily = (Id in Families)

	AvgCols = 0
	if Alns > 0:
		AvgCols = float(SumBases)/Alns

	s = ""
	if IsGene:
		s += "gencvg=%s;" % Cartoon
		s += "gen=%s;" % Id
	elif IsFamily:
		s += "famcvg=%s;" % Cartoon
		s += "fam=%s;" % Id
	else:
		s += "seqcvg=%s;" % Cartoon
		s += "seq=%s;" % Id
	s += "score=%.0f;" % Score
	s += "pctid=%.0f;" % PctId
	s += "alns=%.0f;" % Alns
	s += "avgcols=%.0f;" % AvgCols
	return s

def Report():
	global AlnCount
	global IdToScore

	MeanReadLength = 0
	if AlnCount > 0:
		MeanReadLength = SumReadLength/AlnCount

	s = ""
	if SUMZER_COMMENT != None:
		s += "SUMZER_COMMENT=" + SUMZER_COMMENT + ";"
	s += "totalalns=%d;readlength=%d;" % (AlnCount, MeanReadLength)
	if ExceptCount > 0:
		s += "excepts=%d;" % ExceptCount
	if AlnCount == MAXALNS:
		s += "truncated=yes;"
	else:
		s += "truncated=no;"
	Rep(s)

	Ids = list(IdToAlnCount.keys())

	IdToScore = {}
	for Id in Ids:
		Score = GetScore(Id)
		IdToScore[Id] = Score
	Order = GetOrder(IdToScore)
	Ids = list(IdToScore.keys())

	for k in Order:
		Id = Ids[k]
		if not Id in Families:
			continue
		s = GetLine(Id)
		Rep(s)

	for k in Order:
		Id = Ids[k]
		if not Id in Genes:
			continue
		s = GetLine(Id)
		Rep(s)

	for k in Order:
		Id = Ids[k]
		if Id in Genes or Id in Families:
			continue
		s = GetLine(Id)
		Rep(s)

for Line in sys.stdin:
	if fEcho != None:
		fEcho.write(Line)
	DoAln()
	AlnCount += 1
	if AlnCount >= MAXALNS or ExceptCount >= MAXX:
		Report()
		sys.exit(0)

Report()
