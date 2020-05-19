#!/usr/bin/python3

# Author: Robert C. Edgar
# email robert@drive5.com

import sys
import os
import re

# **************************************************
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# MUST set RAISEX to False for production
# If True, exceptions are raised, which is helpful
# for debugging but fatal for generating a BAM file.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# **************************************************
RAISEX = False

CVG_BINS = 32
COV_PAN_GENOME =    "Cov___Pan"
OTHER_PAN_GENOME = "Other_Pan"
COV_FAM = "Coronaviridae"

SAMInputFileName = sys.argv[1]
MetaFileName = sys.argv[2]
SummaryFileName = sys.argv[3]
SAMOutputFileName = sys.argv[4]

fIn = open(SAMInputFileName)
fSum = open(SummaryFileName, "w")
fOut = open(SAMOutputFileName, "w")
## fTinyhit  = open(TinyhitOutputFileName, "w")

AccToTax = {}
AccToName = {}
AccToLen = {}
AccToFam = {}

#   0     1    2      3
# Acc	Len	Name	Fam
for Line in open(MetaFileName):
	Fields = Line[:-1].split('\t')
	Acc = Fields[0]
	Len = int(Fields[1])
	Name = Fields[2]
	Fam = Fields[3]

	AccToName[Acc] = Name
	AccToLen[Acc] = Len
	AccToFam[Acc] = Fam

AccToSumBases = {}
AccToSumBasesPctId = {}
AccToCoverageVec = {}

AccToCoverageVec[COV_PAN_GENOME] = [0]*CVG_BINS
AccToCoverageVec[OTHER_PAN_GENOME] = [0]*CVG_BINS

CovSeqs = [ None ]*CVG_BINS
OtherSeqs = [ None ]*CVG_BINS

d = {}
Order = []
Keys = []

def CharToProb(c):
	ic = ord(c)
	iq = ic - 33
	if iq < 0:
		return 1.0
	if iq > 60:
		return 0.0
	return 10**(-iq/10.0)

def GetEE(Qual):
	SumP = 0
	L = len(Qual)
	if L == 0:
		return 0.0

	for q in Qual:
		P = CharToProb(q)
		SumP += P
	return SumP

def GetAlnLen(CIGAR):
	Ns = []
	Letters = []

	if CIGAR == "*":
		return 100
	try:
		Len = 0
		n = 0
		for c in CIGAR:
			if c.isdigit():
				n = n*10 + (ord(c) - ord('0'))
			elif c.isalpha() or c == '=':
				if c != 'S' and c != 'H':
					Len += n
				n = 0
	except:
		if RAISEX:
			raise
		Len = 100
	return Len

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

def MakeCartoon(v):
	Max = max(v)
	if Max < 4:
		Max = 4
	s = ""
	for i in v:
		if i == 0:
			s += '_'
		elif i <= Max/4:
			s += '.'
		elif i <= Max/2:
			s += 'o'
		else:
			s += 'O'
	return s 

def AllMatch(CIGAR):
	M = re.match("[0-9]+M", CIGAR)
	return M != None

AccToHits = {}

def AddHit(Acc, TBin, L, PctId):
	try:
		AccToHits[Acc] += 1
		AccToSumBases[Acc] += L
		AccToSumBasesPctId[Acc] += L*PctId
		AccToCoverageVec[Acc][TBin] += 1
	except:
		AccToHits[Acc] = 1
		AccToSumBases[Acc] = L
		AccToSumBasesPctId[Acc] = L*PctId
		AccToCoverageVec[Acc] = [0]*CVG_BINS
		AccToCoverageVec[Acc][TBin] += 1

CovMapped = 0
OtherMapped = 0
SumL = 0
MaxL = 0
for Line in fIn:
	fOut.write(Line)
	# Ignore SAM headers
	# if Line.startswith('@'):
		# continue

	# Wrap everything in a try-except block because
	# the summarizer MUST not crash the pipeline!
	# Better no summary than no BAM file.
	try:
		Fields = Line[:-1].split('\t')

		# ReadLabel = Fields[0]
		Flags = int(Fields[1])
		if (Flags & 0x4) == 0x4:
			continue

		Acc = Fields[2]
		try:
			Fam = AccToFam[Acc]
		except:
			Fam = None

		if Fam == COV_FAM:
			CovMapped += 1
		else:
			OtherMapped += 1

		TPos = int(Fields[3])
		# MAPQ = Fields[4]
		CIGAR = Fields[5]
		# RNEXT = Fields[6]
		# PNEXT = Fields[7]
		# TL = int(Fields[8])
		# SEQ = Fields[9]
		QUAL = Fields[10]
		L = GetAlnLen(CIGAR)
		SumL += L
		if L > MaxL:
			MaxL = L

		try:
			ee = GetEE(QUAL)
			EE = "%.2f" % ee
		except:
			if RAISEX:
				raise
			EE = "-1.0"

		AS = None
		NM = None
		for Field in Fields[11:]:
			if Field.startswith("NM:i:"):
				s = Field.replace("NM:i:", "")
				try:
					NM = int(s)
				except:
					if RAISEX:
						raise
					NM = None
				break

		## if fTinyhit != None:
			## print(Acc + "\t" + str(TPos) + "\t" + str(L) + "\t" + str(NM) + "\t" + EE, file=fTinyhit)

		PctId = 0
		if NM != None and L > 0:
			PctId = (L - NM)*100.0/L

		TL = None
		try:
			TL = AccToLen[Acc]
			TBin = ((TPos-1)*CVG_BINS)//TL
			if TBin < 0:
				TBin = 0
			if TBin >= CVG_BINS:
				TBin = CVG_BINS - 1
		except:
			TBin = 0

		if CIGAR.find("S") < 0:
			AddHit(Acc, TBin, L, PctId)
			if Fam == COV_FAM:
				AddHit(COV_PAN_GENOME, TBin, L, PctId)
			else:
				AddHit(OTHER_PAN_GENOME, TBin, L, PctId)

		if AllMatch(CIGAR):
			SEQ = Fields[9]
			if Fam == COV_FAM:
				if CovSeqs[TBin] == None:
					CovSeqs[TBin] = SEQ
			else:
				if OtherSeqs[TBin] == None:
					OtherSeqs[TBin] = SEQ

	except:
		if RAISEX:
			raise
		pass

def GetCovFract(v):
	N = len(v)
	if N == 0:
		return 0
	n = 0
	for x in v:
		if x > 0:
			n += 1
	return float(n)/N

Accs = list(AccToHits.keys())

AccToCovFract = {}
for Acc in Accs:
	CovFract = GetCovFract(AccToCoverageVec[Acc])
	AccToCovFract[Acc] = CovFract

Order = GetOrder(AccToCovFract)

def IsPan(Acc):
	return Acc.find("_Pan") >= 0

TopCovAcc = None
TopOtherAcc = None
for k in Order:
	Acc = Accs[k]
	try:
		Fam = AccToFam[Acc]
	except:
		Fam = None
	if Fam == None:
		continue

	if Fam == COV_FAM and TopCovAcc == None:
		TopCovAcc = Acc
	if Fam != COV_FAM and TopOtherAcc == None:
		TopOtherAcc = Acc

CovPanCartoon = MakeCartoon(AccToCoverageVec[COV_PAN_GENOME])
OtherPanCartoon = MakeCartoon(AccToCoverageVec[OTHER_PAN_GENOME])

CovPanCovFract = GetCovFract(AccToCoverageVec[COV_PAN_GENOME])
OtherPanCovFract = GetCovFract(AccToCoverageVec[OTHER_PAN_GENOME])

try:
	CovPanSumBases = AccToSumBases[COV_PAN_GENOME]
except:
	CovPanSumBases = 0

try:
	CovPanSumBasesPctId = AccToSumBasesPctId[COV_PAN_GENOME]
except:
	CovPanSumBasesPctId = 0

try:
	OtherPanSumBases = AccToSumBases[OTHER_PAN_GENOME]
except:
	OtherPanSumBases = 0

try:
	OtherPanSumBasesPctId = AccToSumBasesPctId[OTHER_PAN_GENOME]
except:
	OtherPanSumBasesPctId = 0

CovPanCovPct = CovPanCovFract*100.0
CovPanPctId = 0.0
if CovPanSumBasesPctId > 0:
	CovPanPctId = CovPanSumBasesPctId/CovPanSumBases

OtherPanCovPct = OtherPanCovFract*100.0
OtherPanPctId = 0.0
if OtherPanSumBasesPctId > 0:
	OtherPanPctId = OtherPanSumBasesPctId/OtherPanSumBases

try:
	TopCovName = AccToName[TopCovAcc]
except:
	TopCovName = "?"

try:
	TopOtherName = AccToName[TopOtherAcc]
except:
	TopOtherName = "?"

CovScore = CovPanCovPct
OtherScore = OtherPanCovPct

s = "cov_score=%.1f;" % CovScore
s += "mapped=%d;" % CovMapped
s += "pctid=%.1f;" % CovPanPctId
s += "pan=%s;" % CovPanCartoon
s += "topacc=%s;" % TopCovAcc
s += "topname=%s;" % TopCovName
print(s, file=fSum)

s = "other_score=%.1f;" % OtherScore
s += "mapped=%d;" % OtherMapped
s += "pctid=%.1f;" % OtherPanPctId
s += "pan=%s;" % OtherPanCartoon
s += "topacc=%s;" % TopOtherAcc
s += "topname=%s;" % TopOtherName
print(s, file=fSum)

def GetOutputLine(Acc):
	Hits = AccToHits[Acc]
	try:
		Name = AccToName[Acc]
	except:
		Name = "?"
	try:
		Len = AccToLen[Acc]
	except:
		Len = None

	SumBases = AccToSumBases[Acc]
	SumBasesPctId = AccToSumBasesPctId[Acc]
	PctId = 0
	Depth = 0
	if SumBases > 0:
		PctId = float(SumBasesPctId)/SumBases
		if Len != None:
			Depth = float(SumBases)/Len

	CovFract = GetCovFract(AccToCoverageVec[Acc])
	Cartoon = MakeCartoon(AccToCoverageVec[Acc])

	s = "acc=" + Acc + ";"
	s += "pctid=%.1f;" % PctId
	if Len != None:
		s += "cov=%.4f;" % CovFract
		s += "len=" + str(Len) + ";"
		s += "hits=%d;" % Hits
		s += "depth=%.3g;" % Depth
		s += "coverage=" + Cartoon + ";"
	s += "name=%s;" % Name
	return s

for i in Order:
	Acc = Accs[i]
	if IsPan(Acc):
		continue
	s = GetOutputLine(Acc)
	print(s, file=fSum)

for i in range(0, CVG_BINS):
	Seq = CovSeqs[i]
	if Seq != None:
		print(">Cov.%d" % (i+1), file=fSum)
		print(Seq, file=fSum)

for i in range(0, CVG_BINS):
	Seq = OtherSeqs[i]
	if Seq != None:
		print(">Other.%d" % (i+1), file=fSum)
		print(Seq, file=fSum)

fSum.close()
