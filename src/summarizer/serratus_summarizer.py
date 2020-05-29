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

CVG_BINS = 25
COV_FAM = "Coronaviridae"

SAMInputFileName = sys.argv[1]
MetaFileName = sys.argv[2]
SummaryFileName = sys.argv[3]
SAMOutputFileName = sys.argv[4]

SUMZER_COMMENT = os.getenv("SUMZER_COMMENT", None)
if SUMZER_COMMENT != None:
	SUMZER_COMMENT = SUMZER_COMMENT.replace(";", "&")

fIn = open(SAMInputFileName)
fSum = open(SummaryFileName, "w")
fOut = open(SAMOutputFileName, "w")

AccToTax = {}
AccToName = {}
AccToLen = {}
AccToFam = {}
AccToOffset = {}
AccToPanLength = {}

FamToAccToCount = {}
FamToPL = {}

#   0     1    2      3		4		5
# Acc	Len	Name	Fam		Offset	Pan-genome length
Fams = []
for Line in open(MetaFileName):
	Fields = Line[:-1].split('\t')
	Acc = Fields[0]
	Len = int(Fields[1])
	Name = Fields[2]
	Fam = Fields[3]
	if len(Fields) > 4:
		Offset = int(Fields[4])
	else:
		Offset = None
	if len(Fields) > 5:
		PanLength = int(Fields[5])
	else:
		PanLength = None

	AccToName[Acc] = Name
	AccToLen[Acc] = Len
	AccToFam[Acc] = Fam
	AccToOffset[Acc] = Offset
	AccToPanLength[Acc] = PanLength
	FamToPL[Fam] = PanLength

	if Fam not in Fams:
		Fams.append(Fam)
		FamToAccToCount[Fam] = {}

AccToSumBases = {}
AccToSumBasesPctId = {}
AccToCoverageVec = {}
for Fam in Fams:
	AccToCoverageVec[Fam] = [0]*CVG_BINS

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

def GetCvgVec(Acc):
	try:
		return AccToCoverageVec[Acc]
	except:
		return [0]*CVG_BINS

def MakeCartoon(Acc):
	w = GetCvgVec(Acc)
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

def AllMatch(CIGAR):
	M = re.match("[0-9]+M", CIGAR)
	return M != None

AccToGlbs = {}
AccToAlns = {}

def AddHit(Acc, TBin, L, PctId, SoftClipped):
	if Acc == None:
		return

	try:
		AccToAlns[Acc] += 1
	except:
		AccToAlns[Acc] = 1

	try:
		AccToCoverageVec[Acc][TBin] += 1
	except:
		AccToCoverageVec[Acc] = [0]*CVG_BINS
		AccToCoverageVec[Acc][TBin] = 1

	if not SoftClipped:
		try:
			AccToGlbs[Acc] += 1
		except:
			AccToGlbs[Acc] = 1

	try:
		AccToSumBases[Acc] += L
		AccToSumBasesPctId[Acc] += L*PctId
	except:
		AccToSumBases[Acc] = L
		AccToSumBasesPctId[Acc] = L*PctId

def GetBin(Pos, L):
	Bin = ((Pos-1)*CVG_BINS)//L
	if Bin < 0:
		Bin = 0
	if Bin >= CVG_BINS:
		Bin = CVG_BINS - 1
	return Bin

CovMapped = 0
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
		Acc = Acc.split(':')[0]
		try:
			Fam = AccToFam[Acc]
		except:
			Fam = None

		#if Fam == COV_FAM:
		#	CovMapped += 1

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

		SoftClipped = (CIGAR.find("S") >= 0)

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

		PctId = 0
		if NM != None and L > 0:
			PctId = (L - NM)*100.0/L

		try:
			Offset = AccToOffset[Acc]
		except:
			Offset = None

		try:
			PL = AccToPanLength[Acc]
		except:
			PL = None

		try:
			TL = AccToLen[Acc]
		except:
			TL = None

		if TL != None:
			TBin = GetBin(TPos, TL)
			AddHit(Acc, TBin, TL, PctId, SoftClipped)

		try:
			Desc = AccToName[Acc]
			IsComplete = (Desc.find("omplete genome") > 0)
		except:
			IsComplete = False

		PBin = None
		if Fam != None:
			if IsComplete:
				PBin = GetBin(TPos, TL)
			elif Offset != None and PL != None:
				PBin = GetBin(TPos + Offset, PL)
			elif TL != None:
				PBin = GetBin(TPos, TL)
			if PBin != None:
				AddHit(Fam, PBin, PL, PctId, SoftClipped)

		if Fam != None:
			try:
				FamToAccToCount[Fam][Acc] += 1
			except:
				FamToAccToCount[Fam][Acc] = 1

		if PBin != None and AllMatch(CIGAR):
			SEQ = Fields[9]
			if Fam == COV_FAM:
				if CovSeqs[PBin] == None:
					CovSeqs[PBin] = SEQ
			else:
				if OtherSeqs[PBin] == None:
					OtherSeqs[PBin] = SEQ

	except:
		if RAISEX:
			raise
		pass

def GetCovFract(Acc):
	v = GetCvgVec(Acc)

	N = len(v)
	if N == 0:
		return 0
	n = 0
	for x in v:
		if x > 2:
			n += 1
	return float(n)/N

Accs = list(AccToAlns.keys())

AccToCovFract = {}
for Acc in Accs:
	CovFract = GetCovFract(Acc)
	AccToCovFract[Acc] = CovFract

AccOrder = GetOrder(AccToCovFract)

if SUMZER_COMMENT != None:
	s = "SUMZER_COMMENT=" + SUMZER_COMMENT + ";"
	print(s, file=fSum)

def GetOutputLineFam(Fam):
	try:
		Alns = AccToAlns[Fam]
	except:
		Alns = 0

	try:
		Glbs = AccToGlbs[Fam]
	except:
		Glbs = 0

	try:
		Alns = AccToAlns[Fam]
	except:
		Alns = 0

	try:
		PL = FamToPL[Fam]
	except:
		PL = None

	SumBases = AccToSumBases[Fam]
	SumBasesPctId = AccToSumBasesPctId[Fam]
	PctId = 0
	Depth = 0
	if SumBases > 0:
		PctId = float(SumBasesPctId)/SumBases

	CovFract = GetCovFract(Fam)
	Cartoon = MakeCartoon(Fam)

	AccToCount = FamToAccToCount[Fam]
	Order = GetOrder(AccToCount)
	Accs = list(AccToCount.keys())
	TopAcc = Accs[0]
	TopHits = AccToCount[TopAcc]
	try:
		TopName = AccToName[TopAcc]
	except:
		TopName = "?"
	try:
		TopLen = AccToLen[TopAcc]
	except:
		TopLen = None

	Score = 100.0*CovFract
	s = "family=" + Fam + ";"
	s += "score=%.0f;" % Score
	s += "pctid=%.0f;" % PctId
	s += "aln=%d;" % Alns
	s += "glb=%d;" % Glbs
	if PL != None:
		s += "panlen=%d;" % PL
	s += "cvg=" + Cartoon + ";"
	s += "top=" + TopAcc + ";"
	s += "topaln=" + str(TopHits) + ";"
	if TopLen == None:
		s += "toplen=?;"
	else:
		s += "toplen=%d;" % TopLen
	s += "topname=" + TopName + ";"
	return s

def GetOutputLine(Acc):
	if Acc in Fams:
		return GetOutputLineFam(Acc)

	try:
		Glbs = AccToGlbs[Acc]
	except:
		Glbs = 0

	try:
		Alns = AccToAlns[Acc]
	except:
		Alns = 0

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

	CovFract = GetCovFract(Acc)
	Cartoon = MakeCartoon(Acc)

	s = "acc=" + Acc + ";"
	s += "pctid=%.1f;" % PctId
	s += "aln=%d;" % Alns
	s += "glb=%d;" % Glbs
	if Len != None:
		s += "len=%d;" % Len
	try:
		Fam = AccToFam[Acc]
	except:
		Fam = None
	if Len != None:
		s += "cvgpct=%.0f;" % (CovFract*100)
		s += "len=" + str(Len) + ";"
		s += "depth=%.3g;" % Depth
		s += "cvg=" + Cartoon + ";"
	if Fam != None:
		s += "fam=" + Fam + ";"
	s += "name=%s;" % Name
	return s

for i in AccOrder:
	Acc = Accs[i]
	if not Acc in Fams:
		continue
	s = GetOutputLineFam(Acc)
	print(s, file=fSum)

for i in AccOrder:
	Acc = Accs[i]
	if Acc in Fams:
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
