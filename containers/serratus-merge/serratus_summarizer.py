#!/usr/bin/python3

import sys
import os

CVG_BINS = 25
MAXCLIP = 0.1		# Max clipping for alignment to be "global"

THROWX = (os.getenv("SUMZER_THROWX", "YES") != "NO")
SRA = os.getenv("SUMZER_SRA", None)
SUMZER_COMMENT = os.getenv("SUMZER_COMMENT", None)

SAMInputFileName = sys.argv[1]
MetaFileName = sys.argv[2]
SummaryFileName = sys.argv[3]
SAMOutputFileName = sys.argv[4]

fIn = open(SAMInputFileName)
fSum = open(SummaryFileName, "w")
fOut = open(SAMOutputFileName, "w")

AlnCount = 0
SumReadLength = 0
ExceptCount = 0

AccToName = {}
AccToLen = {}
AccToFam = {}
AccToOffset = {}
AccToPanLength = {}
AccToGlbs = {}
AccToAlns = {}
AccToSumBases = {}
AccToSumBasesPctId = {}
AccToCoverageVec = {}

FamToPanLength = {}
FamToAccs = {}

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
	AccToOffset[Acc] = Offset
	AccToPanLength[Acc] = PanLength
	AccToFam[Acc] = Fam
	AccToLen[Acc] = Len

	AccToLen[Fam] = PanLength
	if Fam not in Fams:
		Fams.append(Fam)
		FamToAccs[Fam] = []
	FamToAccs[Fam].append(Acc)

for Fam in Fams:
	AccToCoverageVec[Fam] = [0]*CVG_BINS

def ParseCigar(CIGAR):
	AlnLen = 0
	Clip = 0

	if CIGAR == "*":
		return 100, 0
	try:
		AlnLen = 0
		n = 0
		for c in CIGAR:
			if c.isdigit():
				n = n*10 + (ord(c) - ord('0'))
			elif c.isalpha() or c == '=':
				if c == 'S' or c == 'H':
					Clip += n
				else:
					AlnLen += n
				n = 0
	except:
		if THROWX:
			raise
		AlnLen = 100
		Clip = 0
	return AlnLen, Clip

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
		None

def IncDict(Dict, Key, n=1):
	try:
		Dict[Key] += n
	except:
		Dict[Key] = n

def GetDict(Dict, Key, Default = None):
	try:
		return Dict[Key]
	except:
		return Default

def IncCvgVec(Acc, TBin):
	try:
		AccToCoverageVec[Acc][TBin] += 1
	except:
		AccToCoverageVec[Acc] = [0]*CVG_BINS
		AccToCoverageVec[Acc][TBin] = 1

def CountToSymbol(i):
	#		    012345678901  2^11=2048
	Symbols = "_.:uwaomUWAOM^"
	n = 1
	for c in Symbols:
		if i < n:
			return c
		n *= 2
	return "^"

def MakeCartoon(Acc):
	w = GetCvgVec(Acc)
	if w == None:
		return "~"

	s = ""
	for i in w:
		s += CountToSymbol(i)
	return s 

def GetBin(Pos, TL):
	Bin = ((Pos-1)*CVG_BINS)//TL
	if Bin < 0:
		Bin = 0
	if Bin >= CVG_BINS:
		Bin = CVG_BINS - 1
	return Bin

def AddHit(Acc, TBin, AlnLen, PctId, IsGlobal):
	if Acc == None:
		return

	IncCvgVec(Acc, TBin)
	IncDict(AccToAlns, Acc)
	IncDict(AccToSumBases, Acc, AlnLen)
	IncDict(AccToSumBasesPctId, Acc, AlnLen*PctId)
	if IsGlobal:
		IncDict(AccToGlbs, Acc)

def GetNM(TagFields):
	NM = None
	for Field in TagFields:
		if Field.startswith("NM:i:"):
			s = Field.replace("NM:i:", "")
			try:
				NM = int(s)
			except:
				if THROWX:
					raise
				NM = None
			break
	return NM

for Line in fIn:
	fOut.write(Line)

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
		TPos = int(Fields[3])
		CIGAR = Fields[5]
		SEQ = Fields[9]

		NM = GetNM(Fields[11:])
		AlnLen, Clip = ParseCigar(CIGAR)

		AlnCount += 1
		ReadLength = len(SEQ)
		SumReadLength += ReadLength

		IsGlobal = float(Clip)/float(AlnLen + Clip + 1) >= MAXCLIP

		PctId = 0
		if NM != None and AlnLen > 0:
			PctId = (AlnLen - NM)*100.0/AlnLen

		Desc = GetDict(AccToName, Acc)
		Fam = GetDict(AccToFam, Acc)
		PL = GetDict(AccToPanLength, Acc)
		TL = GetDict(AccToLen, Acc)
		Offset = GetDict(AccToOffset, Acc)

		IsComplete = False
		if Desc != None:
			IsComplete = (Desc.find("omplete genome") > 0)

		if TL != None:
			TBin = GetBin(TPos, TL)
			AddHit(Acc, TBin, AlnLen, PctId, IsGlobal)

		PBin = None
		if Fam != None:
			if IsComplete:
				PBin = GetBin(TPos, TL)
			elif Offset != None and PL != None:
				PBin = GetBin(TPos + Offset, PL)
			elif TL != None:
				PBin = GetBin(TPos, TL)
			if PBin != None:
				AddHit(Fam, PBin, PL, PctId, IsGlobal)

	except:
		if THROWX:
			raise
		pass

def GetPctId(Acc):
	SumBases = GetDict(AccToSumBases, Acc)
	SumBasesPctId = GetDict(AccToSumBasesPctId, Acc)
	PctId = 0
	if SumBases > 0:
		PctId = float(SumBasesPctId)/SumBases
	if PctId < 50:
		PctId = 50
	if PctId > 100:
		PctId = 100
	return PctId

def GetCvgBins(Acc, MinCount):
	v = GetCvgVec(Acc)

	N = len(v)
	if N == 0:
		return 0
	n = 0
	for x in v:
		if x >= MinCount:
			n += 1
	return n

def GetCvgPct(Acc, MinCount):
	n = GetCvgBins(Acc, MinCount)
	Pct = (n*100.0)/CVG_BINS
	return Pct

def GetTypicalBinCount(Acc):
	v = GetCvgVec(Acc)
	if v == None:
		return 0
	
	N = len(v)
	if N < 10:
		return 0

	w = v[:]
	w.sort()
	Sum = 0
	for i in range(2, N-2):
		Sum += w[i]
	Typ = float(Sum)/(N-4)
	return Typ

def GetIdentityWeight(PctId):
	FractId = PctId/100.0
	if FractId < 0.5:
		FractId = 0.5
	if FractId > 1.0:
		FractId = 1.0
	Weight = 1.0/(FractId**3)
	return Weight

def GetScore(Acc):
	Cvg1 = GetCvgBins(Acc, 1)
	Cvg8 = GetCvgBins(Acc, 8)
	RawScore = Cvg8*4 + (Cvg1 - Cvg8)
	PctId = GetPctId(Acc)
	Weight = GetIdentityWeight(PctId)
	Score = RawScore*Weight
	if Score > 100:
		Score = 100
	return Score

def GetDepth(Acc, PctId):
	global AvgReadLength

	TL = GetDict(AccToLen, Acc)
	if TL == None or TL == 0:
		return 0

	PctId = GetPctId(Acc)
	Typ = GetTypicalBinCount(Acc)
	EstimatedTotalBases = AvgReadLength*Typ*CVG_BINS
	Weight = GetIdentityWeight(PctId)
	RawDepth = EstimatedTotalBases/TL
	Depth = RawDepth*Weight
	return Depth

Accs = list(AccToAlns.keys())
AccToScore = {}
for Acc in Accs:
	Score = GetScore(Acc)
	AccToScore[Acc] = Score

AccOrder = GetOrder(AccToScore)
OrderAccs = list(AccToScore.keys())

AvgReadLength = 0
if AlnCount > 0:
	AvgReadLength = float(SumReadLength)/AlnCount

s = ""
if SRA != None:
	s += "sra=%s;" % SRA
s += "readlength=%.0f;" % AvgReadLength
if SUMZER_COMMENT != None:
	s += "SUMZER_COMMENT=" + SUMZER_COMMENT + ";"
print(s, file=fSum)

def GetFamLine(Fam):
	try:
		FamAccs = FamToAccs[Fam]
	except:
		return
	if len(FamAccs) == 0:
		return

	TopAcc = FamAccs[0]
	for k in AccOrder:
		Acc = OrderAccs[k]
		if Acc in FamAccs:
			TopAcc = Acc
			break

	TopScore = AccToScore[TopAcc]
	TopName = GetDict(AccToName, TopAcc, "_missing_name_")
	TopLen = GetDict(AccToLen, TopAcc)
	PL = GetDict(FamToPanLength, Fam)

	s = ""
	if PL != None:
		s += "panlen=%d;" % PL
	s += "top=" + TopAcc + ";"
	s += "topscore=%d;" % TopScore
	if TopLen != None:
		s += "toplen=%d;" % TopLen
	s += "topname=" + TopName + ";"
	return s

def GetOutputLine(Acc):
	Fam = GetDict(AccToFam, Acc);
	Score = GetScore(Acc)
	Alns = GetDict(AccToAlns, Acc)
	Glbs = GetDict(AccToGlbs, Acc, 0)
	PctId = GetPctId(Acc)
	PL = GetDict(AccToPanLength, Acc)
	Depth = GetDepth(Acc, PctId)
	Cartoon = MakeCartoon(Acc)
	Len = GetDict(AccToLen, Acc)

	IsFam = (Acc in Fams)

	s = ""
	if SRA != None:
		s += "sra=" + SRA + ";"

	if IsFam:
		s += "famcvg=" + Cartoon + ";"
	else:
		s += "seqcvg=" + Cartoon + ";"

	if IsFam:
		s += "fam=%s;" % Acc
	else:
		s += "seq=%s;" % Acc

	s += "score=%.0f;" % Score

	s += "pctid=%.0f;" % PctId
	s += "depth=%.1f;" % Depth
	s += "aln=%d;" % Alns
	s += "glb=%d;" % Glbs
	if Len == None:
		s += "len=?;"
	else:
		s += "len=%d;" % Len

	if IsFam:
		s += GetFamLine(Acc)
	else:
		if Fam != None:
			s += "family=%s;" % Fam
		Name = GetDict(AccToName, Acc, "_missing_name_")
		s += "name=%s;" % Name
	return s

for i in AccOrder:
	Acc = OrderAccs[i]
	if not Acc in Fams:
		continue
	s = GetOutputLine(Acc)
	if s != None:
		print(s, file=fSum)

for i in AccOrder:
	Acc = OrderAccs[i]
	if Acc in Fams:
		continue
	s = GetOutputLine(Acc)
	if s != None:
		print(s, file=fSum)

fSum.close()
