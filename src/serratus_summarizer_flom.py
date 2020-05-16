#!/usr/bin/python3

# Author: Robert C. Edgar
# email robert@drive5.com

import sys
import os

# **************************************************
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# MUST set RAISEX to False for production
# If True, exceptions are raised, which is helpful
# for debugging but fatal for generating a BAM file.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# **************************************************
RAISEX = False

COV_BINS = 32
PAN_GENOME = "_mega_"
MIN_LENGTH = 1000

InputFileName = sys.argv[1]
SummaryFileName = sys.argv[2]
OutputFileName = sys.argv[3]
TinyhitFileName = None

if len(sys.argv) > 4:
	TinyhitFileName = sys.argv[4]
	try:
		fTinyhit = open(TinyhitFileName, "w")
	except:
		fTinyhit = None

Dir = os.getenv("SUMZER_DIR", ".")
if not Dir.endswith('/'):
	Dir += "/"

AccLenTaxFileName = Dir + "acc_len_taxid.txt"
TaxDescFileName = Dir + "taxid_desc.txt"

fIn = open(InputFileName)
fSum = open(SummaryFileName, "w")
fOut = None
if OutputFileName != "-":
	fOut = open(OutputFileName, "w")

AccToTax = {}
AccToLen = {}
TaxToDesc = {}
AccToSumBases = {}
AccToSumBasesPctId = {}
AccToCoverageVec = {}

AccToCoverageVec[PAN_GENOME] = [0]*COV_BINS
AccToLen[PAN_GENOME] = 30000

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

for Line in open(AccLenTaxFileName):
	Fields = Line[:-1].split('\t')
	Acc = Fields[0]
	Len = int(Fields[1])
	Tax = Fields[2]
	AccToTax[Acc] = Tax
	AccToLen[Acc] = Len

for Line in open(TaxDescFileName):
	Fields = Line[:-1].split('\t')
	if len(Fields) == 2:
		Tax = Fields[0]
		Desc = Fields[1]
		TaxToDesc[Tax] = Desc

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
		AccToCoverageVec[Acc] = [0]*COV_BINS
		AccToCoverageVec[Acc][TBin] += 1

Mapped = 0
Unmapped = 0
SumL = 0
MaxL = 0
ShortHits = 0
for Line in fIn:
	# Echo stdin to stdout, like tee
	if fOut != None:
		fOut.write(Line)

	# Ignore SAM headers
	if Line.startswith('@'):
		continue

	# Wrap everything in a try-except block because
	# the summarizer MUST not crash the pipeline!
	# Better no summary than no BAM file.
	try:
		Fields = Line[:-1].split('\t')

		# ReadLabel = Fields[0]
		Flags = int(Fields[1])
		if (Flags & 0x4) == 0x4:
			Unmapped += 1
			continue

		Mapped += 1
		Acc = Fields[2]
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

		if fTinyhit != None:
			print(Acc + "\t" + str(TPos) + "\t" + str(L) + "\t" + str(NM) + "\t" + EE, file=fTinyhit)

		PctId = 0
		if NM != None and L > 0:
			PctId = (L - NM)*100.0/L

		TL = None
		try:
			TL = AccToLen[Acc]
			TBin = ((TPos-1)*COV_BINS)//TL
			if TBin < 0:
				TBin = 0
			if TBin >= COV_BINS:
				TBin = COV_BINS - 1
		except:
			TBin = 0

		AddHit(Acc, TBin, L, PctId)
		if TL != None:
			if TL >= MIN_LENGTH:
				AddHit(PAN_GENOME, TBin, L, PctId)
			else:
				ShortHits += 1

	except:
		if RAISEX:
			raise
		pass

if Mapped > 0:
	MeanL = SumL//Mapped
else:
	MeanL = 0

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

TopAcc = "?"
for k in Order:
	Acc = Accs[k]
	if Acc != PAN_GENOME:
		TopAcc = Acc
		break

PanCartoon = MakeCartoon(AccToCoverageVec[PAN_GENOME])
PanCovFract = GetCovFract(AccToCoverageVec[PAN_GENOME])

try:
	PanSumBases = AccToSumBases[PAN_GENOME]
except:
	PanSumBases = 0

try:
	PanSumBasesPctId = AccToSumBasesPctId[PAN_GENOME]
except:
	PanSumBasesPctId = 0

PanCovPct = PanCovFract*100.0
PanPctId = 0.0
if PanSumBasesPctId > 0:
	PanPctId = PanSumBasesPctId/PanSumBases

try:
	TopTax = AccToTax[TopAcc]
except:
	TopTax = "?"
try:
	TopDesc = TaxToDesc[TopTax]
except:
	TopDesc = "?"

Score = PanCovPct

s = "score=%.1f;" % Score
s += "pctid=%.1f;" % PanPctId
s += "pancov=%s;" % PanCartoon
s += "topacc=%s;" % TopAcc
s += "topdesc=%s;" % TopDesc

print(s, file=fSum)

print("pancovpct=%.3g;" % PanCovPct, file=fSum)
print("unmapped=%d;" % Unmapped, file=fSum)
print("mapped=%d;" % Mapped, file=fSum)
print("short_hits=%d;" % ShortHits, file=fSum)
print("mean_aln_length=%d;" % MeanL, file=fSum)
print("max_aln_length=%d;" % MaxL, file=fSum)

def GetOutputLine(Acc):
	Hits = AccToHits[Acc]
	try:
		Tax = AccToTax[Acc]
	except:
		Tax = "?"
	try:
		Desc = TaxToDesc[Tax]
	except:
		Desc = "?"
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
		s += "tax=%s;" % Tax
		s += "coverage=" + Cartoon + ";"
	s += "desc=%s;" % Desc
	return s

for i in Order:
	Acc = Accs[i]
	if Acc == PAN_GENOME:
		continue
	s = GetOutputLine(Acc)
	print(s, file=fSum)

fSum.close()
