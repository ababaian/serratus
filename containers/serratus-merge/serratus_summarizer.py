#!/usr/bin/python3

# Author: Robert C. Edgar
# email robert@drive5.com

import sys
import os

COV_BINS = 32
MIN_COMPLETE_LEN = 25000
PAN_GENOME = "pan_genome"

InputFileName = sys.argv[1]
SummaryFileName = sys.argv[2]
OutputFileName = sys.argv[3]

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
MappedReverse = 0
Unmapped = 0
SumL = 0
MaxL = 0

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
		# CIGAR = Fields[5]
		# RNEXT = Fields[6]
		# PNEXT = Fields[7]
		# TL = int(Fields[8])
		SEQ = Fields[9]
		# QUAL = Fields[10]
		L = len(SEQ)
		SumL += L
		if L > MaxL:
			MaxL = L

		if Acc.lower().find("reverse") >= 0:
			MappedReverse += 1

		AS = None
		NM = None
		for Field in Fields[11:]:
			if Field.startswith("NM:i:"):
				s = Field.replace("NM:i:", "")
				try:
					NM = int(s)
				except:
					NM = None
				break

		PctId = 0
		if NM != None and L > 0:
			PctId = (L - NM)*100.0/L

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
		if TL >= MIN_COMPLETE_LEN:
			AddHit(PAN_GENOME, TBin, L, PctId)

	except:
		pass

MappedReversePct = 0.0
if Mapped > 0:
	MeanL = SumL//Mapped
	MappedReversePct = (MappedReverse*100.0)/Mapped
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
Order = GetOrder(AccToHits)

try:
	PanCov = AccToCoverageVec[PAN_GENOME]
except:
	PanCov = COV_BINS*[0]

PanCovFract = GetCovFract(PanCov)
MaxCovFract1k = 0
for Acc in Accs:
	try:
		Len = AccToLen[Acc]
	except:
		continue
	try:
		CovFract = GetCovFract(AccToCoverageVec[Acc])
	except:
		continue

	if CovFract > MaxCovFract1k:
		MaxCovFract1k = CovFract

PanCovPct = PanCovFract*100.0
MaxCovPct1k = MaxCovFract1k*100.0

Score = 0
if Mapped > 50000:
	Score += 50
elif Mapped > 25000:
	Score += 25
elif Mapped > 10000:
	Score += 15
elif Mapped > 1000:
	Score += 10
elif Mapped > 100:
	Score += 5
elif Mapped > 0:
	Score += 1

Score += (PanCovPct + MaxCovPct1k)/4.0

print("score=%.0f;" % Score, file=fSum)
print("pancovpct=%.3g;" % PanCovPct, file=fSum)
print("maxcov1kpct=%.3g;" % MaxCovPct1k, file=fSum)
print("unmapped=%d;" % Unmapped, file=fSum)
print("mapped=%d;" % Mapped, file=fSum)
print("mapped_reverse=%d;" % MappedReverse, file=fSum)
print("mapped_reverse_pct=%.2f;" % MappedReversePct, file=fSum)
print("mean_aln_length=%d;" % MeanL, file=fSum)
print("max_aln_length=%d;" % MaxL, file=fSum)

for i in Order:
	Acc = Accs[i]
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
	try:
		IsComplete = AccToComplete[Acc]
	except:
		IsComplete = False

	if Acc == PAN_GENOME:
		Desc = "Pan-genome"

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
	s += "hits=%d;" % Hits
	if Len != None:
		s += "len=" + str(Len) + ";"
		s += "depth=%.3g;" % Depth
		s += "pctid=%.1f;" % PctId
		s += "tax=%s;" % Tax
		s += "cov=%.4f;" % CovFract
		s += "coverage=" + Cartoon + ";"
	s += "desc=%s;" % Desc
	print(s, file=fSum)
fSum.close()
