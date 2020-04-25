#!/usr/bin/python2

# Author Robert C. Edgar
# email robert@drive5.com

from __future__ import print_function
import sys
import os

COV_BINS = 16

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

d = {}
Order = []
Keys = []
def Cmp__(i, j):
	global d, Keys

	ki = Keys[i]
	kj = Keys[j]

	ni = d[ki]
	nj = d[kj]

	if ni < nj:
		return 1
	elif ni > nj:
		return -1
	return 0

def GetOrder(Dict):
	global d, Order, Keys

	d = Dict
	Keys = d.keys()
	N = len(Keys)
	Order = range(0, N)
	Order.sort(Cmp__)
	return Order

def MakeCartoon(v):
	Max = max(v)
	s = ""
	for i in v:
		k = (i*4)/Max
		s += "_.oO@"[k]
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
			TBin = ((TPos-1)*COV_BINS)/TL
			if TBin < 0:
				TBin = 0
			if TBin >= COV_BINS:
				TBin = COV_BINS - 1
		except:
			TBin = 0

		try:
			AccToHits[Acc] += 1
			AccToSumBases[Acc] += L
			AccToSumBasesPctId[Acc] += L*PctId
			AccToCoverageVec[Acc][TBin] += 1
		except:
			AccToHits[Acc] = 0
			AccToSumBases[Acc] = L
			AccToSumBasesPctId[Acc] = L*PctId
			AccToCoverageVec[Acc] = [0]*COV_BINS
			AccToCoverageVec[Acc][TBin] += 1
	except:
		pass

MappedReversePct = 0.0
if Mapped > 0:
	MeanL = SumL/Mapped
	MappedReversePct = (MappedReverse*100.0)/Mapped
else:
	MeanL = 0

Accs = AccToHits.keys()
Order = GetOrder(AccToHits)

print("unmapped=%d" % Unmapped, file=fSum)
print("mapped=%d" % Mapped, file=fSum)
print("mapped_reverse=%d" % MappedReverse, file=fSum)
print("mapped_reverse_pct=%.2f" % MappedReversePct, file=fSum)
print("mean_aln_length=%d" % MeanL, file=fSum)
print("max_aln_length=%d" % MaxL, file=fSum)

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

	SumBases = AccToSumBases[Acc]
	SumBasesPctId = AccToSumBasesPctId[Acc]
	PctId = 0
	Depth = 0
	if SumBases > 0:
		PctId = float(SumBasesPctId)/SumBases
		if Len != None:
			Depth = float(SumBases)/Len

	Cartoon = MakeCartoon(AccToCoverageVec[Acc])

	s = "acc=" + Acc + ";"
	s += "hits=%d;" % Hits
	if Len != None:
		s += "len=" + str(Len) + ";"
		s += "depth=%.3g;" % Depth
		s += "pctid=%.1f;" % PctId
		s += "tax=%s;" % Tax
	s += "coverage=" + Cartoon + ";"
	s += "desc=%s;" % Desc
	print(s, file=fSum)
fSum.close()
