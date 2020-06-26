
"""
PART I   Basic Functions
"""


import re
import io


def CountMismatch(ref, alt):
    count = 0
    for i in range(len(ref)):
        if ref[i] != alt[i]:
            count += 1
    return count


def TypeofMutation(ref, alt):
    mark = None
    NormalBase = r"^([ATGCMN])+$"
    SymbolicMa = r"(^\\.)|(\\.$)|(^<)|(>$)|(\\[)|(\\])"
    if re.fullmatch(NormalBase, ref) == None:
        mark = "Unknown"
    elif re.fullmatch(NormalBase, alt) != None:
        if len(ref) ==  len(alt):
            if len(ref) == 1:
                mark = "SNP"
            elif CountMismatch(ref, alt) == 1:
                mark = "SNP"
            else:
                mark = "MNP"
        else:
            if len(ref) < len(alt):
                if len(ref) == 1:
                    if ref[0] == alt[0]:
                        mark = "Insertion"
                    else:
                        mark = "Complex"
                elif (ref[0] == alt[0]) & (alt.endswith(ref[1])):
                    mark = "Insertion"
                else:
                    mark = "Complex"
            else:
                if len(alt) == 1:
                    if alt[0] == ref[0]:
                        mark = "Deletion"
                    else:
                        mark = "Complex"
                elif (ref[0] == alt[0]) & (ref.endswith(alt[1])):
                    mark = "Deletion"
                else:
                    mark = "Complex"
    elif alt == "*":
        mark = "Star"
    elif re.fullmatch(SymbolicMa, alt):
        mark = "Symbolic"
    else:
        mark = "Unknown"

    return mark


def isTransition(ref, alt):
    TransitionPairs = [['A','G'], ['G', 'A'], ['C', 'T'], ['T', 'C']]
    if TypeofMutation(ref, alt) == "SNP":
        if len(ref) == 1:
            return([ref, alt] in TransitionPairs)
        else:
            SNPIndex = [i for i in range(len(ref)) if ref[i] != alt[i]]
            index = SNPIndex[0]
            return([ref[index], alt[index]] in TransitionPairs)
    else:
        return False


def isTransversion(ref, alt):
    if TypeofMutation(ref, alt) == "SNP":
        return(not isTransition(ref, alt))
    else:
        return False


def isPhased(c):
    SymbolPhase = r"^([\.0-9]/).*$"
    SymbolUnphase = r"^([\.0-9]\|).*$"
    if re.fullmatch(SymbolPhase, c) != None:
        return("Phased")
    elif re.fullmatch(SymbolUnphase, c) != None:
        return("Unphased")
    else:
        return("Unknown")


def diploidGtIndex(j, k):
    return( k * (k + 1) / 2 + j )


def diploidGtIndexWithSwap(j, k):
    if (j < k):
        return( diploidGtIndex(j, k) )
    else:
        return( diploidGtIndex(k, j) )


def isCalled(c):
    SymbolCalled = r"^([0-9]).*$"
    SymbolNotCalled = r"^[\.].*$"
    if re.fullmatch(SymbolCalled, c) != None:
        return("Called")
    elif re.fullmatch(SymbolNotCalled, c) != None:
        return("NotCalled")
    else:
        return("Unknown")


def AlleleLen(c):
    SymbolDiploid = r"(^([\.0-9]/).*$)|(^([\.0-9]\|).*$)"
    SymbolHaploid = r"[\.0-9]$"
    if re.fullmatch(SymbolDiploid, c) != None:
        return(2)
    elif re.fullmatch(SymbolHaploid, c) != None:
        return(1)
    else:
        return(0)
    # Return a number of 0, 1, or 2


def alleleRepr(c):
    if isPhased(c) == "Phased":
        j = int(c.split("/")[0])
        k = int(c.split("/")[1])
        return( diploidGtIndex(j,j + k) )
    elif isPhased(c) == "Unphased":
        j = int(c.split("|")[0])
        k = int(c.split("|")[1])
        return( diploidGtIndexWithSwap(j, k) )


smallAllelePair = [[0,0],
                   [0,1], [1,1],
                   [0,2], [1,2], [2,2],
                   [0,3], [1,3], [2,3], [3,3],
                   [0,4], [1,4], [2,4], [3,4], [4,4],
                   [0,5], [1,5], [2,5], [3,5], [4,5], [5,5],
                   [0,6], [1,6], [2,6], [3,6], [4,6], [5,6], [6,6],
                   [0,7], [1,7], [2,7], [3,7], [4,7], [5,7], [6,7], [7,7]]


def allelePairSqrt(i):
    k = int(((8 * i + 1) ** 0.5) / 2 - 0.5)
    j = i - k * (k + 1) / 2
    return([min(k,j), max(k,j)])


def detectAllelePair(i):
    i = int(i)
    if ( i < len(smallAllelePair) ):
        return( smallAllelePair[i] )
    else:
        return( allelePairSqrt(i) )


def allelePair(c):
    if isPhased(c) == "Unphased":
        return( detectAllelePair( alleleRepr(c) ) )
    elif isPhased(c) == "Phased":
        g = detectAllelePair( alleleRepr(c) )
        new_j = int(min(g[0], (g[1] - g[0])))
        new_k = int(max(g[0], (g[1] - g[0])))
        return( [new_j, new_k] )

void = r"(^[\.].*$)"

def isHomRef(c):
    allelelen = AlleleLen(c)
    if re.fullmatch( void, c ) != None:
        return False
    elif allelelen == 0:
        return False
    elif (allelelen == 1) or (allelelen == 2):
        return (alleleRepr(c) == 0)
    else:
        return False


def isHet(c):
    allelelen = AlleleLen(c)
    if re.fullmatch( void, c ) != None:
        return False
    elif (allelelen == 0) or (allelelen == 1):
        return False
    elif allelelen == 2:
        p = allelePair(c)
        return( (alleleRepr(c) > 0) & (p[0] != p[1]) )
    else:
        return False


def isHomVar(c):
    allelelen = AlleleLen(c)
    if re.fullmatch( void, c ) != None:
        return False
    if allelelen == 0:
        return False
    elif allelelen == 1:
        return( alleleRepr(c) > 0 )
    elif allelelen == 2:
        p = allelePair(c)
        return( (alleleRepr(c) > 0) & (p[0] == p[1]) )
    else:
        return False


def SetMissing(GT, AD, DP, GQ, PL, AF):
    if DP == "." or GQ == "." or PL[0] == ".":
        GT = "./."
    elif isCalled(GT) == "Called":
        phaseinfo = isPhased(GT)
        if phaseinfo == "Phased":
            gt = GT.split("/")
        elif phaseinfo == "Unphased":
            gt = GT.split("|")
        else:
            gt = [0,0]
        ref = int(gt[0])
        alt = int(gt[1])
        if ((int(DP)) < 10) or ((int(DP)) > 400):
            GT = "./."
        else:
            if isHomRef(GT) == True:
                if ((int(AD[ref])/int(DP)) < 0.9) or (int(GQ) < 20):
                    GT = "./."
            elif isHomVar(GT) == True:
                if ((int(AD[alt])/int(DP)) < 0.9) or (int(PL[0]) < 20) or ((AF > 0.5) and (int(GQ) < 20)):
                    GT = "./."
            elif isHet(GT) == True:
                if (((int(AD[alt])+int(AD[ref]))/int(DP)) < 0.9) or ((int(AD[alt])/int(DP)) < 0.2) or (int(PL[0]) < 20):
                    GT = "./."
    return(GT)



"""
PART II   Sample QC
"""

"""
To begin with, build an index for every sample
"""


samples = {}
OrderSample = []
for StrIn in open("../work/VariantPass_1.txt", "r"):
    if StrIn.startswith('#CHROM'):
        TempStr = StrIn.split('\n')[0]
        OrderSample = TempStr.split('\t')
        for redunt in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', '']:
            if redunt in OrderSample:
                OrderSample.remove(redunt)
        for patient in OrderSample:
            samples[patient] = {"nTi"     :   0,
                                "nTv"     :   0,
                                "nHet"    :   0,
                                "nHomRef" :   0,
                                "nHomVar" :   0,
                                "nCalled" :   0,
                                "nNotCalled": 0,
                                "nInsertion": 0,
                                "nDeletion":  0,
                                "nStar":      0}
    else:
        TempStr = StrIn.split('\n')[0]
        TempStr = TempStr.split('\t')
        Ref = TempStr[3]
        Alt = TempStr[4].split(",")
        AF_INFO = re.search(r'(?<=AF=)\d*', TempStr[7])
        AF = float(AF_INFO.group(0))
        for i in range(len(OrderSample)):
            ID = OrderSample[i]
            #samples[OrderSample[i]].append(TempStr[i+9])
            Cell = TempStr[i+9].split(":")
            GT = str(Cell[0])
            AD = Cell[1].split(",")
            DP = Cell[2]
            GQ = Cell[3]
            PL = Cell[7].split(",")

            GT = SetMissing(GT, AD, DP, GQ, PL, AF)

            if isCalled(GT) == "Called":
                samples[ID]["nCalled"] += 1
                # Calculate Ti and Tv
                phaseinfo = isPhased(GT)
                if phaseinfo == "Phased":
                    gt = GT.split("/")
                elif phaseinfo == "Unphased":
                    gt = GT.split("|")
                else:
                    gt = [0,0]
                for item in gt:
                    alleleIndex = int(item)
                    if alleleIndex != 0:
                        alt = Alt[alleleIndex-1]
                        typeofmutation = TypeofMutation(Ref, alt)
                        if typeofmutation == "Insertion":
                            samples[ID]["nInsertion"] += 1
                        elif typeofmutation == "Deletion":
                            samples[ID]["nDeletion"] += 1
                        elif typeofmutation == "Star":
                            samples[ID]["nStar"] += 1
                        elif typeofmutation == "SNP":
                            if isTransition(Ref, alt) == True:
                                samples[ID]["nTi"] += 1
                            elif isTransversion(Ref, alt) == True:
                                samples[ID]["nTv"] += 1
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass

                # Calculate nHet and nHom
                if isHet(GT) == True:
                    samples[ID]["nHet"] += 1
                elif isHomRef(GT) == True:
                    samples[ID]["nHomRef"] += 1
                elif isHomVar(GT) == True:
                    samples[ID]["nHomVar"] += 1
                else:
                    pass

            else:
                samples[ID]["nNotCalled"] += 1

# Write out temp calculating results
Header = ['ID', 'nTi', 'nTv',
          'nHet', 'nHomRef', 'nHomVar',
          'nCalled', 'nNotCalled', 'nInsertion',
          'nDeletion', 'nStar', '\n']

with open("../work/SampleQC_1.txt", "w") as FileOut:
    FileOut.write('\t'.join(Header))
    for key in samples.keys():
        StrOut = [key,
                  str(samples[key]["nTi"]),
                  str(samples[key]["nTv"]),
                  str(samples[key]["nHet"]),
                  str(samples[key]["nHomRef"]),
                  str(samples[key]["nHomVar"]),
                  str(samples[key]["nCalled"]),
                  str(samples[key]["nNotCalled"]),
                  str(samples[key]["nInsertion"]),
                  str(samples[key]["nDeletion"]),
                  str(samples[key]["nStar"]),
                  "\n"]
        FileOut.write('\t'.join(StrOut))
FileOut.close()
