"""
Topic: Annotation and QC for CCDG Freeze 1 Data
 Date: 06/27/2018
   By: Yushi F.T.
------------------------------------------------------
General descriptions for objectives of this source code:

* Annotate and perform sample QC on CCDG Freeze 1 Data based on the following
    indicators:
-   1). Ratio of transition and transversion (r_Ti_Tv)
-   2). Ratio of heterozygous calls and homozygous alternate calls (r_Het_HomVar)

* After filtering patients with r_Ti_Tv and r_Het_HomVar outliers, annotate and
    perform variant QC on CCDG Freeze 1 Data based on the following indicators:
-   1). Coverage
-   2). Call rate (fraction of calls non-missing)
-   3). Ratio of heterozygous calls and homozygous alternate calls (r_Het_HomVar)
-   4). Exact P-value of HWE (HWE_P)

* Function of calculating exact p-value of HWE is written with the following
    resources as references:
-   1). Wigginton, J. E., Cutler, D. J., & Abecasis, G. R. (2005). A note on
          exact tests of hardy-weinberg equilibrium. American Journal of Human
          Genetics, 76(5), 887-93.
-   2). Center for Statistical Genetics, University of Michigan, SNP-HWE(R)
          sources. (http://csg.sph.umich.edu/abecasis/Exact/)

* Other algorithms were mainly referred to the software Hail developed by
    Robert's group and Neale lab at the Stanley Center for Psychiatric Research
    of the Broad Institute of MIT and Harvard and the Analytic and Translational
    Genetics Unit of Massachusetts General Hospital. (https://hail.is/index.html)

* This source code includes three parts as:
-   PART I  : Defination for Basic Functions
-   PART II : Sample Annotation and QC
-   PART III: Variant Annotation and QC
------------------------------------------------------
"""

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
        return([ref, alt] in TransitionPairs)
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
        if (int(DP) < 10) or (int(DP) > 400):
            GT = "./."
        else:
            if isHomRef(GT) == True:
                if ((int(AD[ref])/(int(DP)) < 0.9) or (int(GQ) < 20)):
                    GT = "./."
            elif isHomVar(GT) == True:
                if ((int(AD[alt])/int(DP)) < 0.9) or (int(PL[0]) < 20) or ((AF > 0.5) and (int(GQ) < 20)):
                    GT = "./."
            elif isHet(GT) == True:
                if (((int(AD[alt])+int(AD[ref]))/int(DP)) < 0.9) or ((int(AD[alt])/int(DP)) < 0.2) or (int(PL[0]) < 20):
                    GT = "./."
    return(GT)

"""
PART II   Variant QC
"""

### Calculate and write out row by row
originalCount = 0
passCount = 0
with open("../work/VariantQC.txt", "w") as FileOut:
    header = ['CHROM', 'POS', 'ID',
              'REF', 'ALT', 'QUAL',
              'Coverage', 'nHet', 'nHomRef', 'nHomVar', 'nCalled', 'nNotCalled','\n']
    FileOut.write( '\t'.join(header) )
    for StrIn in open("../data/originaldata.vcf", "r"):
        if StrIn.startswith('chr'):
            TempStr = StrIn.split('\n')[0]
            TempStr = TempStr.split('\t')
            mark = TempStr[6]
            if mark == "PASS":
                passCount += 1
                originalCount += 1
                nHet = 0
                nHomRef = 0
                nHomVar = 0
                nCalled = 0
                nNotCalled = 0
                Ref = TempStr[3]
                Alt = TempStr[4].split(",")
                AN_INFO = re.search(r'(?<=AN=)\d*', TempStr[7])
                DP_INFO = re.search(r'(?<=DP=)\d*', TempStr[7])
                AF_INFO = re.search(r'(?<=AF=)\d*', TempStr[7])
                AN_info = int(AN_INFO.group(0))
                DP_info = int(DP_INFO.group(0))
                AF = float(AF_INFO.group(0))

                for i in range((len(TempStr) - 9)):
                    Cell = TempStr[i+9].split(":")
                    GT = str(Cell[0])
                    AD = Cell[1].split(",")
                    DP = Cell[2]
                    GQ = Cell[3]
                    PL = Cell[7].split(",")

                    GT = SetMissing(GT, AD, DP, GQ, PL, AF)

                    if isCalled(GT) == "Called":
                        if isHomRef(GT) == True:
                            nHomRef += 1
                            nCalled += 1
                        elif isHet(GT) == True:
                            nHet += 1
                            nCalled += 1
                        elif isHomVar(GT) == True:
                            nHomVar += 1
                            nCalled += 1
                        else:
                            pass
                    else:
                        nNotCalled += 1

                # Calculate coverage
                Coverage = DP_info / (AN_info / 2)

                record = TempStr[0:6]
                record.append( str(Coverage) )
                record.append( str(nHet) )
                record.append( str(nHomRef) )
                record.append( str(nHomVar) )
                record.append( str(nCalled) )
                record.append( str(nNotCalled) )
                record.append( "\n" )
                StrOut = "\t".join(record)
                FileOut.write( StrOut )

            else:
                originalCount += 1
        else:
            pass

FileOut.close()


with open("../work/Number_VariantPass.txt", "w") as FileOut:
    TempStr = ['Number of Variant PASS:', str(passCount)+'\n'+'Number of Original Variants:', str(originalCount)+'\n']
    StrOut = '\t'.join(TempStr)
    FileOut.write(StrOut)
FileOut.close()
