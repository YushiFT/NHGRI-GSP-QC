"""
Topic: Annotation and QC for CCDG Freeze 1 Data
 Date: 05/23/2018
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



"""
PART II   Variant QC
"""

### Calculate and write out row by row
with open("../work/VariantQC.txt", "w") as FileOut:

    header = ['CHROM', 'POS', 'ID',
              'REF', 'ALT', 'QUAL',
              'Coverage', 'nHet', 'nHomRef', 'nHomVar', 'nCalled', 'nNotCalled','\n']

    FileOut.write( '\t'.join(header) )

    for StrIn in open("../work/VariantPass.txt", "r"):
        if StrIn.startswith('#CHROM'):
            pass
        else:
            TempStr = StrIn.split('\n')[0]
            TempStr = TempStr.split('\t')
            nHet = 0
            nHomRef = 0
            nHomVar = 0
            nCalled = 0
            nNotCalled = 0
            Ref = TempStr[3]
            Alt = TempStr[4].split(",")
            AN_INFO = re.search(r'(?<=AN=)\d*', TempStr[7])
            DP_INFO = re.search(r'(?<=DP=)\d*', TempStr[7])
            AN_info = int(AN_INFO.group(0))
            DP_info = int(DP_INFO.group(0))

            for i in range((len(TempStr) - 9)):
                Cell = TempStr[i+9].split(":")
                GT = str(Cell[0])
                #AD = Cell[1].split(",")
                #DP = int(Cell[2])
                #GQ = Cell[3]
                #PL = Cell[7].split(",")

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

FileOut.close()



"""
PART III   Sample QC
"""

"""
To begin with, build an index for every sample
"""


samples = {}
OrderSample = []
for StrIn in open("../work/VariantPass.txt", "r"):
    if StrIn.startswith('#CHROM'):
        TempStr = StrIn.split('\n')[0]
        OrderSample = TempStr.split('\t')
        for redunt in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            OrderSample.remove(redunt)
        for patient in OrderSample:
            samples[patient] = {"nTi"     :   0,
                                "nTv"     :   0,
                                "nHet"    :   0,
                                "nHomRef" :   0,
                                "nHomVar" :   0,
                                "nCalled" :   0,
                                "nNotCalled": 0}
    else:
        pass


### Now we deal with each cell (i.e. each specific sample)
for StrIn in open("../work/VariantPass.txt", "r"):
    if StrIn.startswith('#CHROM'):
        pass
    else:
        TempStr = StrIn.split('\n')[0]
        TempStr = TempStr.split('\t')
        Ref = TempStr[3]
        Alt = TempStr[4].split(",")
        for i in range(len(OrderSample)):
            ID = OrderSample[i]
            #samples[OrderSample[i]].append(TempStr[i+9])
            Cell = TempStr[i+9].split(":")
            GT = str(Cell[0])
            #AD = Cell[1].split(",")
            #DP = int(Cell[2])
            #GQ = Cell[3]
            #PL = Cell[7].split(",")

            if isCalled(GT) == "Called":
                samples[ID]["nCalled"] += 1
                # Calculate Ti and Tv
                if isPhased(GT) == "Phased":
                    gt = GT.split("/")
                elif isPhased(GT) == "Unphased":
                    gt = GT.split("|")
                else:
                    gt = [0,0]
                for item in gt:
                    alleleIndex = int(item)
                    if alleleIndex != 0:
                        alt = Alt[alleleIndex-1]
                        if isTransition(Ref, alt) == True:
                            samples[ID]["nTi"] += 1
                        elif isTransversion(Ref, alt) == True:
                            samples[ID]["nTv"] += 1
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
          'nCalled', 'nNotCalled', '\n']

with open("../work/SampleQC.txt", "w") as FileOut:
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
                  "\n"]
        FileOut.write('\t'.join(StrOut))
FileOut.close()
