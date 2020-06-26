"""
Filter out the header lines of VariantQC files
"""

for i in range(2, 19):
    FileOutPath = "../temp/small_" + str(i) + ".txt"
    FileInPath = "../temp/VariantQC_" + str(i) + ".txt"
    with open( FileOutPath, "w" ) as FileOut:
        for StrIn in open( FileInPath, "r" ):
            TempStr = StrIn.split('\t')
            StrOut = TempStr[7:12]
            StrOut.append('\n')
            FileOut.write('\t'.join(StrOut))
    FileOut.close()
