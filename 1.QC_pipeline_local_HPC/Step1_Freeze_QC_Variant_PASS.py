"""
Step 1: Filter out variants without PASS mark
"""
originalCount = 0
passCount = 0
with open("../work/VariantPass.txt", "w") as FileOut:
    for StrIn in open("../data/originaldata.vcf", "r"):
        if StrIn.startswith('#') or StrIn.startswith('"#'):
            if StrIn.startswith('#CHROM'):
                FileOut.write(StrIn)
            else:
                pass
        else:
            TempStr = StrIn.split('\n')[0]
            TempStr = TempStr.split('\t')
            mark = TempStr[6]
            if mark == "PASS":
                FileOut.write(StrIn)
                passCount += 1
                originalCount += 1
            else:
                originalCount += 1
FileOut.close()


with open("../work/Number_VariantPass.txt", "w") as FileOut:
    TempStr = ['Number of Variant PASS', str(passCount), '\n',
               'Number of Original Variants', str(originalCount), '\n']
    StrOut = '\t'.join(TempStr)
    FileOut.write(StrOut)
FileOut.close()
