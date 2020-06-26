"""
Objective: Sparse original data into small files with 30 samples each
"""

SmallFileNumber = 5
FileList = [open("../work/VariantPass_%i.txt" %i, "w") for i in range(1,SmallFileNumber+1)]

originalCount = 0
passCount = 0


for StrIn in open("../data/originaldata.vcf", "r"):
    if StrIn.startswith('#') or StrIn.startswith('"#') or StrIn.startswith('\n'):
        if StrIn.startswith('#CHROM'):
            TempStr = StrIn.split('\n')[0]
            TempStr = TempStr.split('\t')


            for j in range( SmallFileNumber ):
                if j != (SmallFileNumber-1):
                    StrOut = TempStr[0:9]
                    FileOut = FileList[j]
                    start = 9 + j*30
                    end = 39 + j*30
                    for i in range(start, end):
                        StrOut.append(TempStr[i])
                    StrOut.append('\n')
                    FileOut.write('\t'.join(StrOut))
                else:
                    StrOut = TempStr[0:9]
                    FileOut = FileList[j]
                    start = 9 + j*30
                    end = len(TempStr)
                    for i in range(start, end):
                        StrOut.append(TempStr[i])
                    StrOut.append('\n')
                    FileOut.write('\t'.join(StrOut))

        else:
            pass
    else:
        TempStr = StrIn.split('\n')[0]
        TempStr = TempStr.split('\t')
        mark = TempStr[6]
        if mark == 'PASS':

            for j in range( SmallFileNumber ):
                if j != (SmallFileNumber-1):
                    StrOut = TempStr[0:9]
                    FileOut = FileList[j]
                    start = 9 + j*30
                    end = 39 + j*30
                    for i in range(start, end):
                        StrOut.append(TempStr[i])
                    StrOut.append('\n')
                    FileOut.write('\t'.join(StrOut))
                else:
                    StrOut = TempStr[0:9]
                    FileOut = FileList[j]
                    start = 9 + j*30
                    end = len(TempStr)
                    for i in range(start, end):
                        StrOut.append(TempStr[i])
                    StrOut.append('\n')
                    FileOut.write('\t'.join(StrOut))

            passCount += 1
            originalCount += 1

        else:
            originalCount += 1

for file in FileList:
    file.close()


with open("../work/Number_VariantPass.txt", "w") as FileOut:
    TempStr = ['Number of Variant PASS', str(passCount)+'\n'+'Number of Original Variants', str(originalCount)+'\n']
    StrOut = '\t'.join(TempStr)
    FileOut.write(StrOut)
FileOut.close()
