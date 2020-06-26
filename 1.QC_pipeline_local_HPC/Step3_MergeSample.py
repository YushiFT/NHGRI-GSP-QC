"""
Merge all samples together
"""


with open("../work/SampleQC.txt", "w") as FileOut:
    for StrIn in open("../work/SampleQC_1.txt", "r"):
        FileOut.write(StrIn)

    for i in range(2,151):
        FilePath = "../work/SampleQC_" + str(i) + ".txt"
        First = True
        for StrIn in open(FilePath, "r"):
            if First == True:
                First = False
                pass
            else:
                FileOut.write(StrIn)
FileOut.close()
