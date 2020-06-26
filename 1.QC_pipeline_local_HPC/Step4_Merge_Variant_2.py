"""
Write out merge variant result
"""

with open("../temp/VariantQC.txt", "w") as FileOut:

    Header = ['CHROM', 'POS', 'ID',
              'REF', 'ALT', 'QUAL',
              'Coverage', 'nHet', 'nHomRef',
              'nHomVar', 'nCalled', 'nNotCalled', '\n']
    FileOut.write('\t'.join(Header))

    First = True
    for StrIn in open("../temp/temp.txt"):
        if First == True:
            First = False
            pass
        else:
            TempStr = StrIn.split('\t')

            nHet = int(TempStr[7])
            for i in range(13,115,6):
                nHet += int(TempStr[i])

            nHomRef = int(TempStr[8])
            for i in range(14, 116, 6):
                nHomRef += int(TempStr[i])

            nHomVar = int(TempStr[9])
            for i in range(15, 117, 6):
                nHomVar += int(TempStr[i])

            nCalled = int(TempStr[10])
            for i in range(16, 118, 6):
                nCalled += int(TempStr[i])

            nNotCalled = int(TempStr[11])
            for i in range(17, 119, 6):
                nNotCalled += int(TempStr[i])

            StrOut = TempStr[0:7]
            StrOut.append( str(nHet) )
            StrOut.append( str(nHomRef) )
            StrOut.append( str(nHomVar) )
            StrOut.append( str(nCalled) )
            StrOut.append( str(nNotCalled) )
            StrOut.append( '\n' )
            FileOut.write( '\t'.join(StrOut) )
FileOut.close()
