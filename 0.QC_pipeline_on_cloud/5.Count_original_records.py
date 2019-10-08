import hail as hl
import random
import timeit
import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. Initialization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start = timeit.default_timer()

hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')


#chrom = str(sys.argv[1])
chrom = str(22)
# define input files
#filelist = ['gs://rcstorage/genotype/gnarly_chr22.' + str(i) + '.variant_filtered.vcf.gz' for i in range(1670)]
new_file = ['gs://rcstorage/variantslist/CCDG.sites_only.filtered.vcf.gz']
lcr_file = 'gs://rcstorage/multidata/LCR-hs38.bed'

# define output file
num_file_out = ( '../YushiFTang/count/' + chrom + 'count.txt')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. Import VCF containing Filter information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("For chromosome 22, the number of original records is: ")
vds_pre = hl.import_vcf(new_file, force_bgz=True, reference_genome='GRCh38')
vds_pre = vds_pre.filter_rows(vds_pre.locus.contig == 'chr22', keep=True)


num01 = vds_pre.count()
print(num01)



print("For chromosome 22, the number of records fail to pass before spliting multi is: ")
vds_fail_pass_pre = vds_pre.filter_rows(hl.len(vds_pre.filters)==0, keep=False)
num02 = vds_fail_pass_pre.count()
print(num02)



print("For chromosome 22, the number of records pass before spliting multi is: ")
vds_pass_pre = vds_pre.filter_rows(hl.len(vds_pre.filters)==0, keep=True)
num03 = vds_pass_pre.count()
print(num03)



print("removing lcr")
print("after removing lcr, the number of variants is: ")
lcr = hl.import_bed(lcr_file, reference_genome='GRCh38')
vds_pass_pre = vds_pass_pre.filter_rows(hl.is_defined(lcr[vds_pass_pre.locus]), keep=False)
num04 = vds_pass_pre.count()
print(num04)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. Split multi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("spliting multi...")
print("after spliting multi, the number of variants with PASS label is: ")
vds_pass_pre = hl.split_multi_hts(vds_pass_pre)
num05 = vds_pass_pre.count()
print(num05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. Import VCF
#     Combine all VCF chunks for one chromosome and import as vds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#print("importing vcf...")
#print("the number of records before spliting multi is: ")

#vds = hl.import_vcf(filelist, force_bgz=True, reference_genome='GRCh38')

#num1 = vds.count()

#print(num1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III-I. Count variants without PASS in Filter column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#print("filtering variants without pass...")
#print("the number of records fail to pass before spliting multi is: ")

#vds_fail_pass = vds.filter_rows(hl.len(vds.filters)==0, keep=False)

#num2 = vds_fail_pass.count()

#print(num2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III-II. Count variants with PASS in Filter column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#print("filtering variants without pass...")
#print("the number of records fail to pass before spliting multi is: ")

#vds_pass = vds.filter_rows(hl.len(vds.filters)==0, keep=True)

#num3 = vds_pass.count()

#print(num3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. Split multi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#print("spliting multi...")
#vds = hl.split_multi_hts(vds.select_entries(vds.GT))

#num2 = vds.count()

#print("the number of records after spliting multi is: ")

#print(num2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print Runtime
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stop = timeit.default_timer()

print("runtime: " + str(stop - start) + " seconds")
