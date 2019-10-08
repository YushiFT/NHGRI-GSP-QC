import hail as hl
import random
import timeit
import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. Initialization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start = timeit.default_timer()

hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')


chrom = str(1)

# define input files
vds_splitmulti_file = 'gs://rcstorage/matrixtable/' + chrom + '/splitmulti.vds'
lcr_file = 'gs://rcstorage/multidata/LCR-hs38.bed'
variant_list_file = 'gs://rcstorage/qced/'+chrom+'/qced_'+chrom+'_variant_list.txt'


# define output files
sample_qc_info_postqc_file = 'gs://rcstorage/qced/' + chrom + '/sample_qc_info_postqc_revisegt.txt'


print("importing vds files...")
vds = hl.read_matrix_table(vds_splitmulti_file)
num0 = vds.count()
print(num0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. Remove LCR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("removing lcr...")
lcr = hl.import_bed(lcr_file, reference_genome='GRCh38')
vds = vds.filter_rows(hl.is_defined(lcr[vds.locus]), keep=False)
num1 = vds.count()
print(num1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III. Annotate variants with PASS or FAIL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("annotating variants...")

variant_list = hl.import_table(variant_list_file)

variant_list_post = variant_list.add_index()
variant_list_post = variant_list_post.key_by('idx')
vds_post = vds.add_row_index()

vds_post = vds_post.annotate_rows(**variant_list_post[vds_post.row_idx])

num2 = vds_post.count()
print(num2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. Filtering variants without PASS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("removing variants...")
vds_post = vds_post.filter_rows(vds_post.label == 'PASS', keep=True)

num3 = vds_post.count()
print(num3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# V. Annotate variants again to update variants' annotations
#     Annotate variants with specific indicators, e.g. called rate.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#vds_post = hl.variant_qc(vds_post)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VI. Revising GT (set missing if):
#      DP > 400 or DP < 10
#      Heterozygous AND [ (AD alt/DP) < 20% OR PL ref < 20 OR (AD ref + AD alt) / DP < 90% ]
#      Homozygous ref AND [ GQ < 20 OR (AD ref / DP) < 0.9 ]
#      Homozygous alt AND [ PL ref < 20 OR (AD alt / DP) < 0.9 ]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("revising individual genotype...")

vds_post = vds_post.filter_entries(
           (vds_post.DP > 400) |
           (vds_post.DP < 10) |
           (vds_post.GT.is_hom_ref() & ((vds_post.AD[0] / vds_post.DP < 0.9) | (vds_post.GQ < 20))) |
           (vds_post.GT.is_hom_var() & ((vds_post.AD[1] / vds_post.DP < 0.9) | (vds_post.PL[0] < 20) | (vds_post.GQ < 20))) |
           (vds_post.GT.is_het() & ( ((vds_post.AD[0] + vds_post.AD[1]) / vds_post.DP < 0.9) | (vds_post.AD[1] / vds_post.DP < 0.20) | (vds_post.PL[0] < 20))), keep = False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VI. Performing sample QC on remaining variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("sample QC...")
vds_post = hl.sample_qc(vds_post)

print("writing sample QC results...")
vds_post.cols().flatten().export(sample_qc_info_postqc_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print Runtime
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stop = timeit.default_timer()

print("runtime: " + str(stop - start) + " seconds")
