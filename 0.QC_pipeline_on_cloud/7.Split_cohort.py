import hail as hl
import random
import timeit
import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I. Initialization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

start = timeit.default_timer()

hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')


chrom = str(21)

# define input files
vds_splitmulti_file = 'gs://rcstorage/matrixtable/' + chrom + '/splitmulti.vds'
lcr_file = 'gs://rcstorage/multidata/LCR-hs38.bed'
variant_list_file = 'gs://rcstorage/qced/'+chrom+'/qced_'+chrom+'_variant_list.txt'
mark_file = 'gs://rcstorage/multidata/sample_marker.tsv'

# define output files
qced_vcf_baylor = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_baylor.vcf.bgz'
qced_vcf_baylor_aric = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_baylor_aric.vcf.bgz'
qced_vcf_baylor_sol = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_baylor_sol.vcf.bgz'
qced_vcf_broad = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_broad.vcf.bgz'
qced_vcf_broad_aflmu = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_broad_aflmu.vcf.bgz'
qced_vcf_broad_taichi = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_broad_taichi.vcf.bgz'
qced_vcf_broad_virgo = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_broad_virgo.vcf.bgz'
qced_vcf_nygc = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_nygc.vcf.bgz'
qced_vcf_nygc_galaii = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_nygc_galaii.vcf.bgz'
qced_vcf_nygc_ssc = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_nygc_ssc.vcf.bgz'
qced_vcf_washu = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_washu.vcf.bgz'
qced_vcf_main = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_main.vcf.bgz'

#qced_vcf_washu_cleveland = 'gs://rcstorage/qced_cohort/' + chrom + '/ccdgf2_qced_washu_cleveland.vcf.bgz'


print("importing vds files...")
vds = hl.read_matrix_table(vds_splitmulti_file)
num0 = vds.count()
print(num0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. Mark samples with center or studies
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("annotating center...")
table = hl.import_table(mark_file, impute=True).key_by('Sample')
vds = vds.annotate_cols(**table[vds.s])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III. Remove LCR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("removing lcr...")
lcr = hl.import_bed(lcr_file, reference_genome='GRCh38')
vds = vds.filter_rows(hl.is_defined(lcr[vds.locus]), keep=False)
num1 = vds.count()
print(num1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. Annotate variants with PASS or FAIL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("annotating variants...")

variant_list = hl.import_table(variant_list_file)

variant_list_post = variant_list.add_index()
variant_list_post = variant_list_post.key_by('idx')
vds = vds.add_row_index()

vds = vds.annotate_rows(**variant_list_post[vds.row_idx])

num2 = vds.count()
print(num2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# V. Filtering variants without PASS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("removing variants...")
vds = vds.filter_rows(vds.label == 'PASS', keep=True)

num3 = vds.count()
print(num3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VI. Downcoding
#     Only keep GT for individual column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("downcoding...")
vds = vds.select_entries(vds.GT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VII. Split into small cohorts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds_baylor = vds.filter_cols(vds.Mark == 'BAYLOR', keep=True)
vds_baylor_aric = vds.filter_cols(vds.Mark == 'BAYLOR_ARIC', keep=True)
vds_baylor_sol = vds.filter_cols(vds.Mark == 'BAYLOR_SOL', keep=True)
vds_broad = vds.filter_cols(vds.Mark == 'BROAD', keep=True)
vds_broad_aflmu = vds.filter_cols(vds.Mark == 'BROAD_AFLMU', keep=True)
vds_broad_taichi = vds.filter_cols(vds.Mark == 'BROAD_TAICHI', keep=True)
vds_broad_virgo = vds.filter_cols(vds.Mark == 'BROAD_VIRGO', keep=True)
vds_nygc = vds.filter_cols(vds.Mark == 'NYGC', keep=True)
vds_nygc_galaii = vds.filter_cols(vds.Mark == 'NYGC_GALAII', keep=True)
vds_nygc_ssc = vds.filter_cols(vds.Mark == 'NYGC_SSC', keep=True)
vds_washu = vds.filter_cols((vds.Mark == 'WASHU') | (vds.Mark == 'WASHU_CLEVELAND'), keep=True)
#vds_washu_cleveland = vds.filter_cols(vds.Mark == 'WASHU_CLEVELAND', keep=True)
vds_main = vds.filter_cols((vds.Mark == 'BAYLOR') | (vds.Mark == 'BROAD') | (vds.Mark == 'NYGC') | (vds.Mark == 'WASHU') | (vds.Mark == 'WASHU_CLEVELAND'), keep=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VIII. Write output VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print("writing out vcf...")
hl.export_vcf(vds_baylor, qced_vcf_baylor)
hl.export_vcf(vds_baylor_aric, qced_vcf_baylor_aric)
hl.export_vcf(vds_baylor_sol, qced_vcf_baylor_sol)
hl.export_vcf(vds_broad, qced_vcf_broad)
hl.export_vcf(vds_broad_aflmu, qced_vcf_broad_aflmu)
hl.export_vcf(vds_broad_taichi, qced_vcf_broad_taichi)
hl.export_vcf(vds_broad_virgo, qced_vcf_broad_virgo)
hl.export_vcf(vds_nygc, qced_vcf_nygc)
hl.export_vcf(vds_nygc_galaii, qced_vcf_nygc_galaii)
hl.export_vcf(vds_nygc_ssc, qced_vcf_nygc_ssc)
hl.export_vcf(vds_washu, qced_vcf_washu)
#hl.export_vcf(vds_washu_cleveland, qced_vcf_washu_cleveland)
hl.export_vcf(vds_main, qced_vcf_main)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print Runtime
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stop = timeit.default_timer()

print("runtime: " + str(stop - start) + " seconds")
