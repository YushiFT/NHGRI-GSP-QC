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
chrom = str(1)

# define input files
vds_splitmulti_file = 'gs://rcstorage/matrixtable/' + chrom + '/splitmulti.vds'
pop_file = 'gs://rcstorage/population/ccdgf2_predicted_ethnicity_PC1-15.tsv'
lcr_file = 'gs://rcstorage/multidata/LCR-hs38.bed'

# define output files
qced_vds_file = 'gs://rcstorage/qced/' + chrom + '/ccdgf2_qced.vds'
variant_qc_table_file = 'gs://rcstorage/qced/' + chrom + '/ccdgf2_variant_qc_table.txt.gz'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. Remove LCR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("importing vds files")
vds = hl.read_matrix_table(vds_splitmulti_file)

print("removing lcr")
lcr = hl.import_bed(lcr_file, reference_genome='GRCh38')
vds = vds.filter_rows(hl.is_defined(lcr[vds.locus]), keep=False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III. Annotate samples with population
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("annotating ethnicity")
table = hl.import_table(pop_file, impute=True).key_by('Sample')
vds = vds.annotate_cols(**table[vds.s])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IV. Annotate variants first
#     Annotate variants with specific indicators, e.g. called rate.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hl.variant_qc(vds)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# V. Revising GT (set missing if):
#      DP > 400 or DP < 10
#      Heterozygous AND [ (AD alt/DP) < 20% OR PL ref < 20 OR (AD ref + AD alt) / DP < 90% ]
#      Homozygous ref AND [ GQ < 20 OR (AD ref / DP) < 0.9 ]
#      Homozygous alt AND [ PL ref < 20 OR (AD alt / DP) < 0.9 ]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("revising individual genotype")

vds = vds.filter_entries(
      (vds.DP > 400) |
      (vds.DP < 10) |
      (vds.GT.is_hom_ref() & ((vds.AD[0] / vds.DP < 0.9) | (vds.GQ < 20))) |
      (vds.GT.is_hom_var() & ((vds.AD[1] / vds.DP < 0.9) | (vds.PL[0] < 20) | ((vds.variant_qc.AF[0] > 0.5) & (vds.GQ < 20)))) |
      (vds.GT.is_het() & ( ((vds.AD[0] + vds.AD[1]) / vds.DP < 0.9) |
           (vds.AD[1] / vds.DP < 0.20) |
           (vds.PL[0] < 20))), keep = False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VI. Annotate variants again to update variants' annotations
#     Annotate variants with specific indicators, e.g. called rate.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hl.variant_qc(vds)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VII. Calculating HWE P-value and population specific called_rate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("calculating HWE p-value")
custompops = ['AFR', 'AMR', 'EAS', 'EUR', 'FIN', 'PUR', 'SAS']

hwelist = [hl.agg.filter(vds.Pop == pop, hl.agg.hardy_weinberg_test(vds.GT))['p_value'] for pop in custompops]
vds = vds.annotate_rows(
   hwe=hl.Struct(**dict(zip(custompops, hwelist))),
   lowestphwe = hl.min(hl.array(hwelist))
 )

print("calculating population-specific called rate")
#callratelist = [hl.agg.fraction(hl.agg.filter(vds.Pop == pop, hl.is_defined(vds.GT))) for pop in custompops]
callratelist = [hl.agg.filter(vds.Pop == pop, hl.agg.fraction(hl.is_defined(vds.GT))) for pop in custompops]

vds = vds.annotate_rows(
   callrate=hl.Struct(**dict(zip(custompops, callratelist))),
   lowestcallrate = hl.min(hl.array(callratelist))
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VIII. Annotate filtered variants
#    Remove variants with specific indicators, e.g. called rate.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print("filtering variants")

def MAF(AF_list):
    AF_list = AF_list.append(1-hl.sum(AF_list))
    return(hl.sorted(AF_list)[-2])

# define conditions and thresholds
step0 = vds.variant_qc.AC[1] > 0
step1 = hl.len(vds.filters) == 0
step2 = vds.info.QD > 4
step3 = vds.lowestcallrate >= 0.98
step4 = ( (vds.lowestphwe >= 1e-5) & (MAF(vds.info.AF) > 0.01) ) | ( (vds.lowestphwe >= 1e-6) & (MAF(vds.info.AF) <= 0.01) )

qc_struct = hl.Struct(step0 = step0, step1 = step1, step2 = step2, step3=step3, step4=step4)

qccumul = hl.Struct(step0 = step0,
                    step1 = step0 & step1,
                    step2 = step0 & step1 & step2,
                    step3 = step0 & step1 & step2 & step3,
                    step4 = step0 & step1 & step2 & step3 & step4)

# Annotate rows with PASS indicator
vds = vds.annotate_rows(qc = qc_struct, qccum = qccumul, qcpass = qccumul.step4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IX. Downcoding
#     Only keep GT for individual column
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#print("downcoding...")
#vds = vds.select_entries(vds.GT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X. Write output VDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#print("writing out vds...")
#vds.write(qced_vds_file, overwrite=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# XIII. write variants table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("write variants table...")
vds.rows().select('info', 'qc', 'qccum', 'qcpass', 'variant_qc', 'filters', 'hwe', 'lowestphwe', 'callrate', 'lowestcallrate' ).flatten().export(variant_qc_table_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print Runtime
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stop = timeit.default_timer()

print("runtime: " + str(stop - start) + " seconds")
