from hail import *
hc = HailContext(tmp_dir = '/nvme/tmp/shusson')

vds = hc.read('/nvme/marpin/MGRB_phase2_hs37d5x/GATK3_fastq2gvcf-hs37d5x-1.0_split-combine_genotype_vqsr_hail/MGRB.phase2.SNPtier12.match.vqsr.minrep.locusannot.WGStier12.unrelated.tgp.hrc.gnomad.dbsnp.clinvar.cato.eigen.vep.vds')

vds = vds.split_multi()

vds = vds.annotate_variants_expr([
        'va.allele = va.alleles[va.aIndex-1]',
        'va.extra.wasSplit = va.wasSplit',
        'va.extra.nHomVarFemale = gs.filter(s => sa.pheno.isFemale).filter(g => g.isHomVar()).count()',
        'va.extra.nHomVarMale = gs.filter(s => !sa.pheno.isFemale).filter(g => g.isHomVar()).count()',
        'va.extra.nHomRefFemale = gs.filter(s => sa.pheno.isFemale).filter(g => g.isHomRef()).count()',
        'va.extra.nHomRefMale = gs.filter(s => !sa.pheno.isFemale).filter(g => g.isHomRef()).count()'
    ]).annotate_variants_expr([
        'va.extra.nHomVar = va.extra.nHomVarMale + va.extra.nHomVarFemale',
        'va.extra.nHomRef = va.extra.nHomRefMale + va.extra.nHomRefFemale'
    ])

vds = vds.filter_variants_expr('va.allele.metrics.allele_counts.total.toInt() <= 0', keep=False)

vds = vds.annotate_variants_expr('''
    va.extra.nHet = {
        total: 
            if (v.inXNonPar())
                (2 * (va.allele.metrics.allele_counts.male - va.extra.nHomVarMale)) + (va.allele.metrics.allele_counts.female - (2 * va.extra.nHomVarFemale))
            else if (v.inYNonPar())
                0
            else
                va.allele.metrics.allele_counts.total - (2 * va.extra.nHomVar)   
            }
    ''')

vds.write('/nvme/tmp/shusson/sgc.vds')
