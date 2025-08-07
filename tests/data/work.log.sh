../../bin/basevar caller --mapq=10 --min-af=0.05 --batch-count=1 --thread=1 --regions=CHROMOSOME_I:900-1200 --pop-group=sample_group.info --output-vcf vz.vcf --output-cvg t.cvg -R ce.fa.gz -I range.bam -I range.cram -I range.bam --filename-has-samplename

#chr11:5246595-5248428,chr13:32890633-32972781,chr16:222869-227506,chr17:41197764-41276135
# long options
../../bin/basevar caller --min-BQ=0 --min-mapq=10 --min-af=0.05 --batch-count=20 --thread=4 --regions=chr11:5246595-5248428,chr17:41197764-41276135 --pop-group=sample_group.info --output-vcf tt.vcf.gz -R ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz -L bam90.list bam100/00alzqq6jw.bam bam100/09t3r9n2rg.bam bam100/0fkpl1p55b.bam bam100/13dg1gvsfk.bam bam100/17phildszl.bam bam100/1dbpgqt0dq.bam bam100/1kyws27hoc.bam bam100/1ych8rmufr.bam bam100/4e56w6ezsx.bam bam100/51rwla2fps.bam > log2

# short options
ulimit -n 100000 && ../../bin/basevar caller -Q 0 -q 10 -m 0.05 -B 20 -t 4 -r chr11:5246595-5248428,chr17:41197764-41276135 -G sample_group.info --output-vcf tt.vcf -R ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz -L bam90.list bam100/00alzqq6jw.bam bam100/09t3r9n2rg.bam bam100/0fkpl1p55b.bam bam100/13dg1gvsfk.bam bam100/17phildszl.bam bam100/1dbpgqt0dq.bam bam100/1kyws27hoc.bam bam100/1ych8rmufr.bam bam100/4e56w6ezsx.bam bam100/51rwla2fps.bam

# For NA12878
../../bin/basevar caller -Q 0 -q 10 -m 0.05 -B 20 -t 4 --output-vcf tt.vcf -R ../../learn_from/GLIMPSE/tutorial/reference_genome/hs38DH.chr22.fa.gz ../../learn_from/GLIMPSE/tutorial/NA12878_1x_bam/NA12878.bam
