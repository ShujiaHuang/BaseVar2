../../bin/basevar caller -f ce.fa.gz -o test_p01.vcf -B 50 -t 2 range.bam
../../bin/basevar caller --min-mapq=10 --min-af=0.05 --batch-count=1 --thread=1 --regions=CHROMOSOME_I:900-1200 --pop-group=sample_group.info --output vz.vcf -f ce.fa.gz --filename-has-samplename range.bam range.cram range.bam

#chr11:5246595-5248428,chr13:32890633-32972781,chr16:222869-227506,chr17:41197764-41276135
# short options
../../bin/basevar caller -Q 0 -q 10 -m 0.05 -B 10 -t 4 -r chr11:5246595-5248428,chr17:41197764-41276135 -G sample_group.info -o tt.vcf -f ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz -L bam90.list bam100/00alzqq6jw.bam bam100/09t3r9n2rg.bam bam100/0fkpl1p55b.bam bam100/13dg1gvsfk.bam bam100/17phildszl.bam bam100/1dbpgqt0dq.bam bam100/1kyws27hoc.bam bam100/1ych8rmufr.bam bam100/4e56w6ezsx.bam bam100/51rwla2fps.bam

# Test basevar caller by using synthetic data
python synthetic/generate.py --output-dir synthetic_normal_test --skip-cram 2>&1
../../bin/basevar caller -f synthetic_normal_test/ref/mini_ref.fa -o out2.vcf -G synthetic_normal_test/samples_group.info synthetic_normal_test/bam/*.bam
python synthetic/evaluate.py --vcf out2.vcf --truth synthetic_normal_test/ground_truth_variants.tsv

../../bin/basevar caller -f synthetic/ref/mini_ref.fa -o out.vcf -G synthetic/samples_group.info synthetic/bam/*.bam
#../../bin/basevar caller -f synthetic/ref/mini_ref.fa -o out_c.vcf -G synthetic/samples_group.info synthetic/cram/*.cram
python synthetic/evaluate.py --vcf out.vcf --truth synthetic/ground_truth_variants.tsv





# subsam
../../bin/basevar subsam -i tt.vcf~ -o t.vcf 5ffp4ybnks 9jikb1nr7d
le tt.vcf\~ | cut -f 1-10,18-18 |le


# bcftools test
bcftools mpileup -Q 0 -q 10 -f ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz bam100/00alzqq6jw.bam bam100/09t3r9n2rg.bam bam100/0fkpl1p55b.bam bam100/13dg1gvsfk.bam bam100/17phildszl.bam bam100/1dbpgqt0dq.bam bam100/1kyws27hoc.bam bam100/1ych8rmufr.bam bam100/4e56w6ezsx.bam bam100/51rwla2fps.bam | bcftools call -mv | le

# For NA12878
../../bin/basevar caller -Q 0 -q 10 -m 0.05 -B 10 -t 4 --output tt.vcf -f ../../learn_from/GLIMPSE/tutorial/reference_genome/hs38DH.chr22.fa.gz ../../learn_from/GLIMPSE/tutorial/NA12878_1x_bam/NA12878.bam
samtools tview --reference ../../learn_from/GLIMPSE/tutorial/reference_genome/hs38DH.chr22.fa.gz ../../learn_from/GLIMPSE/tutorial/NA12878_1x_bam/NA12878.bam

# motif
../../bin/basevar motif -o ttt -f ~/Projects/BaseVar/tests/data/hg19.NC_012920.fasta.gz -L bam90.list bam100/00alzqq6jw.bam bam100/09t3r9n2rg.bam bam100/0fkpl1p55b.bam bam100/13dg1gvsfk.bam bam100/17phildszl.bam bam100/1dbpgqt0dq.bam bam100/1kyws27hoc.bam bam100/1ych8rmufr.bam bam100/4e56w6ezsx.bam bam100/51rwla2fps.bam -t 10 > log




