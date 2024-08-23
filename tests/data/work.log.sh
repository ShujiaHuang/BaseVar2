../../tests/basevar basetype --mapq=10 --min-af=0.05 --batch-count=50 --thread=8 --regions=ref1,ref2:5,ref3:10-20 --pop-group=sample_group.info --output-vcf t.vcf --output-cvg t.cvg --smart-rerun -R tinyfasta.fa -L bam100.list

../../tests/basevar basetype --mapq=10 --min-af=0.05 --batch-count=15 --thread=1 --regions=ref1,ref2:5,ref3:10-20 --pop-group=sample_group.info --output-vcf t.vcf --output-cvg t.cvg --smart-rerun -R tinyfasta.fa -L bam90.list -I bam100/17phildszl.bam -I bam100/1dbpgqt0dq.bam -I bam100/1kyws27hoc.bam

../../tests/basevar basetype --mapq=10 --min-af=0.05 --batch-count=15 --thread=1 --regions=CHROMOSOME_I:900-1200 --pop-group=sample_group.info --output-vcf t.vcf.gz --output-cvg t.cvg.gz --smart-rerun -R ce.fa.gz -I range.bam
../../tests/basevar basetype --mapq=10 --min-af=0.05 --batch-count=15 --thread=1 --regions=CHROMOSOME_I:900-1200,CHROMOSOME_IV:900-1200 --pop-group=sample_group.info --output-vcf t.vcf.gz --output-cvg t.cvg.gz --smart-rerun -R ce.fa.gz -I range.bam



