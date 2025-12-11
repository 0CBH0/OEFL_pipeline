# Handle OE, some intermediate files need to be copied.
cellranger count --id=OE --transcriptome=~/library/refdata-gex-mm10-2020-A --fastqs ~/raw_data/OE_10X --sample=OE --localcores=32 --localmem=128
Rscript calcReads_oe.R
nextflow run epi2me-labs/wf-single-cell --expected_cells 15000 --bam OEfulllength.bam --kit multiome:v1 --ref_genome_dir ~/library/refdata-gex-mm10-2020-A
python3 isoquant.py --reference ~/library/refdata-gex-mm10-2020-A/fasta/genome.fa --genedb ~/library/refdata-gex-mm10-2020-A/genes/genes.gtf --bam oe_full_length/output/OEfulllength/tagged.bam --data_type nanopore --read_group tag:CB -o oe_full_length/output/OEfulllength/isoquant_bam --threads 24
python3 iso_proc_oe.py

# Handle OSN, some intermediate files need to be copied.
cellranger count --id=OSN --transcriptome=~/library/refdata-gex-mm10-2020-A --fastqs ~/raw_data/OSN_10X --sample=OSN --localcores=32 --localmem=128
Rscript calcReads_osn.R
nextflow run epi2me-labs/wf-single-cell --expected_cells 15000 --bam OSNfulllength.bam --kit multiome:v1 --ref_genome_dir ~/library/refdata-gex-mm10-2020-A
python3 isoquant.py --reference ~/library/refdata-gex-mm10-2020-A/fasta/genome.fa --genedb ~/library/refdata-gex-mm10-2020-A/genes/genes.gtf --bam oe_full_length/output/OSNfulllength/tagged.bam --data_type nanopore --read_group tag:CB -o oe_full_length/output/OSNfulllength/isoquant_bam --threads 24
python3 iso_proc_osn.py

