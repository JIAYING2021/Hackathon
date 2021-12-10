/* 
 * Pipeline input parameters 
 */
params.reads = "$baseDir/*_{1,2}.fastq.gz"
params.genome = "$baseDir/ref.fa"
params.annot = "$baseDir/Homo_sapiens.GRCh38.101.chr.gtf"

log.info """\
         R N A S E Q   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         genome       : ${params.genome}
         annotation   : ${params.annot}
         """
         .stripIndent()

/* 
 * Define the `index` process that create a binary index 
 * given the genome file
 */
process index {
    
    input:
    path genome from ref
    
    output:
    path 'index' into index_ch

    script: 
    """
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir index --genomeFastaFiles $genome
    """
}

/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */

Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch } 

/* 
 * Maps each read-pair by using STAR mapper tool
 */

process mapping {
    tag "$pair_id"
         
    input:
    path indexdir from index_ch
    tuple pair_id, path(reads) from read_pairs_ch
    
    output:
    tuple val(pair_id), file("${pair_id}.bam") into {bam_ch; bam_ch2}
 
    script:
    """
    STAR --outSAMstrandField intronMotif \
    --outFilterMismatchNmax 4 \
    --outFilterMultimapNmax 10 \
    --genomeDir $indexdir \
    --readFilesCommand zcat
    --readFilesIn ${reads[0]} ${reads[1]} \
    --runThreadN ${task.cpus} \
    --outSAMunmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM 36949420170 \
    > ${pair_id}.bam
    """
}

/* 
 * Bam file indexing
 */

process bamindex {
    tag "BAM file indexing on $bam_id"

    input:
    val bam_id from bam_ch

    output:
    file '*bam.bai' into bai_ch
    
    script:
    """
    samtools index ${bam_id}.bam ${bam_id}.bam.bai
    """
}  

/*
 *  Quantification of gene expression levels by using the "subread" tool
 */
process quantification {
    tag "Quantification of gene expression"
    
    input:
    path annot from params.annot
    file bam_file from bam_ch2.flatten()
    file bai from bai_ch
    
    output:
    file("count.txt") into count_matrix
    
    script:
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a $annot -o count.txt $bam_file 
    """
} 

/*
 *  Extraction of gene count
 */

process extraction {
   tag "Extraction of gene counts"
    input:
    file count from count_matrix

    output:
    file("counts.txt") into counts
    
    script:
    """
    sed -i "1d" $count > count1.txt
    cat count1.txt | cut -f 1, 7-14 > counts.txt
    """ 
}

/*
 *  Identification of differentially expressed genes (DEGs) using DESeq2 in R studio
 */

workflow.onComplete { 
 log.info ( workflow.success ? "\n Done! Congratulations! \n" : "Oops ... something went wrong!" )
}

