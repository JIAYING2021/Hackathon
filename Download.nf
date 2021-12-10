/* 
 *  Download FastQ files 
 */
Channel.of("SRR628582".."SRR628589")
       .set{ids}

process getfastq {
    publishDir "/home/ubuntu/Hackathon021", mode: 'copy'
    
    input:
    val SRAID from ids
    
    script:
    """
    wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1  
    fastq-dump --gzip --split-files ./${SRAID}.sra
    """
}

/* 
 *  Download Genome file
 */
Channel.of(1..22, "MT", "X", "Y")
       .set{chrs}

process getgenome {
    input:
    val chr from chrs

    output:
    file '*.fa.gz' into chromosomes
    
    script:
    """
    wget -O ${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    """
}

/* 
 *  Merge genome chromosome file
 */

process mergegenome {
    publishDir "/home/ubuntu/Hackathon021", mode: 'copy'


    input:
    file chromosome from chromosomes.flatten()

    script:
    """
    gunzip -c $chromosome > ref.fa

    """
}


/*
 *  Download human genome gnnotation file (gtf format)
 */

process getannot {
    
    script:
    """
    wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.101.chr.gtf.gz
    """
}
