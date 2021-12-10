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
    input:
    file chr from chromosomes.flatten()

    script:
    """
    gunzip -c $chr > ref.fa
    """
}
