
/*
 *  Download Genome file
 */
Channel.of(1..22, "Mt", "X", "Y")
       .set{chrs}

process getgenome {
    input:
    val chr from chrs

    script:
    """
    wget -O ${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    gunzip -c *.fa.gz > ref.fa
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
