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
