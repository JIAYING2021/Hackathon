/*
 *  Download FastQ files
 */
Channel.of("SRR628582".."SRR628589")
       .set{ids}

process getfastq {
    input:
    val SRAID from ids

    script:
    """
    wget -O ${SRAID}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${SRAID}/${SRAID}.1
    fastq-dump --gzip --split-files ./${SRAID}.sra
    """
}
