"""Call peaks using MACS2 or Epic2"""

rule call_peaks:
    input:
        bam = "align/{genome}/{aligner}/{id}/aligned.bam"
    output:
        bed = "peaks/{genome}/{aligner}/{peakcaller}/{id}/peaks.bed"
    params:
        gs = lambda wc: config['genome_size'][wc.genome],
        nomodel = lambda wc: '--nomodel' if wc.peakcaller == 'macs2' else ''
    log:
        "logs/call_peaks/{genome}_{aligner}_{peakcaller}_{id}.log"
    conda:
        "envs/peak.yml"
    shell:
        """
        if [ "{wildcards.peakcaller}" = "macs2" ]; then
            macs2 callpeak -t {input.bam} -f BAMPE -g {params.gs} {params.nomodel} -n {wildcards.id} --outdir $(dirname {output.bed})
            mv $(dirname {output.bed})/{wildcards.id}_peaks.narrowPeak {output.bed}
        else
            epic2 -t {input.bam} -c {input.bam} -o {output.bed}
        fi
        """ > {log} 2>&1
