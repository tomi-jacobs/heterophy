process GENERATE_REPORT {
    tag "Generating report"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path inputs   // coverage tables, SNP stats, comparison JSON

    output:
    path "heterophy_report.html", emit: report

    script:
    """
    generate_report.py \\
        --outdir       . \\
        --title        "${params.report_title}" \\
        --pipeline_ver "1.0.0"
    """
}
