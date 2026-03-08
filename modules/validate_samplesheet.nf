/*
========================================================================================
    MODULE: validate_samplesheet
========================================================================================
*/
process VALIDATE_SAMPLESHEET {
    tag "Validating samplesheet"
    label 'process_low'

    conda "envs/python.yml"
    container "heterophy/python:1.0.0"

    input:
    path samplesheet

    output:
    path "validated_samplesheet.csv", emit: validated

    script:
    """
    validate_samplesheet.py \\
        --input ${samplesheet} \\
        --output validated_samplesheet.csv
    """
}
