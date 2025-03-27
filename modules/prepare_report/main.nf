process prepare_report {

    publishDir("$params.localdir"),  mode: 'copy'
    publishDir("$params.outdir"),  mode: 'copy'
    tag 'Preparing report'

    input:
    path report
    path report_utils

    output:
    path 'report.ipynb'
    path 'report_utils.py'

    script:
    """
    cp $report 'report.ipynb'
    cp $report_utils 'report_utils.py'
    """
}