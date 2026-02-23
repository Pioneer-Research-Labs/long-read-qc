process aggregate_sample_sheets {
    publishDir "$params.outdir", mode: 'copy'
    tag 'Aggregating sample sheets'

    input:
    path sample_sheet_map

    output:
    path 'aggregated_sample_sheet.csv'

    script:
    """
    aggregate_sample_sheets.py $sample_sheet_map
    """
}