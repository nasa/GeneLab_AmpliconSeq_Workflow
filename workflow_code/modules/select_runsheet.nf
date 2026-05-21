process SELECT_RUNSHEET {
    tag "Selecting runsheet for target region: ${target_region}"

    input:
        path("runsheets/*")
        val(target_region)
        val(specify_runsheet)

    output:
        path("*.csv"), emit: runsheet

    script:
    def specify = specify_runsheet ? "--specify-runsheet ${specify_runsheet}" : ""
    """
    select_runsheet.py --input-dir runsheets --target ${target_region} ${specify}
    """
}