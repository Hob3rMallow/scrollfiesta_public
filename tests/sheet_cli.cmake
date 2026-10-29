# Removed experiment commands must fail as usage errors, before creating runs.
file(MAKE_DIRECTORY "${WORK_DIR}")
execute_process(COMMAND "${SHEET_ASSEMBLE}" --help
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT rc EQUAL 0 OR NOT "${out}${err}" MATCHES "configs/default.json")
    message(FATAL_ERROR "Public help failed: ${rc}: ${out}${err}")
endif()
if("${out}${err}" MATCHES "--[a-z-]+-case|--flatten-probe|--sparse-benchmark")
    message(FATAL_ERROR "Experimental commands remain in public help")
endif()
function(expect_usage)
    execute_process(COMMAND "${SHEET_ASSEMBLE}" ${ARGN}
                    RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
    if(NOT rc EQUAL 2 OR EXISTS "${WORK_DIR}/unexpected")
        message(FATAL_ERROR "Removed/unknown option was not rejected before IO: ${ARGN}: ${rc}")
    endif()
endfunction()
foreach(option --bottom-up-case --bottom-up-repair-case --required-processes-case
               --required-processes-refine-case --required-processes-aligned-case)
    expect_usage("${option}" missing.asr "${WORK_DIR}/unexpected" 1800 100 1)
endforeach()
expect_usage(--bottom-up-attempt-case missing.asr "${WORK_DIR}/unexpected" 0 1)
expect_usage(--repair-sequence-case missing.asr "${WORK_DIR}/unexpected" missing.plan 900)
expect_usage(--repair-case missing.bin "${WORK_DIR}/unexpected" 1)
expect_usage(--min-piece-case missing.asr "${WORK_DIR}/unexpected" 100)
expect_usage(--flatten-probe missing.vmesh 0)
expect_usage(--unknown)
expect_usage(missing-pile "${WORK_DIR}/unexpected" --stop-after continuity)
expect_usage(missing-pile "${WORK_DIR}/unexpected" --stop-after source)
# Exercise the shipped config through the real parser. Missing input must be
# reported after config loading, including when the default is implicit.
foreach(explicit_config FALSE TRUE)
    set(config_arg)
    if(explicit_config)
        set(config_arg "${SOURCE_DIR}/configs/default.json")
    endif()
    execute_process(COMMAND "${SHEET_ASSEMBLE}" "${WORK_DIR}/missing-pile"
                    "${WORK_DIR}/empty-run" ${config_arg} --stop-after clean
                    WORKING_DIRECTORY "${SOURCE_DIR}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 3 OR NOT "${out}${err}" MATCHES "no per-cube meshes")
        message(FATAL_ERROR "Public config did not load: ${rc}: ${out}${err}")
    endif()
endforeach()
