# ---------------------------------------------------------------------------------------------
# Check for mandatory environment variables

# BOARD
if(NOT DEFINED BOARD)
  message(FATAL_ERROR "BOARD variable is not defined. Use -DBOARD=<board>\n\tValid boards: alveo_u200")
endif()

# Define a list of supported boards
set(VALID_BOARDS "alveo_u200")

set(SELECTED_BOARD "")

# Loop through the list of devices and check if one is enabled
foreach(VALID_BOARD ${VALID_BOARDS})

    if(BOARD STREQUAL VALID_BOARD)
        set(SELECTED_BOARD ${VALID_BOARD})
        break()
    endif()

endforeach()

if(SELECTED_BOARD STREQUAL "")
    message(FATAL_ERROR "Unknown board: ${BOARD}. Valid devices: ${VALID_BOARDS}")
endif()

# FPGA_CLOCK
if (NOT DEFINED FPGA_CLOCK)
  message(STATUS "FPGA_CLOCK variable is not defined. Using default value 200. Use -DFPGA_CLOCK=<clockMHz> to use a different value.")
  set(FPGA_CLOCK "200")
endif()

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Set basic configuration

set(TARGET_NAME "wfa_edit_alignment_fpga")
set(SOURCE_FILE "wfa_edit_alignment_fpga.c")

set(CMAKE_C_COMPILER "clang")
set(CMAKE_C_FLAGS "${CFLAGS} -fompss-2 -fompss-fpga-wrapper-code")
set(CMAKE_C_LINK_FLAGS "${LDFLAGS}")

set(EMULATION_FLAGS "-DFPGA_EMU")

set(DESIGN_FLAGS "-g -fompss-fpga-hls-tasks-dir ${CMAKE_BINARY_DIR}")

set(AIT_FLAGS__ "--name=${TARGET_NAME} --board=${BOARD} -c=${FPGA_CLOCK}")
set(AIT_FLAGS_DESIGN__ "--to_step=design")
set(AIT_FLAGS_D__ "--debug_intfs=both")

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Check for optional environment variables

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Set final AIT flags

# AIT flags
set(AIT_FLAGS_ "-fompss-fpga-ait-flags \"${AIT_FLAGS__}\"")
set(AIT_FLAGS_DESIGN_ "-fompss-fpga-ait-flags \"${AIT_FLAGS__} ${AIT_FLAGS_DESIGN__}\"")
set(AIT_FLAGS_DD_ "-fompss-fpga-ait-flags \"${AIT_FLAGS__} ${AIT_FLAGS_DESIGN__} ${AIT_FLAGS_D__}\"")
set(AIT_FLAGS_DB_ "-fompss-fpga-ait-flags \"${AIT_FLAGS__} ${AIT_FLAGS_D__}\"")

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# FPGA Emulation Target

set(PROGRAM_EMU "${TARGET_NAME}-emu")

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# FPGA Targets

# Host code targets
set(PROGRAM_P "${TARGET_NAME}-p")
set(PROGRAM_I "${TARGET_NAME}-i")
set(PROGRAM_D "${TARGET_NAME}-d")

# Design targets
set(PROGRAM_SEQ "${TARGET_NAME}-seq")
set(DESIGN_P "design-p")
set(DESIGN_I "design-i")
set(DESIGN_D "design-d")

# Bitstream targets
set(BITSTREAM_P "bitstream-p")
set(BITSTREAM_I "bitstream-i")
set(BITSTREAM_D "bitstream-d")

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# FPGA emulation compilations

# FPGA emulation
add_executable(${PROGRAM_EMU} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${PROGRAM_EMU} PROPERTIES COMPILE_FLAGS "${EMULATION_FLAGS}")

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Host code compilations

# Host code-p
add_executable(${PROGRAM_P} EXCLUDE_FROM_ALL ${SOURCE_FILE})

# Host code-i
add_executable(${PROGRAM_I} EXCLUDE_FROM_ALL ${SOURCE_FILE})

# Host code-d
add_executable(${PROGRAM_D} EXCLUDE_FROM_ALL ${SOURCE_FILE})

# Host code-seq
add_executable(${PROGRAM_SEQ} EXCLUDE_FROM_ALL ${SOURCE_FILE})

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Design compilations

# Design-p
add_executable(${DESIGN_P} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${DESIGN_P} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_DESIGN_}")
add_custom_command(
    TARGET ${DESIGN_P}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${DESIGN_P}
)

# Design-i
add_executable(${DESIGN_I} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${DESIGN_I} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_DESIGN_}")
add_custom_command(
    TARGET ${DESIGN_I}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${DESIGN_I}
)

# Design-d
add_executable(${DESIGN_D} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${DESIGN_D} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_DD_}")
add_custom_command(
    TARGET ${DESIGN_D}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${DESIGN_D}
)

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Bitstream compilations

# Bitstream-p
add_executable(${BITSTREAM_P} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${BITSTREAM_P} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_}")
add_custom_command(
    TARGET ${BITSTREAM_P}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${BITSTREAM_P}
)


# Bitstream-i
add_executable(${BITSTREAM_I} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${BITSTREAM_I} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_}")
add_custom_command(
    TARGET ${BITSTREAM_I}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${BITSTREAM_I}
)

# Bitstream-d
add_executable(${BITSTREAM_D} EXCLUDE_FROM_ALL ${SOURCE_FILE})
set_target_properties(${BITSTREAM_D} PROPERTIES COMPILE_FLAGS "${AIT_FLAGS_DB_}")
add_custom_command(
    TARGET ${BITSTREAM_D}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/${BITSTREAM_D}
)

# ---------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------
# Cleaning targets

# Extra clean
add_custom_target(extra_clean
    COMMAND ${CMAKE_COMMAND} -E rm -fv *.o ${TARGET_NAME}-? ${TARGET_NAME}_hls_automatic_clang.cpp ait_extracted.json
    COMMAND ${CMAKE_COMMAND} -E rm -rfv ${TARGET_NAME}_ait
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}

)

# ---------------------------------------------------------------------------------------------