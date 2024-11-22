# ---------------------------------------------------------------------------------------------
# Set basic configuration

set(TARGET_NAME "wfa_edit_alignment_cpu")
set(SOURCE_FILE "wfa_edit_alignment_cpu.c")

set(CMAKE_C_COMPILER "clang")
set(CMAKE_C_FLAGS "${CFLAGS} -Wall -Wextra -Werror")
set(CMAKE_C_LINK_FLAGS "${LDFLAGS}")

# ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# CPU Targets

add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${SOURCE_FILE})

# ---------------------------------------------------------------------------------------------