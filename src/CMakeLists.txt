cmake_minimum_required(VERSION 3.00)
project(WFA-EditDistance)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_VERBOSE_MAKEFILE ON)

# If variable DEVICE is not defined, show an error message
if(NOT DEFINED DEVICE)
  message(FATAL_ERROR "DEVICE variable is not defined. Use -DDEVICE=<device>\n\tValid devices: CPU, FPGA")
endif()

# Define a list of supported devices
set(VALID_DEVICES "CPU" "FPGA")

# Counter for the number of valid devices
set(VALID_DEVICES_COUNT 0)

foreach(DEV IN LISTS DEVICE)
    # Loop through the list of devices and check if one is enabled
    foreach(VALID_DEVICE ${VALID_DEVICES})

        # Variable to hold which of the valid devices is the device
        set(SELECTED_DEVICE "")

        if(DEV STREQUAL VALID_DEVICE)
            set(SELECTED_DEVICE ${VALID_DEVICE})
            message(STATUS "Device ${DEV} is valid.")
            math(EXPR VALID_DEVICES_COUNT "${VALID_DEVICES_COUNT} + 1")
            message(STATUS "Configuring device ${DEV}...")
            add_subdirectory(${DEV})
            message(STATUS "Device ${DEV} correctly configured.")
            break()
        endif()

    endforeach()

    # If device is not valid
    if(SELECTED_DEVICE STREQUAL "")
        message(STATUS "Unknown device: ${DEV}. Valid devices: ${VALID_DEVICES}")
    endif()

endforeach()

# If there isn't a single valid device
if(VALID_DEVICES_COUNT EQUAL 0)
    message(FATAL_ERROR "There isn't any valid device in devices ${DEVICE}. Valid devices: ${VALID_DEVICES}")
endif()