cmake_minimum_required(VERSION 3.25)
project(PMTCalib)

#----------------------------------------------------------------------------
# Find dependencies
find_package(ROOT REQUIRED COMPONENTS MathMore Minuit2)
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)
install(FILES ${HEADERS} DESTINATION include)

#----------------------------------------------------------------------------
# Generate ROOT dictionary for use in PyROOT notebooks.
add_library(PMTCalib SHARED ${SOURCES} ${G__PMTCalib})
target_include_directories(PMTCalib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
)
target_link_libraries(PMTCalib PUBLIC
    ROOT::Core ROOT::MathMore ROOT::Minuit2
    ${ROOT_COMPONENT_LIBRARIES}
    PkgConfig::FFTW3
)
set_target_properties(PMTCalib PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
)
   
ROOT_GENERATE_DICTIONARY(G__PMTCalib
    ${HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h
    MODULE PMTCalib
)

# ------------------------------------------------------------------------------
# Post-build copy (the only reliable way under Ninja Multi-Config)
# ROOT always places the dictionary here:
#   ${CMAKE_CURRENT_BINARY_DIR}/libPMTCalib_rdict.pcm
#   ${CMAKE_CURRENT_BINARY_DIR}/libPMTCalib.rootmap
#
# We mirror them next to the library:
#   $<TARGET_FILE_DIR:PMTCalib>/
#
add_custom_command(TARGET PMTCalib POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/libPMTCalib_rdict.pcm"
        "$<TARGET_FILE_DIR:PMTCalib>/libPMTCalib_rdict.pcm"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/libPMTCalib.rootmap"
        "$<TARGET_FILE_DIR:PMTCalib>/libPMTCalib.rootmap"
)

install(TARGETS PMTCalib LIBRARY DESTINATION lib)
install(FILES
    "$<TARGET_FILE_DIR:PMTCalib>/libPMTCalib_rdict.pcm"
    "$<TARGET_FILE_DIR:PMTCalib>/libPMTCalib.rootmap"
    DESTINATION lib
)

# ------------------------------------------------------------------------------
# Optionally compile examples
if(DEFINED PMTCalib_Examples)
    file(GLOB EXAMPLES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cc)
    foreach(example IN LISTS EXAMPLES)
        get_filename_component(EX_NAME "${example}" NAME_WE)
        add_executable(${EX_NAME} ${test} ${example})
        target_include_directories(${EX_NAME} PRIVATE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            ${ROOT_INCLUDE_DIRS}
        )
        target_link_libraries(${EX_NAME} PRIVATE
            PMTCalib
            ${ROOT_LIBRARIES}
            ${ROOT_COMPONENT_LIBRARIES}
        )
        install(TARGETS ${EX_NAME} RUNTIME DESTINATION bin)
    endforeach(example IN LIST EXAMPLES)
    
endif(DEFINED PMTCalib_Examples)

#----------------------------------------------------------------------------
# Export include dir variable for downstream use (optional)
set(PMT_CALIB_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "PMTCalib include dir")
