cmake_minimum_required(VERSION 3.25)
project(PMTCalib)

#----------------------------------------------------------------------------
# Find dependencies
find_package(ROOT REQUIRED)
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/*.h)
#remove LinkDef.h from headers list
list(FILTER HEADERS EXCLUDE REGEX "LinkDef.h$")
install(FILES ${HEADERS} DESTINATION include)

#----------------------------------------------------------------------------
# Generate ROOT dictionary for use in PyROOT notebooks.
ROOT_GENERATE_DICTIONARY(G__PMTCalib
    ${HEADERS}
    LINKDEF ${PROJECT_SOURCE_DIR}/LinkDef.h
    MODULE PMTCalib
    OPTIONS
    -I${PROJECT_SOURCE_DIR}/src
)

install(FILES
    $<TARGET_FILE_DIR:PMT_CALIB_LIB>/libPMTCalib_rdict.pcm
    $<TARGET_FILE_DIR:PMT_CALIB_LIB>/libPMTCalib.rootmap
    DESTINATION lib
)

#----------------------------------------------------------------------------
# Create a shared library from the dictionary
add_library(PMT_CALIB_LIB SHARED ${SOURCES} G__PMTCalib.cxx)
target_include_directories(PMT_CALIB_LIB PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
)
target_link_libraries(PMT_CALIB_LIB PUBLIC
    ROOT::Core ROOT::Hist ROOT::MathMore
    ${ROOT_COMPONENT_LIBRARIES}
    PkgConfig::FFTW3
)
install(TARGETS PMT_CALIB_LIB LIBRARY DESTINATION lib)

#----------------------------------------------------------------------------
# Export include dir variable for downstream use (optional)
set(PMT_CALIB_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "PMTCalib include dir")