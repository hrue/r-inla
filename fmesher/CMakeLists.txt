## Usage:
##   1) Create a build directory and go there
##   2) Run "cmake /path/to/fmesher/"
##      Use "ccmake" instead of "cmake" to interactively set configuration variables
##   3) Run "make fmesher" and "make doc_fmesher"

cmake_minimum_required(VERSION 3.9)

option(FMESHER_USE_OPENMP "Use OpenMP" OFF)
option(FMESHER_BUILD_WITH_X "Build fmesher with X11 support" OFF)
option(FMESHER_BUILD_WITH_GSL "Build fmesher with GSL support" ON)
option(FMESHER_BUILD_FOR_DEBUGGING "Build fmesher with DEBUG activated" OFF)

if (FMESHER_USE_OPENMP)
  FIND_PACKAGE( OpenMP)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

if(FMESHER_BUILD_WITH_GSL)
  FIND_PACKAGE(GSL)
  if(GSL_FOUND)
    message("GSL FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GSL_EXE_LINKER_FLAGS}")
  else()
    message("GSL NOT FOUND; Spherical harmonics will not be available")
    ## Turn off the only use of gsl in fmesher:
    add_definitions(-DNO_SPHERICAL_HARMONICS=1)
  endif()
else()
  ## Turn off the only use of gsl in fmesher:
  add_definitions(-DNO_SPHERICAL_HARMONICS=1)
endif()

if(FMESHER_BUILD_FOR_DEBUGGING)
  message("DEBUG ACTIVATED")
  add_definitions(-DDEBUG)
endif()

SET_PROPERTY( SOURCE xtmpl.c x11utils.cc PROPERTY COMPILE_FLAGS -Wno-unused-result )
SET_PROPERTY( SOURCE xtmpl.c x11utils.cc PROPERTY COMPILE_FLAGS -Wno-pedantic )

if (CMAKE_SYSTEM_NAME MATCHES "^Linux$")
  add_definitions(-DLINUX=1)
endif()
if (CMAKE_SYSTEM_NAME MATCHES "^Darwin$")
  add_definitions(-DMAC=1)
endif()
if (CMAKE_SYSTEM_NAME MATCHES "^Windows$")
  add_definitions(-DWINDOWS=1)
endif()

if((FMESHER_BUILD_WITH_X) AND NOT (CMAKE_SYSTEM_NAME MATCHES "^Windows$"))
  message("Building with X11")
  set(FMESHER_X_FILES src/xtmpl.c src/x11utils.cc)
  set(FMESHER_LIBS ${FMESHER_LIBS} X11)
  add_definitions(-DFMESHER_USE_X)
else()
  message("Building without X11")
  ## Turn off using the X11 drawning feature:
  add_definitions(-DFMESHER_NO_X)
endif()

# Define executables #
add_library(fmesherlib STATIC
  src/predicates.cc
  src/mesh.cc
  src/meshc.cc
  src/vector.cc
  src/ioutils.cc
  src/cmdline.c
  src/basis.cc
  src/locator.cc
  src/trees.cc
  ${FMESHER_X_FILES})


## fmesher program for inla
add_executable(fmesher
  src/fmesher.cc)
add_executable(treetest
  src/treetest.cc)
target_link_libraries(fmesher fmesherlib ${FMESHER_LIBS} ${GSL_LIBRARIES})
target_link_libraries(treetest fmesherlib ${FMESHER_LIBS} ${GSL_LIBRARIES})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc_fmesher
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/doc/latex
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

#add_subdirectory(test)
