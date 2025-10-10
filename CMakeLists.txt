cmake_minimum_required(VERSION 3.20)

project(gulls LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_WARNINGS "Enable compiler warnings" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(GSL REQUIRED)
find_package(CFITSIO REQUIRED)

set(SRC_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(CLASSES_DIR  "${SRC_DIR}/classes")
set(HEADERS_DIR  "${SRC_DIR}/headers")

set(COMMON_SOURCES
  "${SRC_DIR}/gulls.cpp"
  "${SRC_DIR}/strfns.cpp"
  "${SRC_DIR}/readParamfile.cpp"
  "${SRC_DIR}/readObservatoryfile.cpp"
  "${SRC_DIR}/buildWorld.cpp"
  "${CLASSES_DIR}/random.cpp"
  "${SRC_DIR}/info.cpp"
  "${SRC_DIR}/timeSequencer.cpp"
  "${SRC_DIR}/astroFns.cpp"
  "${SRC_DIR}/mathFns.cpp"
  "${CLASSES_DIR}/wittFSPL.cpp"
  "${CLASSES_DIR}/image.cpp"
  "${CLASSES_DIR}/psf.cpp"
  "${SRC_DIR}/readStarfields.cpp"
  "${SRC_DIR}/readSLList.cpp"
  "${SRC_DIR}/readPlanets.cpp"
  "${CLASSES_DIR}/parallax.cpp"
  "${CLASSES_DIR}/ephem.cpp"
  "${CLASSES_DIR}/coords.cpp"
  "${CLASSES_DIR}/zodiacalLight.cpp"
  "${CLASSES_DIR}/rbf.cpp"
  "${SRC_DIR}/photometry.cpp"
  "${SRC_DIR}/buildEvent.cpp"
  "${CLASSES_DIR}/VBMicrolensingLibrary.cpp"
)

set(CURVEGEN_FORTRAN
  "${SRC_DIR}/track.f"
  "${SRC_DIR}/coeff.f"
  "${SRC_DIR}/zroots.f"
  "${SRC_DIR}/findImages.f"
  "${SRC_DIR}/inside.f"
  "${SRC_DIR}/indexx.f"
  "${SRC_DIR}/limb.f"
  "${SRC_DIR}/uniform.f"
  "${SRC_DIR}/gould.f"
  "${SRC_DIR}/area.f"
  "${SRC_DIR}/mag.f"
)

set(FSBLCOMMON_SOURCES
  "${SRC_DIR}/pllxLightcurveGenerator.cpp"
  "${SRC_DIR}/outputLightcurve.cpp"
  "${SRC_DIR}/backupGenerator.cpp"
  "${CLASSES_DIR}/fs.cpp"
  "${CLASSES_DIR}/zroots2.cpp"
  "${CLASSES_DIR}/binaryMag.cpp"
  "${SRC_DIR}/pllxLightcurveFitter.cpp"
  "${SRC_DIR}/pllxLightcurveFitter_FS.cpp"
  ${CURVEGEN_FORTRAN}
)

add_library(gulls_core STATIC
  ${COMMON_SOURCES}
  ${FSBLCOMMON_SOURCES}
)

target_include_directories(gulls_core
  PUBLIC
    "${SRC_DIR}"
    "${HEADERS_DIR}"
    "${CLASSES_DIR}"
)

target_link_libraries(gulls_core
  PUBLIC
    CFITSIO::CFITSIO
    GSL::gsl
    GSL::gslcblas
)

if(ENABLE_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(gulls_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-deprecated-declarations>)
  elseif(MSVC)
    target_compile_options(gulls_core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4>)
  endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

function(add_gulls_executable target_name)
  add_executable(${target_name} ${ARGN})
  target_link_libraries(${target_name}
    PRIVATE
      gulls_core
  )
  if(UNIX)
    target_link_libraries(${target_name} PRIVATE m)
    if(NOT APPLE)
      find_library(DL_LIBRARY dl)
      if(DL_LIBRARY)
        target_link_libraries(${target_name} PRIVATE ${DL_LIBRARY})
      endif()
    endif()
  endif()
endfunction()

add_gulls_executable(gulls_std.x
  "${SRC_DIR}/detectionCuts.cpp"
  "${SRC_DIR}/standardPlanet.cpp"
)

add_gulls_executable(gulls_croin.x
  "${SRC_DIR}/detectionCuts.cpp"
  "${SRC_DIR}/croinPlanet.cpp"
)

add_gulls_executable(gullsFish.x
  "${SRC_DIR}/fisherDetCuts.cpp"
  "${SRC_DIR}/fisher.cpp"
  "${SRC_DIR}/croinPlanet.cpp"
  "${CLASSES_DIR}/fisherinversion.cpp"
)