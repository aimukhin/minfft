cmake_minimum_required(VERSION 2.8)
project(MINFFT)

if (WIN32 AND MINGW)
  option(MINFFT_STATIC_BUILD "Build full static libraries on MinGW/Win32" ON)
else()
  set(MINFFT_STATIC_BUILD false)
endif()

option(DISABLE_LINK_WITH_M "Disables linking with m library to build with clangCL from MSVC" OFF)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(cmake/target_optimizations.cmake)
include(cmake/target_fast_math.cmake)
include(cmake/compiler_warnings.cmake)

# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "CMAKE_BUILD_TYPE not specified: defaulting to 'Release'.")
endif(NOT CMAKE_BUILD_TYPE)

######################################################

# Special MINGW stuff here
if(MINFFT_STATIC_BUILD)
  # see https://cmake.org/pipermail/cmake/2012-September/051970.html
  # see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc -static")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS} -static-libgcc -s")
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -static-libgcc -static-libstdc++ -s")
endif()


######################################################

# link with math library?
if(DISABLE_LINK_WITH_M)
  set(MATHLIB "")
else()
  if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    message(WARNING "INFO: detected MSVC: will not link math lib m")
    set(MATHLIB "")
  else()
    message(STATUS "INFO: detected NO MSVC: will link math lib m")
    set(MATHLIB "m")
  endif()
endif()

######################################################

set(MINFFT_VARIANTS float double extended)

foreach(X IN LISTS MINFFT_VARIANTS)
  add_library(MinFFT_${X} STATIC minfft.c minfft.h)

  # message(STATUS "set lib specs for ${X}")
  target_activate_c_compiler_warnings(MinFFT_${X})
  target_set_c_fast_math(MinFFT_${X})
  target_compile_options(MinFFT_${X} PRIVATE "${OPT_FAST_MATH}")

  target_include_directories(MinFFT_${X} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(MinFFT_${X} ${MATHLIB})
endforeach()

target_compile_definitions(MinFFT_float    PUBLIC MINFFT_SINGLE)
target_compile_definitions(MinFFT_extended PUBLIC MINFFT_EXTENDED)

######################################################

