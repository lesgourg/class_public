cmake_minimum_required(VERSION 3.16)
project(class_public_tt LANGUAGES C CXX)

# ---------- Build setup ----------
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE) # Debug/Release/RelWithDebInfo/MinSizeRel
endif()

# Safe, per-language flags (fixed: each flag is its own generator-expression)
# C flags
add_compile_options(
        $<$<COMPILE_LANGUAGE:C>:-fPIC>
        $<$<COMPILE_LANGUAGE:C>:-g>
)
# C++ flags
add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:-fPIC>
        $<$<COMPILE_LANGUAGE:CXX>:-g>
        $<$<COMPILE_LANGUAGE:CXX>:-fpermissive>
        $<$<COMPILE_LANGUAGE:CXX>:-Wno-write-strings>
)

# Optimize only in Release (mirrors OPTFLAG = -O3)
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
endif()

# Pthreads (your OMPFLAG was -pthread)
find_package(Threads REQUIRED)

# Optional OpenMP (off by default; enable if you used -fopenmp)
option(WITH_OPENMP "Enable OpenMP" OFF)
if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# Define __CLASSDIR__ like in the Makefile
add_compile_definitions(__CLASSDIR__="${CMAKE_SOURCE_DIR}")

# ---------- Include dirs ----------
include_directories(
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/external/RecfastCLASS
        ${CMAKE_SOURCE_DIR}/external/heating
        ${CMAKE_SOURCE_DIR}/external/Halofit
        ${CMAKE_SOURCE_DIR}/external/HMcode
)

# HyRec toggle (matches HYREC behavior in your Makefile)
option(WITH_HYREC "Build with HyRec" ON)
set(HYREC_DIR ${CMAKE_SOURCE_DIR}/external/HyRec2020 CACHE PATH "HyRec path")
if(WITH_HYREC)
    add_compile_definitions(HYREC)
    include_directories(${HYREC_DIR})
endif()

# ---------- Helper to handle file location differences ----------
function(append_first_existing OUT_VAR)
    foreach(cand IN LISTS ARGN)
        if(EXISTS "${CMAKE_SOURCE_DIR}/${cand}")
            list(APPEND ${OUT_VAR} "${CMAKE_SOURCE_DIR}/${cand}")
            set(${OUT_VAR} "${${OUT_VAR}}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    message(FATAL_ERROR "None of these files exist: ${ARGN}")
endfunction()

# ---------- Sources ----------
# Files compiled as C++ in your Makefile (".opp" rule)
set(CXX_SOURCES "")
append_first_existing(CXX_SOURCES "tools/arrays.c")
append_first_existing(CXX_SOURCES "tools/hyperspherical.c")
append_first_existing(CXX_SOURCES "source/perturbations.c")
append_first_existing(CXX_SOURCES "source/primordial.c")
append_first_existing(CXX_SOURCES "source/transfer.c")
append_first_existing(CXX_SOURCES "source/harmonic.c")
append_first_existing(CXX_SOURCES "external/HMcode/hmcode.c")

# Pure C sources
set(C_SOURCES "")
# tools/
append_first_existing(C_SOURCES "tools/growTable.c")
append_first_existing(C_SOURCES "tools/dei_rkck.c")
append_first_existing(C_SOURCES "tools/sparse.c")
append_first_existing(C_SOURCES "tools/evolver_rkck.c")
append_first_existing(C_SOURCES "tools/evolver_ndf15.c")
append_first_existing(C_SOURCES "tools/parser.c")
append_first_existing(C_SOURCES "tools/quadrature.c")
append_first_existing(C_SOURCES "tools/common.c" "source/common.c")
append_first_existing(C_SOURCES "source/trigonometric_integrals.c" "tools/trigonometric_integrals.c")

# source/
append_first_existing(C_SOURCES "source/input.c")
append_first_existing(C_SOURCES "source/background.c")
append_first_existing(C_SOURCES "source/thermodynamics.c")
append_first_existing(C_SOURCES "source/fourier.c")
append_first_existing(C_SOURCES "source/lensing.c")
append_first_existing(C_SOURCES "source/distortions.c")
append_first_existing(C_SOURCES "source/output.c")
append_first_existing(C_SOURCES "main/class.c")

# externals (always C in your Makefile)
append_first_existing(C_SOURCES "external/RecfastCLASS/wrap_recfast.c")
append_first_existing(C_SOURCES "external/heating/injection.c")
append_first_existing(C_SOURCES "external/heating/noninjection.c")
append_first_existing(C_SOURCES "external/Halofit/halofit.c")

# HyRec (optional)
if(WITH_HYREC)
    append_first_existing(C_SOURCES "external/HyRec2020/hyrectools.c")
    append_first_existing(C_SOURCES "external/HyRec2020/helium.c")
    append_first_existing(C_SOURCES "external/HyRec2020/hydrogen.c")
    append_first_existing(C_SOURCES "external/HyRec2020/history.c")
    append_first_existing(C_SOURCES "external/HyRec2020/wrap_hyrec.c")
    append_first_existing(C_SOURCES "external/HyRec2020/energy_injection.c")
endif()

# Treat the ".opp" equivalents as C++
foreach(src IN LISTS CXX_SOURCES)
    set_source_files_properties(${src} PROPERTIES LANGUAGE CXX)
endforeach()

# ---------- Targets ----------
# Static library libclass.a (exclude main/class.c from the lib)
set(LIB_SOURCES ${C_SOURCES} ${CXX_SOURCES})
list(REMOVE_ITEM LIB_SOURCES "${CMAKE_SOURCE_DIR}/main/class.c")

add_library(class STATIC ${LIB_SOURCES})
target_link_libraries(class PRIVATE Threads::Threads m)
if(WITH_OPENMP)
    target_link_libraries(class PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
endif()

# Executable "class" linking the library
add_executable(class_cli "${CMAKE_SOURCE_DIR}/main/class.c")
target_link_libraries(class_cli PRIVATE class)
set_target_properties(class_cli PROPERTIES OUTPUT_NAME "class")

# ---------- Optional tests (uncomment if you use them) ----------
# add_executable(test_loops     test/test_loops.c)              ; target_link_libraries(test_loops PRIVATE class)
# add_executable(test_loops_omp test/test_loops_omp.c)          ; target_link_libraries(test_loops_omp PRIVATE class)
# add_executable(test_harmonic  test/test_harmonic.c)           ; target_link_libraries(test_harmonic PRIVATE class)
# add_executable(test_transfer  test/test_transfer.c)           ; target_link_libraries(test_transfer PRIVATE class)
# add_executable(test_fourier   test/test_fourier.c)            ; target_link_libraries(test_fourier PRIVATE class)
# add_executable(test_perturbations test/test_perturbations.c)  ; target_link_libraries(test_perturbations PRIVATE class)
# add_executable(test_thermodynamics test/test_thermodynamics.c); target_link_libraries(test_thermodynamics PRIVATE class)
# add_executable(test_background    test/test_background.c)     ; target_link_libraries(test_background PRIVATE class)
# add_executable(test_hyperspherical test/test_hyperspherical.c)
# target_sources(test_hyperspherical PRIVATE test/test_hyperspherical.c tools/hyperspherical.c)
# set_source_files_properties(${CMAKE_SOURCE_DIR}/tools/hyperspherical.c PROPERTIES LANGUAGE C)
# target_link_libraries(test_hyperspherical PRIVATE m)

# ---------- Python 'classy' target (pip install) ----------
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_custom_target(classy
            COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER}
            ${Python3_EXECUTABLE} -m pip install .
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python
            USES_TERMINAL
            COMMENT "Installing Python package 'classy' with pip"
    )
endif()

# ---------- (Optional) install rules ----------
# install(TARGETS class_cli RUNTIME DESTINATION bin)
# install(TARGETS class     ARCHIVE DESTINATION lib)
