cmake_minimum_required(VERSION 3.2)

project(BOSS_MODULE LANGUAGES Fortran C CXX VERSION 0.1)
    if((DEFINED ENV{FC}) AND ("$ENV{FC}" MATCHES "^.*ifort$"))
        set(FC $ENV{FC} CACHE STRING "Fortran compiler")
        message("FC set to ${FC}")
    else()
        set(FC "gfortran" CACHE STRING "Fortran compiler")
        message("Defaulting to gfortran")
    endif()

    set(QUIET_MAKE "--no-print-directory")

    if(APPLE)
        set(LAPACK_LINK "-framework Accelerate")
    else()
        if(FC MATCHES "^.*ifort$")
            set(BLA_VENDOR Intel10_64lp)
        endif()
        find_package(LAPACK REQUIRED)
        list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
        get_filename_component(LAPACK_LIB_DIR ${LAPACK_LIB} DIRECTORY)
        set(LAPACK_LINK "-L${LAPACK_LIB_DIR} -llapack -lblas")
        message(STATUS "LAPACK_LINK: ${LAPACK_LIB} ${LAPACK_LIB_DIR}")
    endif()

    add_custom_target(cuba1.4
                      COMMENT "Running make lib for cuba1.4"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} CFLAGS=-fPIC ./configure
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} make ${QUIET_MAKE} lib
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cuba-1.4/
                      )
    add_custom_target(cuba1.4-clean
                      COMMENT "Running make distclean for cuba1.4"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} make ${QUIET_MAKE} distclean
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cuba-1.4/
                      )

    add_custom_target(boss_module
                      COMMENT "Running make for boss module"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} LAPACK_LINK=${LAPACK_LINK} make ${QUIET_MAKE}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_interface
                      DEPENDS cuba1.4)
    add_custom_target(boss_module-clean
                      COMMENT "Running make clean for  boss module"
                      COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${FC} LAPACK_LINK=${LAPACK_LINK} make ${QUIET_MAKE} clean
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_interface
                      DEPENDS cuba1.4-clean)
