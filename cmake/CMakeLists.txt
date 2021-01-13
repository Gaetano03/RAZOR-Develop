#project(utopiae-orbit-dependencies
#        LANGUAGES CXX
#        VERSION 1.0.0)

# Include external libraries
include(ExternalProject)

set(EIGEN_INCLUDE_DIR /usr/local/Cellar/eigen/3.3.7/include/eigen3/Eigen)
include_directories (${EIGEN_INCLUDE_DIR})

# When gsl is required
if (ENABLE_GSL)
    find_package(GSL REQUIRED)
    list(APPEND MANDATORY_LIBRARIES "${GSL_LIBRARIES}")
    include_directories(${GSL_INCLUDE_DIRS})
endif()

# Smart math is required
include(cmake/external/smart-math.cmake)
list (APPEND MANDATORY_LIBRARIES "${SMART_MATH_STATIC_LIBRARY}")
list (APPEND MANDATORY_LIBRARIES "${SMART_MATH_LIBRARY}")
include_directories (${SMART_MATH_INCLUDE_DIR})

# Smart uq is optional
if (ENABLE_UQ)

    # Add compiler flag
    add_definitions(-D__USE_SMARTUQ)

    include(cmake/external/smart-uq.cmake)
    list (APPEND MANDATORY_LIBRARIES "${SMART_UQ_STATIC_LIBRARY}")
    list (APPEND MANDATORY_LIBRARIES "${SMART_UQ_LIBRARY}")
    include_directories (${SMART_UQ_INCLUDE_DIR})
endif()

# When HDF5 is required
if (ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C CXX REQUIRED)
    list(APPEND MANDATORY_LIBRARIES "${HDF5_CXX_LIBRARIES}")
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

# Pagmo is optional
if (ENABLE_PAGMO)

    # Add compiler flag
    add_definitions(-D__USE_PAGMO)

    include(cmake/external/pagmo.cmake)
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_STATIC_LIBRARY}")
    list (APPEND MANDATORY_LIBRARIES "${PAGMO_LIBRARY}")

    include_directories (${PAGMO_INCLUDE_DIR})

endif()


