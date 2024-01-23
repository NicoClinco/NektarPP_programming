cmake_minimum_required(VERSION 2.8)
PROJECT(Nektar-SESSION-READER)

# Nektar++ requires C++11. Try to infer this for older CMake versions (less than
# 3.1.0)
IF ("${CMAKE_VERSION}" VERSION_LESS "3.1")
    IF (NOT MSVC)
        INCLUDE(CheckCXXCompilerFlag)
        CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

        IF (COMPILER_SUPPORTS_CXX11)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        ELSE()
            MESSAGE(FATAL_ERROR "Nektar++ requires a compiler with C++11 support.")
        ENDIF()
    ELSEIF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        MESSAGE(FATAL_ERROR "Nektar++ requires full C++11 support: please upgrade to Visual Studio 2013 or later")
    ENDIF()
ELSE()
    SET(CMAKE_CXX_STANDARD 11)
    SET(CMAKE_CXX_STANDARD_REQUIRED ON)
    SET(CMAKE_CXX_EXTENSIONS OFF)
ENDIF()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/app)

# STANDARD BLOCK TO FIND THE LIBRARY #
FIND_PACKAGE(Nektar++ PATHS "/home/nico/nektar-v5.3.0/build")

INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${NEKTAR++_TP_INCLUDE_DIRS})
LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS})
LINK_DIRECTORIES(${NEKTAR++_TP_LIBRARY_DIRS})
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
#######################################

ADD_EXECUTABLE(ReadSession ReadXML.cpp)
TARGET_LINK_LIBRARIES(ReadSession ${NEKTAR++_LIBRARIES})

