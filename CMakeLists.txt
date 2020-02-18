cmake_minimum_required(VERSION 2.6)
project (ICUCockpitFeatureExtraction)

# You can tweak some common (for all subprojects) stuff here. For example:
if(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 ${SSE_FLAGS} -mavx") #-Wall -Wextra -fno-strict-aliasing -Wno-format-extra-args -Wno-unknown-pragmas 
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s") 
endif() 

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE   ON)

if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()

# All submodules have common dependencies on OpenCV. We find the packages here
find_package(OpenCV REQUIRED)

find_package(xml2 REQUIRED)
if (xml2_FOUND)
        include_directories(${XML2_INCLUDE_DIR})
endif (xml2_FOUND)


if (OpenCV_FOUND)
	include_directories(${OpenCV_INCLUDE_DIRS})
endif (OpenCV_FOUND)

find_package(dlib REQUIRED)
if (dlib_FOUND)
	include_directories(${DLIB_INCLUDE_DIR})
endif (dlib_FOUND)

# NOTE: The order matters! The most independent ones should go first.
add_subdirectory(modules/icucockpit_utils)
add_subdirectory(modules/icucockpit_libs) # icu_feature_extraction library
add_subdirectory(modules/icucockpit_service) 
