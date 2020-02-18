cmake_minimum_required(VERSION 2.6)

project(icucockpit_fe_service C CXX)

file(GLOB SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

include_directories(${icucockpit_libs_INCLUDE_DIRS})
include_directories(${icucockpit_utils_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} icucockpit_utils icucockpit_libs ${OpenCV_LIBRARIES} ${LibXml2_LIBRARIES} ${dlib_LIBRARIES})
