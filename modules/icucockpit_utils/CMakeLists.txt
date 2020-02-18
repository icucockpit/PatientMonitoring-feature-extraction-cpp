cmake_minimum_required (VERSION 2.6)
project (icucockpit_utils)
# The version number.
set (ICUCOCKPIT__VERSION_MAJOR 1)
set (ICUCOCKPIT__VERSION_MINOR 0)
 

include_directories(include)
file(GLOB SOURCES "src/*.cpp")
file(GLOB INCLUDES "include/*.cpp")

add_library(${PROJECT_NAME} SHARED ${SOURCES}) 
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBRARIES})

# Expose C's definitions (in this case only the ones of XXX transitively)
# to other subprojects through cache variable.
set(${PROJECT_NAME}_DEFINITIONS 
    CACHE INTERNAL "${PROJECT_NAME}: Definitions" FORCE)

# Expose C's public includes (including the ones of C's dependencies transitively)
# to other subprojects through cache variable.
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


install(TARGETS ${PROJECT_NAME}
	ARCHIVE DESTINATION /usr/local/lib
 	LIBRARY DESTINATION /usr/local/lib
 	COMPONENT library
)
install(DIRECTORY include/${PROJECT_NAME} DESTINATION /usr/local/include COMPONENT dev)
