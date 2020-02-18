cmake_minimum_required (VERSION 2.6)
project (icucockpit_libs)
# The version number.
set (ICUCOCKPIT_VERSION_MAJOR 1)
set (ICUCOCKPIT_VERSION_MINOR 0)

include_directories(include)
find_package(LibXml2 REQUIRED)
if (LibXml2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
endif (LibXml2_FOUND)

find_package(Boost COMPONENTS date_time REQUIRED)
if (Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif (Boost_FOUND)

include_directories(${icucockpit_utils_INCLUDE_DIRS})

file(GLOB SOURCES "src/*.cpp")
#file(GLOB INCLUDES "include/*.h")
add_library(${PROJECT_NAME} SHARED ${SOURCES}) 

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBRARIES} ${LibXml2_LIBRARIES} ${Boost_LIBRARIES})

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
