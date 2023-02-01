cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(PeterProcessor)
set(CMAKE_CXX_STANDARD 17)


find_package(ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES)
include(ilcsoft_default_settings)

# Remove this ugly warning
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-effc++")
message(STATUS "C++ flags are set to: " ${CMAKE_CXX_FLAGS})


include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCE_DIR}/src/PeterProcessor.cpp)

### DEPENDENCIES ###
find_package(Marlin REQUIRED)
find_package(MarlinUtil REQUIRED)

FOREACH(pkg Marlin MarlinUtil)
   include_directories(${${pkg}_INCLUDE_DIRS})
   target_link_libraries(${PROJECT_NAME} ${${pkg}_LIBRARIES})
ENDFOREACH()

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)
