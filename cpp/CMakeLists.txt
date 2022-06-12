cmake_minimum_required(VERSION 3.14)
project(HaloPathMath VERSION 2.1.1)

set(CMAKE_VERBOSE_MAKEFILE OFF)
option(DEBUG "Default OFF" OFF)
option(BUILD_TEST "Default ON" ON)

# ==================================================================================================
# Compiler related
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  # using regular Clang or AppleClang
  set(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL} -dead_strip")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL} --gc-sections")
endif()

if(MSVC)
  add_definitions(-DNOMINMAX)
  set(CMAKE_DEBUG_POSTFIX "d")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fpermissive -pthread -frtti")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -D__STDC_FORMAT_MACROS -ftime-trace")
endif()

# ==================================================================================================
# OS related

# ==================================================================================================
# Basic
set(CMAKE_CXX_STANDARD 17)            # C++17...
set(CMAKE_CXX_STANDARD_REQUIRED ON)   #...is required...
set(CMAKE_CXX_EXTENSIONS OFF)         #...without compiler extensions like gnu++11
set(CMAKE_C_STANDARD 11)              # C11...

# Flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -ffunction-sections -fdata-sections -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffunction-sections -fdata-sections -fvisibility=hidden")

include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result AND "${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "IPO supported. Enable inter-procedural optimization.")
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Directories
set(BUILD_DIR "${CMAKE_SOURCE_DIR}/build")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BUILD_DIR}/${CMAKE_BUILD_TYPE}/lib/${ARCH}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BUILD_DIR}/${CMAKE_BUILD_TYPE}/lib/${ARCH}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BUILD_DIR}/${CMAKE_BUILD_TYPE}/bin/${ARCH}")

set(PROJ_SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(PROJ_TEST_DIR "${CMAKE_SOURCE_DIR}/test")
set(PROJ_ROOT "${CMAKE_SOURCE_DIR}")

# ==================================================================================================
# Eigen3
find_package(Eigen3 3.3 REQUIRED)

# ==================================================================================================
# Add some definitions
if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux")
  add_compile_definitions(OS_LINUX) 
elseif(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")
  add_compile_definitions(OS_MAC)
elseif(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
  add_compile_definitions(OS_WIN)
endif()


# ==================================================================================================
# Set submodules
set(proj_src
	"${PROJ_SRC_DIR}/core/geo.cpp"
	"${PROJ_SRC_DIR}/sph_grid.cpp"
	"${PROJ_SRC_DIR}/util/log.cpp")

set(proj_include "${PROJ_SRC_DIR}")

# ==================================================================================================
# Add obj-lib and corresponding target definitions
add_library(proj_lib OBJECT ${proj_src})
target_include_directories(proj_lib PUBLIC
  "${proj_include}")
target_link_libraries(proj_lib PUBLIC Eigen3::Eigen)

# ==================================================================================================
# Add subdirectory
add_subdirectory(src)
if(BUILD_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

