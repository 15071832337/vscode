# stop if cmake version is below 3.0
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
#set(nlohmann_json_DIR "/opt/homebrew/opt/nlohmann-json/lib/cmake/nlohmann_json")
#find_package(nlohmann_json_DIR 3.11.3 REQUIRED)
# project name and enable C++ support
project(rcbkdipole CXX C)

# Optimization on (for debug flags, change to Debug)
set(CMAKE_BUILD_TYPE Release)

find_package(GSL REQUIRED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
)

add_subdirectory(src)

set_target_properties(amplitude PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
	ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

##########
enable_language(Fortran)
add_library(KKP_Fortran STATIC
    src/fragmentation_kkp.f
)
# 设置C++17标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 显式指定ROOT路径
set(ROOT_DIR "/Users/suyoupeng/alice/sw/osx_arm64/ROOT/v6-32-06-alice7-local1")

#set(ROOT_DIR $ENV{ROOTSYS})


# 查找依赖包
find_package(ROOT REQUIRED PATHS ${ROOT_DIR} NO_DEFAULT_PATH)
find_package(GSL REQUIRED)

# 打印调试信息
#message(STATUS "ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")
#message(STATUS "ROOT_LIBRARIES = ${ROOT_LIBRARIES}")
#message(STATUS "GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIRS}")
#message(STATUS "GSL_LIBRARIES = ${GSL_LIBRARIES}")

# 先构建基础库
add_library(Interpolator STATIC src/interpolation.cpp)
add_library(DataFile STATIC src/datafile.cpp)
add_library(Tools STATIC src/tools.cpp)
add_library(KKP STATIC src/kkp.cpp)
add_library(FRAG STATIC src/fragmentation.cpp)
add_library(GBW STATIC src/GBW.cpp)

# 设置基础库的包含路径
target_include_directories(Interpolator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(DataFile PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
#target_include_directories(GBW PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# 构建AmplitudeLib库
add_library(AmplitudeLib STATIC src/amplitudelib.cpp)

# 设置AmplitudeLib的属性和链接
target_include_directories(AmplitudeLib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(AmplitudeLib PRIVATE
    Interpolator
    DataFile
	Tools
    ${ROOT_LIBRARIES}
)

# 构建主程序
add_executable(CGCTest src/CGCTest.cpp)


# 设置主程序的包含路径
target_include_directories(CGCTest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${ROOT_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
)



# 链接所有必需的库
target_link_libraries(CGCTest
    AmplitudeLib
    Interpolator
    DataFile
	Tools
	FRAG
	KKP        
    KKP_Fortran   
    ${ROOT_LIBRARIES}
    ${GSL_LIBRARIES}
    GBW
)





