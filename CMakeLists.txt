# how to write CMakeList.txt: https://subingwen.cn/cmake/CMake-prime
cmake_minimum_required (VERSION 3.25)

MESSAGE(STATUS "Operation system is ${CMAKE_SYSTEM}")

# Project information
project (BaseVar VERSION 1.1.0)
set(CMAKE_CXX_STANDARD 17)  # C++ 17

# 指定可执行程序的输出路径
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR}/htslib)      # 包含 htslib 的头文件
# link_directories(${CMAKE_SOURCE_DIR}/htslib/)      # 链接编译好的 htslib 静态链接库路径，是个目录
link_libraries(${PROJECT_SOURCE_DIR}/htslib/libhts.a)  # 链接编译好的 htslib 静态链接库 htslib.a，是个文件

# 查找当前目录下的所有源文件, 并将名称保存到 DIR_SRCS 变量
# aux_source_directory(src DIR_SRCS)
# file(GLOB main_head ${CMAKE_SOURCE_DIR}/src/*.h)  # 我注意到这个已经不需要添加了
file(GLOB SRC_LIST ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)

# 生成可执行程序
add_executable(basevar ${SRC_LIST})
IF (APPLE)
    target_link_libraries(basevar ${COMMON_SYS_LIB})
ELSEIF (LINUX)
    target_link_libraries(basevar ${COMMON_SYS_LIB} ssl crypto)
ELSE ()
    MESSAGE(STATUS "Only support for Linux and MacOS operator and not for ${CMAKE_SYSTEM_NAME}")
ENDIF ()
 