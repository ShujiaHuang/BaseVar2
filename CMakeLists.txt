cmake_minimum_required (VERSION 3.12)

# Project information
project (BaseVar VERSION 1.1.0)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/htslib)
link_directories(${CMAKE_SOURCE_DIR}/htslib)

# 查找当前目录下的所有源文件, 并将名称保存到 DIR_SRCS 变量
aux_source_directory(src DIR_SRCS)

# 生成可执行程序
add_executable(basevar ${DIR_SRCS})



