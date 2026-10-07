cmake_minimum_required(VERSION 3.16)

project(tgv-denoise LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(tgv-denoise
    src/tgv_denoise.cpp
)

target_include_directories(tgv-denoise
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)