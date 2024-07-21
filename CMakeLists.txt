cmake_minimum_required(VERSION 3.10)
project(MyProject)

# 设置编译选项
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src/3rdparty/mdspan)

add_executable(App src/main.cpp)
target_link_libraries(App PUBLIC mdspan)
target_include_directories(App PUBLIC src/3rdparty/mdspan/include)