cmake_minimum_required(VERSION 3.2)
project(sword LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(sword_optimize_for_portability "Build sword with -msse4.1" OFF)

if (sword_optimize_for_portability)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

add_executable(sword
    src/chain.cpp
    src/database_alignment.cpp
    src/database_search.cpp
    src/evalue.cpp
    src/hash.cpp
    src/kmers.cpp
    src/main.cpp
    src/reader.cpp
    src/score_matrix.cpp
    src/utils.cpp
    src/writer.cpp
    vendor/opal/src/opal.cpp)

include_directories(vendor/opal/src)

if (NOT TARGET thread_pool)
    add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
endif()

target_link_libraries(sword thread_pool)

install(TARGETS sword DESTINATION ${CMAKE_INSTALL_BINDIR})
