cmake_minimum_required(VERSION 3.14)
project(pqc-test C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O2")

include(FetchContent)

# Fetch canokey-crypto for SHA3 implementation
FetchContent_Declare(
    canokey_crypto
    GIT_REPOSITORY https://github.com/canokeys/canokey-crypto.git
    GIT_TAG master
)
FetchContent_GetProperties(canokey_crypto)
if(NOT canokey_crypto_POPULATED)
    FetchContent_Populate(canokey_crypto)
endif()

# Build only sha3.c (+ memzero.c dependency) from canokey-crypto
add_library(sha3 STATIC
    ${canokey_crypto_SOURCE_DIR}/src/sha3.c
    ${canokey_crypto_SOURCE_DIR}/src/memzero.c
)
target_include_directories(sha3 PUBLIC ${canokey_crypto_SOURCE_DIR}/include)

# ML-DSA-65 library
add_library(ml-dsa-65 STATIC ml-dsa-65.c se_impl.c)
target_include_directories(ml-dsa-65 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ml-dsa-65 PUBLIC sha3)

# Test executable
add_executable(test_ml_dsa test_main.c)
target_link_libraries(test_ml_dsa ml-dsa-65)

# Debug test
add_executable(debug_test debug_test.c)
target_link_libraries(debug_test ml-dsa-65)
