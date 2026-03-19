cmake_minimum_required(VERSION 3.14)
project(pqc-test C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O2 -ffunction-sections -fdata-sections")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,--gc-sections")

include(FetchContent)

# Fetch canokey-crypto for SHA3 implementation
FetchContent_Declare(
    canokey_crypto
    GIT_REPOSITORY https://github.com/canokeys/canokey-crypto.git
    GIT_TAG master
    SOURCE_SUBDIR DO_NOT_RUN_CMAKE
)
# Fetch mldsa-native reference implementation for cross-check testing
FetchContent_Declare(
    mldsa_native
    GIT_REPOSITORY https://github.com/pq-code-package/mldsa-native.git
    GIT_TAG v1.0.0-alpha
    SOURCE_SUBDIR DO_NOT_RUN_CMAKE
)
FetchContent_MakeAvailable(canokey_crypto mldsa_native)

# Build only sha3.c (+ memzero.c dependency) from canokey-crypto
add_library(sha3 STATIC
    ${canokey_crypto_SOURCE_DIR}/src/sha3.c
    ${canokey_crypto_SOURCE_DIR}/src/memzero.c
)
target_include_directories(sha3 PUBLIC ${canokey_crypto_SOURCE_DIR}/include)

add_library(mldsa-native-65 STATIC
    ${mldsa_native_SOURCE_DIR}/mldsa/mldsa_native.c
)
target_include_directories(mldsa-native-65 PUBLIC
    ${mldsa_native_SOURCE_DIR}/mldsa
)
target_compile_definitions(mldsa-native-65 PRIVATE
    MLD_CONFIG_PARAMETER_SET=65
    MLD_CONFIG_NAMESPACE_PREFIX=mldsa
)

# Find mbedtls (required by se_impl.c)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(MBEDTLS mbedtls)
endif()
if(MBEDTLS_FOUND)
  set(MBEDTLS_LIB ${MBEDTLS_STATIC_LIBRARIES})
  message(STATUS "Using mbedtls from pkg-config: ${MBEDTLS_LIB}")
else()
  find_path(MBEDTLS_INCLUDE_DIRS mbedtls/bignum.h
    HINTS /opt/homebrew/include /opt/homebrew/anaconda3/include /usr/local/include /usr/include)
  find_library(MBEDTLS_CRYPTO_LIB mbedcrypto
    HINTS /opt/homebrew/lib /opt/homebrew/anaconda3/lib /usr/local/lib /usr/lib)
  if(MBEDTLS_INCLUDE_DIRS AND MBEDTLS_CRYPTO_LIB)
    message(STATUS "Using mbedtls from search path: ${MBEDTLS_CRYPTO_LIB}")
    set(MBEDTLS_FOUND TRUE)
    set(MBEDTLS_LIB ${MBEDTLS_CRYPTO_LIB})
  else()
    message(FATAL "mbedtls not found!")
  endif()
endif()

# ML-DSA-65 library
add_library(ml-dsa-65 STATIC ml-dsa-65.c se_impl.c)
target_include_directories(ml-dsa-65 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ml-dsa-65 PUBLIC sha3)
target_include_directories(ml-dsa-65 PRIVATE ${MBEDTLS_INCLUDE_DIRS})
target_link_libraries(ml-dsa-65 PUBLIC ${MBEDTLS_LIB})
# Set RPATH so runtime can find the shared library
# get_filename_component(MBEDTLS_LIB_DIR "${MBEDTLS_CRYPTO_LIB}" DIRECTORY)
# set(CMAKE_INSTALL_RPATH "${MBEDTLS_LIB_DIR}")
# set(CMAKE_BUILD_RPATH "${MBEDTLS_LIB_DIR}")

# Test executable
add_executable(test_ml_dsa
    test_main.c
    ${mldsa_native_SOURCE_DIR}/test/notrandombytes/notrandombytes.c
)
target_include_directories(test_ml_dsa PRIVATE
    ${canokey_crypto_SOURCE_DIR}/include
    ${mldsa_native_SOURCE_DIR}/mldsa
)
target_link_libraries(test_ml_dsa ml-dsa-65 mldsa-native-65)

# Debug test
add_executable(debug_test debug_test.c)
target_link_libraries(debug_test ml-dsa-65)

