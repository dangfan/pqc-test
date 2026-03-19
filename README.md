# pqc-test

A small C project for building and validating an `ML-DSA-65` (FIPS 204) implementation.

## Overview

This repository includes:

- `ml-dsa-65.c`, `ml-dsa-65.h`: core `ML-DSA-65` implementation (`keygen/sign/verify`)
- `se_impl.c`: software backend for `se.h` primitives (including Montgomery multiplication)
- `test_main.c`: regression-style functional tests
- `debug_test.c`: lightweight debug/smoke executable with progress output
- `CMakeLists.txt`: build configuration (fetches `sha3` from `canokey-crypto`)

## Requirements

- macOS or Linux
- CMake `>= 3.14`
- A C11-compatible compiler
- Git (required by CMake `FetchContent`)

## Build

From the project root:

```bash
cmake -S . -B build
cmake --build build -j
```

## Run Tests

### Regression test suite

```bash
./build/test_ml_dsa
```

`test_ml_dsa` runs multiple correctness checks, including:

- KeyGen + Sign + Verify round-trip
- Deterministic signing reproducibility
- Context-string positive/negative verification cases

### Debug / smoke run

```bash
./build/debug_test
```

`debug_test` is intended for quick developer feedback and prints step/timing progress such as:

- `KeyGen...`
- `Sign...`
- `Verify...`

## `test_ml_dsa` vs `debug_test`

- `test_ml_dsa`
  - Focus: correctness and regression checks
  - Contains explicit pass/fail validations across multiple scenarios
  - Best for routine validation before/after code changes

- `debug_test`
  - Focus: fast iteration and quick sanity checks
  - Provides concise progress/timing output
  - Useful during debugging sessions

## Project Layout

```text
.
├── CMakeLists.txt
├── ml-dsa-65.c
├── ml-dsa-65.h
├── se.h
├── se_impl.c
├── test_main.c
├── debug_test.c
└── README.md
```
