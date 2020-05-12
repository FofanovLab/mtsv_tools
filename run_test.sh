#!/bin/bash -e
# build statically linked binary with Rust
C_INCLUDE_PATH=$PREFIX/include LIBRARY_PATH=$PREFIX/lib cargo test --manifest-path $PREFIX
