cmake_minimum_required(VERSION 3.5)

project(FirstCurve LANGUAGES C)

find_package(PkgConfig REQUIRED)
pkg_check_modules(gmp REQUIRED IMPORTED_TARGET gmp)

add_executable(FirstCurve main.c curve.c curve.h)


target_link_libraries(FirstCurve PRIVATE PkgConfig::gmp)
