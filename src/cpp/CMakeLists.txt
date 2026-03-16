cmake_minimum_required(VERSION 3.14)
project(starfield_plugin CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(starfield STATIC
    src/catalog.cpp
    src/star_tracker.cpp
)
target_include_directories(starfield PUBLIC include)

add_executable(test_catalog tests/test_catalog.cpp)
target_link_libraries(test_catalog starfield)

add_executable(test_star_tracker tests/test_star_tracker.cpp)
target_link_libraries(test_star_tracker starfield)

enable_testing()
add_test(NAME catalog COMMAND test_catalog)
add_test(NAME star_tracker COMMAND test_star_tracker)

add_executable(test_sixdof tests/test_sixdof.cpp)
target_include_directories(test_sixdof PRIVATE include)
add_test(NAME starfield_sixdof COMMAND test_sixdof)
