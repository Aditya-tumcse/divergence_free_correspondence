find_package(Catch2 REQUIRED)

add_executable(unit_tests
    test_utilities.cpp
    test_numerics.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        divergence_free_corr
        Catch2::Catch2WithMain
        project_options
)

include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
