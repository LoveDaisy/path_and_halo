# ==================================================================================================
# Google Test
find_package(GTest)
if(NOT GTest_FOUND)
  if("${GTEST_DIR}" STREQUAL "")
    message(STATUS "[test] GTEST_DIR is NOT set! Use github source.")
    set(GTEST_REPO "https://github.com/google/googletest.git")
  else()
    message(STATUS "[test] GTEST_DIR is set! Use local gtest at ${GTEST_DIR}")
    set(GTEST_REPO "file://${GTEST_DIR}/.git")
  endif()
  
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY ${GTEST_REPO}
    GIT_TAG release-1.11.0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# ==================================================================================================
add_library(proj_test_lib OBJECT ${proj_src})
target_include_directories(proj_test_lib PUBLIC ${proj_include})
target_compile_definitions(proj_test_lib PUBLIC FOR_TEST)
target_compile_definitions(proj_test_lib PUBLIC VERBOSE_LOG)


add_executable(unit_test
  "${PROJ_TEST_DIR}/test_ad.cpp"
  "${PROJ_TEST_DIR}/test_geo.cpp"
  "${PROJ_TEST_DIR}/test_sph_grid.cpp")
target_include_directories(unit_test
  PRIVATE ${PROJ_TEST_DIR})
target_link_libraries(unit_test
  PRIVATE GTest::gtest GTest::gtest_main ${CMAKE_THREAD_LIBS_INIT} proj_test_lib)


add_test(NAME "UnitTest" COMMAND unit_test)

