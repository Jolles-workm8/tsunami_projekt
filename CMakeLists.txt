cmake_minimum_required(VERSION 3.1)

project(tsunami_lab LANGUAGES CXX)

find_package(Catch2 QUIET)
# this time, just to demonstrate, we're fetching Catch2 in case it is not
# already available on the system
if (NOT Catch2_FOUND)
  include(FetchContent)

  FetchContent_Declare(
      catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG        v2.13.3
  )

  FetchContent_MakeAvailable(catch2)
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
endif()

include(CTest)
include(Catch)
include_directories(./src)

add_executable(tsunami_lab main.cpp)
