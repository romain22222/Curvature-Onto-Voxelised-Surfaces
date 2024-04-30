cmake_minimum_required(VERSION 3.11)
PROJECT(DGtalTest)

#Required in DGtal
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/recipes)

include(FetchContent)

FetchContent_Declare(polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
)

FetchContent_MakeAvailable(polyscope)

set (SRCS
        main.cpp
)


add_executable(project "${SRCS}")
target_link_libraries(project polyscope)

FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

TARGET_LINK_LIBRARIES(project ${DGTAL_LIBRARIES})
