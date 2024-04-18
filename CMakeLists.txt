cmake_minimum_required(VERSION 3.11)
PROJECT(DGtalTest)

#Required in DGtal
set(CMAKE_CXX_STANDARD 11)
add_subdirectory("/usr/local/polyscope" "polyscope")
set(SRCS
        main.cpp
)

add_executable(project "${SRCS}")
target_link_libraries(project polyscope)

FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

TARGET_LINK_LIBRARIES(project ${DGTAL_LIBRARIES})