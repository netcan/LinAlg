project(Matrix)
set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES Matrix.h Matrix.cpp)
add_library(matrix ${SOURCE_FILES})
add_executable(main main.cpp)
target_link_libraries(main matrix)
