cmake_minimum_required(VERSION 3.9)
project(LinAlg VERSION 0.1 DESCRIPTION "A Linear Algebra library for studying by netcan")
set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES Matrix.cpp)
set(HEADER_FILES Matrix.h)

add_library(linalg SHARED ${SOURCE_FILES})
set_target_properties(linalg PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 0
        PUBLIC_HEADER ${HEADER_FILES})

configure_file(linalg.pc.in linalg.pc @ONLY)
include(GNUInstallDirs)
install(TARGETS linalg
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/linalg.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)


add_executable(main main.cpp)
target_link_libraries(main linalg)
