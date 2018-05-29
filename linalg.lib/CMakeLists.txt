set(CMAKE_CXX_STANDARD 14)

set(LIB_SOURCE_FILES src/LinAlg.cpp)
set(LIB_HEADER_FILES include/LinAlg.h)

include_directories(include)
add_library(linalg SHARED ${LIB_SOURCE_FILES})
set_target_properties(linalg PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 0
        PUBLIC_HEADER ${LIB_HEADER_FILES})

configure_file(linalg.pc.in linalg.pc @ONLY)
include(GNUInstallDirs)
install(TARGETS linalg
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/linalg.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

add_executable(main src/main.cpp)
target_link_libraries(main linalg)
