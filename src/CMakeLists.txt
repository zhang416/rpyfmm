FILE(GLOB sources *.cc)

include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${dashmm_SOURCE_DIR}/include
)

link_directories(${HPX_LIBRARY_DIRS})

add_library(rpyfmm STATIC ${sources})
target_link_libraries(rpyfmm ${HPX_LDFLAGS})

install(TARGETS rpyfmm DESTINATION lib)
