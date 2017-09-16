add_executable(demo EXCLUDE_FROM_ALL demo.cc)

include_directories(
  ${HPX_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
  ${dashmm_SOURCE_DIR}/include
)

target_link_libraries(demo PUBLIC rpyfmm dashmm ${HPX_LDFLAGS})
