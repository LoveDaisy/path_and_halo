add_executable(render_halo
  "${PROJ_SRC_DIR}/render_main.cpp")
target_link_libraries(render_halo PRIVATE proj_lib)