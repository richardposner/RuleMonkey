# examples/CMakeLists.txt — wired up only when RULEMONKEY_BUILD_EXAMPLES=ON.
#
# This directory is not part of the engine library; it exists to show
# how a host project can link librulemonkey via add_subdirectory() and
# how the public-API headers are consumed.

add_executable(rm_embed_example embed.cpp)
target_link_libraries(rm_embed_example PRIVATE rulemonkey)
