# add_executable(TEST_load_store test_loads_store.cpp)
# target_link_libraries(TEST_load_store pfpds)

# add_executable(TEST_register_gravity test_register_gravity.cpp)
# target_link_libraries(TEST_register_gravity pfpds)

# add_executable(TEST_updates test_updates.cpp)
# target_link_libraries(TEST_updates pfpds)

add_executable(TEST_cuda test_cuda.cu)
target_link_libraries(TEST_cuda PUBLIC pfpds CUDA::cudart)