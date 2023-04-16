add_library(pfpds
        particle.cpp utils.cpp world.cpp user_def/gravity_particle.cpp kernels.cu device.cpp
)
target_link_libraries(pfpds PUBLIC CUDA::cudart)