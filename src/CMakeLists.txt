add_library(pfpds
        particle.cpp utils.cpp world.cpp user_def/gravity_particle.cpp
        dense_scheduler.cpp
)
target_link_libraries(pfpds PUBLIC MPI::MPI_CXX)
