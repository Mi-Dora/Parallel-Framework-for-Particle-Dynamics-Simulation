add_library(neo_particles
        dense_scheduler.cpp particle.cpp user_def/gravity_particle.cpp utils.cpp cutoff_scheduler.cpp centroid_scheduler.cpp
)
target_link_libraries(neo_particles PUBLIC MPI::MPI_CXX)
