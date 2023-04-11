#include "user_particle.h"
#include "utils.h"

using namespace std;

double user_particle::ke=8.9875517873681764e9;

vector<double> user_particle::compute_acceleration(user_particle& other) {
    double force;
    vector<double> orientation_vec = get_distance_vec(other.position, position);
    double len = get_length(orientation_vec);
    l2_norm(orientation_vec, len); // normalized the orientation_vec
    force = ke*charge*other.charge/len/len; // Coulomb's law
    acceleration[0] += force * orientation_vec[0] / mass;
    acceleration[1] += force * orientation_vec[1] / mass;
    acceleration[2] += force * orientation_vec[2] / mass;
}
