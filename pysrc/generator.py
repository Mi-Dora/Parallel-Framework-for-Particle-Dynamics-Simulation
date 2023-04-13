import numpy as np


def gravity_generator(n_particle, n_dim, pos_range, mass_range):
    """
    :param n_particle: int
    :param n_dim: int
    :param pos_range: [float, float]
    :param mass_range: [float, float]
    :return: void
    """
    header = "{} {}".format(n_particle, n_dim)
    mu = 0
    sigma = 0.1
    id = np.linspace(1, n_particle, n_particle).reshape((n_particle, 1))
    pos = np.random.uniform(pos_range[0], pos_range[1], (n_particle, n_dim))
    vel = np.random.normal(mu, sigma, (n_particle, n_dim))
    acc = np.random.normal(mu, sigma, (n_particle, n_dim))
    n_feat = np.ones((n_particle, 1))
    mass = np.random.uniform(mass_range[0], mass_range[1], (n_particle, 1))
    mat = np.concatenate((id, pos, vel, acc, n_feat, mass), axis=1)
    np.savetxt("particle_init.txt", mat, delimiter=' ', newline='\n', header=header, comments='')


if __name__ == "__main__":
    n_particle = 1000
    n_dim = 3
    pos_range = [0, 100]
    mass_range = [0, 0.01]
    gravity_generator(n_particle, n_dim, pos_range, mass_range)
