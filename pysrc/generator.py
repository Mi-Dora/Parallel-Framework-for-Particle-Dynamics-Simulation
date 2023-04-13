import numpy as np
import argparse


def gravity_generator(n_particle, n_dim, pos_range, mass_range, output_filename):
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
    np.savetxt(output_filename, mat, delimiter=' ', newline='\n', header=header, comments='')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_particle", type=int, default=1000)
    parser.add_argument("--n_dim", type=int, default=3)
    parser.add_argument("--pos_range", type=list, default=[-100, 100, -100, 100])
    parser.add_argument("--mass_range", type=list, default=[0, 0.01])
    parser.add_argument("--output_filename", type=str, required=True)
    parser.parse_args()
    
    n_particle = 1000
    n_dim = 3
    if len(parser.pos_range) == 2:
        x0, x1 = parser.pos_range[0], parser.pos_range[1]
        parser.pos_range = [x0, x1, x0, x1]
    assert len(parser.pos_range) == 4
    gravity_generator(parser.n_particle, parser.n_dim, parser.pos_range, parser.mass_range, parser.output_filename)
