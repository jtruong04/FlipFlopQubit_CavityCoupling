import numpy as np

gamma_e = 27.970                                     # GHz/T
gamma_n = 0.01723                                    # GHz/T
delta_gamma = -0.002                                 # Dimensionless
hyperfine = 0.1170                                   # GHz
dielectric = 11.7                                    # Dimensionless
Delta = gamma_e * delta_gamma / (gamma_e + gamma_n)  # Dimensionless
hbar = 1.0 / (2.0 *np.pi)                            # Dimensionless