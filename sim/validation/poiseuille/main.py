# ---------------------------------------
# Analytical Solution
# ---------------------------------------


def get_analytic(_mesh):
    Wz_a = np.zeros(NN, dtype="float64")
    Psi_a = np.zeros(NN, dtype="float64")
    vx_a = np.zeros(NN, dtype="float64")

    for i in range(_mesh.num_nodes):
        vx_a[i] = 0
        psi_a[i] = 0
        Wz_a[i] = 0

    return Wz_a, Psi_a, vx_a
