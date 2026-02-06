import numpy as np
from scipy.interpolate import interp1d
import warnings

# ========== Physical constants (using SI units) ==========
c = 2.99792458e8  # m/s
G = 6.67430e-11  # G in m^3/kg/s^2
Msun_turn = 1.98847e30  # Solar mass (kg)


# ========== TOV equations ==========
def TOV_eqns(r, M, p, pother, rho, G, c):
    """TOV equations, returns dp/dr and dM/dr"""
    if r <= 0:
        return 0, 0

    if M <= 0:
        M = 1e-20

    term1 = -(G * M * rho) / (r ** 2 * c ** 2)
    term2 = (1 + p / rho)
    term3 = (1 + 4 * np.pi * (r ** 3) * (p + pother) / (M * c ** 2))
    term4 = (1 - 2 * G * M / (r * c ** 2)) ** (-1)

    dpdr = term1 * term2 * term3 * term4
    dMdr = (4 * np.pi * r ** 2 * rho) / c ** 2

    return dpdr, dMdr


# ========== TOV equation solver function ==========
def solve_TOV(p1_center, p2_center, eos1, eos2, G, c):
    """Solve the coupled TOV equation system"""
    r_min = 1e-5
    M1_initial = 1
    M2_initial = 1

    dr_initial = 13
    max_steps = int(1e10)

    r = r_min
    M1 = M1_initial
    M2 = M2_initial
    p1 = p1_center
    p2 = p2_center

    p2_zero = (p2 == 0)
    p1_zero = (p1 == 0)

    R1 = 0
    R2 = 0

    # Create interpolation functions
    if len(eos1) > 1:
        sort_idx1 = np.argsort(eos1[:, 1])
        eos1_sorted = eos1[sort_idx1]
        rho_from_p1 = interp1d(eos1_sorted[:, 1], eos1_sorted[:, 0],
                               kind='linear', bounds_error=False, fill_value=0)
    else:
        rho_from_p1 = lambda p: eos1[0, 0] if len(eos1) > 0 else 0

    if len(eos2) > 1:
        sort_idx2 = np.argsort(eos2[:, 1])
        eos2_sorted = eos2[sort_idx2]
        rho_from_p2 = interp1d(eos2_sorted[:, 1], eos2_sorted[:, 0],
                               kind='linear', bounds_error=False, fill_value=0)
    else:
        rho_from_p2 = lambda p: eos2[0, 0] if len(eos2) > 0 else 0

    step = 0
    while (not p1_zero or not p2_zero):
        step += 1

        if p2 < 1e-100 and p1 < 1e-100:
            break

        if p2_zero and step == 1:
            M2_new = M2
        if p1_zero and step == 1:
            M1_new = M1

        if step == 1:
            dr = dr_initial

        if step > 1 and p1 < 1e-30:
            dr = r / 1000

        M = M1 + M2

        rho1 = float(rho_from_p1(p1))
        rho2 = float(rho_from_p2(p2))

        if not p1_zero and not p2_zero:
            k1_p1, k1_M1 = TOV_eqns(r, M, p1, p2, rho1, G, c)
            k1_p2, k1_M2 = TOV_eqns(r, M, p2, p1, rho2, G, c)

            k2_p1, k2_M1 = TOV_eqns(r + dr / 2, M, p1 + dr / 2 * k1_p1, p2 + dr / 2 * k1_p2, rho1, G, c)
            k2_p2, k2_M2 = TOV_eqns(r + dr / 2, M, p2 + dr / 2 * k1_p2, p1 + dr / 2 * k1_p1, rho2, G, c)

            k3_p1, k3_M1 = TOV_eqns(r + dr / 2, M, p1 + dr / 2 * k2_p1, p2 + dr / 2 * k2_p2, rho1, G, c)
            k3_p2, k3_M2 = TOV_eqns(r + dr / 2, M, p2 + dr / 2 * k2_p2, p1 + dr / 2 * k2_p1, rho2, G, c)

            k4_p1, k4_M1 = TOV_eqns(r + dr, M, p1 + dr * k3_p1, p2 + dr * k3_p2, rho1, G, c)
            k4_p2, k4_M2 = TOV_eqns(r + dr, M, p2 + dr * k3_p2, p1 + dr * k3_p1, rho2, G, c)

            p1_new = p1 + dr / 6 * (k1_p1 + 2 * k2_p1 + 2 * k3_p1 + k4_p1)
            M1_new = M1 + dr / 6 * (k1_M1 + 2 * k2_M1 + 2 * k3_M1 + k4_M1)
            p2_new = p2 + (dr / 6) * (k1_p2 + 2 * k2_p2 + 2 * k3_p2 + k4_p2)
            M2_new = M2 + (dr / 6) * (k1_M2 + 2 * k2_M2 + 2 * k3_M2 + k4_M2)

            if p1_new <= 0:
                p1_zero = True
                R1 = r
                p1_new = 0
                M1_new = M1

            if p2_new <= 0:
                p2_zero = True
                R2 = r
                p2_new = 0
                M2_new = M2

        elif not p1_zero and p2_zero:
            k1_p1, k1_M1 = TOV_eqns(r, M, p1, 0, rho1, G, c)
            k2_p1, k2_M1 = TOV_eqns(r + dr / 2, M, p1 + dr / 2 * k1_p1, 0, rho1, G, c)
            k3_p1, k3_M1 = TOV_eqns(r + dr / 2, M, p1 + dr / 2 * k2_p1, 0, rho1, G, c)
            k4_p1, k4_M1 = TOV_eqns(r + dr, M, p1 + dr * k3_p1, 0, rho1, G, c)

            p1_new = p1 + dr / 6 * (k1_p1 + 2 * k2_p1 + 2 * k3_p1 + k4_p1)
            M1_new = M1 + dr / 6 * (k1_M1 + 2 * k2_M1 + 2 * k3_M1 + k4_M1)

            if p1_new <= 0:
                p1_zero = True
                R1 = r
                p1_new = 0
                M1_new = M1

            p2_new = 0

        elif p1_zero and not p2_zero:
            k1_p2, k1_M2 = TOV_eqns(r, M, p2, 0, rho2, G, c)
            k2_p2, k2_M2 = TOV_eqns(r + dr / 2, M, p2 + dr / 2 * k1_p2, 0, rho2, G, c)
            k3_p2, k3_M2 = TOV_eqns(r + dr / 2, M, p2 + dr / 2 * k2_p2, 0, rho2, G, c)
            k4_p2, k4_M2 = TOV_eqns(r + dr, M, p2 + dr * k3_p2, 0, rho2, G, c)

            p2_new = p2 + (dr / 6) * (k1_p2 + 2 * k2_p2 + 2 * k3_p2 + k4_p2)
            M2_new = M2 + (dr / 6) * (k1_M2 + 2 * k2_M2 + 2 * k3_M2 + k4_M2)

            if p2_new <= 0:
                p2_zero = True
                R2 = r
                p2_new = 0
                M2_new = M2

            p1_new = 0

        r = r + dr
        p1 = p1_new
        p2 = p2_new
        M1 = M1_new
        M2 = M2_new

        if step > max_steps:
            warnings.warn("Maximum steps reached, exiting early")
            break

    M1_final = M1
    M2_final = M2

    if not p1_zero:
        R1 = r
    if not p2_zero:
        R2 = r

    return M1_final, M2_final, R1, R2


# ========== Main program ==========
def main():
    print("Solving coupled TOV equations system - Nuclear matter and dark matter (using MeV/fm^3 units)")

    try:
        eos1 = np.loadtxt('eos1.txt')
        eos2 = np.loadtxt('0.1 300.txt')
    except FileNotFoundError as e:
        print(f"Error: EOS data file not found: {e}")
        return

    if eos1.size == 0 or eos2.size == 0:
        raise ValueError('EOS data loading failed, please check file paths and contents')

    print("\nPlease input parameters:")
    rho1_min = float(input('Minimum central energy density of nuclear matter (MeV/fm^3): '))
    rho1_max = float(input('Maximum central energy density of nuclear matter (MeV/fm^3): '))
    rho1_step = float(input('Step size for nuclear matter central energy density (0 means keep as minimum): '))

    rho2_min = float(input('Minimum central energy density of dark matter (MeV/fm^3): '))
    rho2_max = float(input('Maximum central energy density of dark matter (MeV/fm^3): '))
    rho2_step = float(input('Step size for dark matter central energy density (0 means keep as minimum): '))

    conversion_factor = 1.602e32
    eos1 = eos1 * conversion_factor
    eos2 = eos2 * conversion_factor
    rho1_min = rho1_min * conversion_factor
    rho1_max = rho1_max * conversion_factor
    rho1_step = rho1_step * conversion_factor
    rho2_min = rho2_min * conversion_factor
    rho2_max = rho2_max * conversion_factor
    rho2_step = rho2_step * conversion_factor

    if rho1_step == 0:
        rho1_values = [rho1_min]
    else:
        rho1_values = np.arange(rho1_min, rho1_max + rho1_step, rho1_step)

    if rho2_step == 0:
        rho2_values = [rho2_min]
    else:
        rho2_values = np.arange(rho2_min, rho2_max + rho2_step, rho2_step)

    results = []

    for rho1_center in rho1_values:
        for rho2_center in rho2_values:
            print(
                f'Calculating rho1={rho1_center / conversion_factor:.2e} MeV/fm^3, rho2={rho2_center / conversion_factor:.2e} MeV/fm^3...')

            idx1 = np.argmin(np.abs(eos1[:, 0] - rho1_center))
            p1_center = eos1[idx1, 1]

            idx2 = np.argmin(np.abs(eos2[:, 0] - rho2_center))
            p2_center = eos2[idx2, 1]

            if p1_center is None or p2_center is None:
                warnings.warn(f'Cannot find initial pressure for rho1={rho1_center:.2e} or rho2={rho2_center:.2e}, skipping')
                continue

            try:
                M1_final, M2_final, R1, R2 = solve_TOV(p1_center, p2_center, eos1, eos2, G, c)

                M1_final_Msun = M1_final / Msun_turn
                M2_final_Msun = M2_final / Msun_turn
                R1_km = R1 * 1e-3
                R2_km = R2 * 1e-3

                results.append([
                    rho1_center / conversion_factor,
                    rho2_center / conversion_factor,
                    M1_final_Msun,
                    M2_final_Msun,
                    R1_km,
                    R2_km
                ])

                print(
                    f'Result: M1={M1_final_Msun:.3f} M☉, M2={M2_final_Msun:.3f} M☉, R1={R1_km:.2f} km, R2={R2_km:.2f} km')

            except Exception as e:
                warnings.warn(f'Error calculating rho1={rho1_center:.2e}, rho2={rho2_center:.2e}: {e}')

    if not results:
        raise ValueError('No results generated, please check input parameters and EOS data')

    results_array = np.array(results)

    print(
        '\nResults table: [rho1_center (MeV/fm^3), rho2_center (MeV/fm^3), M1_final (M☉), M2_final (M☉), R1 (km), R2 (km)]')
    print(results_array)

    try:
        header = 'rho1_center(MeV/fm^3)\trho2_center(MeV/fm^3)\tM1(Msun)\tM2(Msun)\tR1(km)\tR2(km)'
        np.savetxt('TOV_results.txt', results_array,
                   delimiter='\t', fmt='%.6e', header=header)
        print('Results successfully saved to TOV_results.txt')
    except Exception as e:
        print(f'Unable to write to file, please view results in command line window: {e}')
        print(results_array)

    return results_array


# ========== Run main program ==========
if __name__ == "__main__":
    results = main()
    print("\nProgram execution completed!")