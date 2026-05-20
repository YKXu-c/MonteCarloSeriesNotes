import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import subprocess
import json
import sys
from pathlib import Path


def kagome_kondo_ising():
    """Draw kagome lattice (Co sublattice) with Tb localized spins and c-f hybridization."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))

    # Kagome lattice vectors
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3) / 2])

    # Kagome basis: 3 sites per unit cell at midpoints of the triangular lattice bonds
    # Site A: (a1)/2, Site B: (a2)/2, Site C: (a1+a2)/2
    basis = [a1 / 2, a2 / 2, (a1 + a2) / 2]

    # Generate kagome sites for a 4x4 patch
    co_sites = []
    n_cells = 4
    for i in range(n_cells):
        for j in range(n_cells):
            origin = i * a1 + j * a2
            for b in basis:
                co_sites.append(origin + b)
    co_sites = np.array(co_sites)

    # Triangular lattice sites (Tb positions) — at the hexagonal centers
    # These are at the original triangular lattice points
    tb_sites = []
    for i in range(n_cells + 1):
        for j in range(n_cells + 1):
            tb_sites.append(i * a1 + j * a2)
    # Also add sites shifted by (2a1+a2)/3 and (a1+2a2)/3 for the hexagonal centers
    hex_centers = []
    for i in range(n_cells):
        for j in range(n_cells):
            # Up triangles center
            c1 = (i * a1 + j * a2 + (i + 1) * a1 + j * a2 + i * a1 + (j + 1) * a2) / 3
            c2 = ((i + 1) * a1 + j * a2 + (i + 1) * a1 + (j + 1) * a2 + i * a1 + (j + 1) * a2) / 3
            hex_centers.extend([c1, c2])
    hex_centers = np.array(hex_centers)

    # Plot Co kagome bonds (nearest-neighbor connections on kagome)
    bond_len = 0.5 + 0.01  # slightly generous threshold
    for i, si in enumerate(co_sites):
        for j, sj in enumerate(co_sites):
            if j > i:
                d = np.linalg.norm(si - sj)
                if d < bond_len:
                    ax.plot([si[0], sj[0]], [si[1], sj[1]], 'b-', lw=1.2, alpha=0.5)

    # Plot Co sites (kagome lattice)
    ax.scatter(co_sites[:, 0], co_sites[:, 1], s=80, c='steelblue', zorder=5,
               edgecolors='navy', linewidths=0.8, label=r'Co (kagome)')

    # Plot Tb localized spin sites at hexagonal centers
    ax.scatter(hex_centers[:, 0], hex_centers[:, 1], s=180, c='orangered', zorder=6,
               edgecolors='darkred', linewidths=0.8, label=r'Tb (localized)', marker='o')

    # Draw Ising spin arrows on Tb sites
    arrow_scale = 0.15
    for site in hex_centers:
        # Alternate up/down to show ferrimagnetic/antiferromagnetic ordering
        direction = 1 if (hash(tuple(site.round(3))) % 2 == 0) else -1
        ax.annotate('', xy=(site[0], site[1] + direction * arrow_scale),
                     xytext=(site[0], site[1] - direction * arrow_scale),
                     arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5))

    # Draw c-f hybridization (dashed lines between nearest Co-Tb pairs)
    for tb in hex_centers:
        for co in co_sites:
            d = np.linalg.norm(tb - co)
            if d < 0.35:
                ax.plot([tb[0], co[0]], [tb[1], co[1]], '--', color='gray',
                        lw=0.7, alpha=0.4, zorder=3)

    ax.set_xlim(-0.3, 4.2)
    ax.set_ylim(-0.3, 3.8)
    ax.set_aspect('equal')
    ax.axis('off')

    # Legend
    co_patch = mpatches.Patch(color='steelblue', label='Co (kagome sublattice)')
    tb_patch = mpatches.Patch(color='orangered', label='Tb (localized spins)')
    hybrid_line = plt.Line2D([0], [0], linestyle='--', color='gray', label='c-f hybridization')
    ax.legend(handles=[co_patch, tb_patch, hybrid_line], loc='lower right', fontsize=10,
              framealpha=0.9, edgecolor='gray')

    # Title
    ax.set_title('Kondo-Ising Model on Kagome Lattice\n'
                 r'$H = -t\sum_{\langle ij\rangle\sigma} c_{i\sigma}^\dagger c_{j\sigma}'
                 r' + J_K\sum_i S_i^z s_i^z - J\sum_{\langle ij\rangle} S_i^z S_j^z$',
                 fontsize=13, pad=15)

    plt.tight_layout()
    outpath = Path('pics/kagome_kondo_ising.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


def kagome_kondo_ising_v2():
    """Isometric bilayer Kondo-Ising schematic on kagome lattice, side-view perspective.

    Like FIG.1 of the AM paper: isometric 3D view with Co kagome layer on top
    and Tb localized spin layer below, connected by vertical Kondo coupling J_K.
    The kagome lattice is drawn in pseudo-3D by applying an oblique projection.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Oblique projection: map (x, y, z) -> (x + 0.35*z, y + 0.35*z)
    # Layer separation in the z-direction
    layer_sep = 2.8  # larger gap: Kondo-Ising model, not real crystal

    def proj(x, y, z):
        return x + 0.3 * z, y + 0.35 * z

    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3) / 2])
    basis = [a1 / 2, a2 / 2, (a1 + a2) / 2]
    n_cells = 3

    # Generate Co kagome sites (top layer, z=1)
    co_sites_2d = []
    for i in range(n_cells):
        for j in range(n_cells):
            origin = i * a1 + j * a2
            for b in basis:
                co_sites_2d.append(origin + b)
    co_sites_2d = np.array(co_sites_2d)
    co_sites = np.array([proj(s[0], s[1], layer_sep) for s in co_sites_2d])

    # Generate Tb triangular lattice sites (bottom layer, z=0)
    tb_sites_2d = []
    for i in range(n_cells + 1):
        for j in range(n_cells + 1):
            tb_sites_2d.append(i * a1 + j * a2)
    tb_sites_2d = np.array(tb_sites_2d)
    tb_sites = np.array([proj(s[0], s[1], 0) for s in tb_sites_2d])

    # --- Draw Tb layer first (bottom, behind) ---
    # Ising bonds between nearest Tb sites
    tb_bond_len = 1.0 + 0.05
    for i_s, si in enumerate(tb_sites_2d):
        for j_s, sj in enumerate(tb_sites_2d):
            if j_s > i_s:
                d = np.linalg.norm(si - sj)
                if d < tb_bond_len:
                    p1 = proj(si[0], si[1], 0)
                    p2 = proj(sj[0], sj[1], 0)
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], '-', color='#b2182b',
                            lw=1.0, alpha=0.3, zorder=1)

    # Tb sites with AFM arrows (blue up, gold down)
    arrow_len = 0.22
    for idx, site in enumerate(tb_sites_2d):
        px, py = proj(site[0], site[1], 0)
        i_idx = int(round(site[0]))
        j_idx = int(round((2 * site[1] - site[0]) / np.sqrt(3)))
        direction = 1 if (i_idx + j_idx) % 2 == 0 else -1
        color = '#2166ac' if direction == 1 else '#d4a017'
        ax.annotate('', xy=(px, py + direction * arrow_len),
                     xytext=(px, py - direction * arrow_len * 0.3),
                     arrowprops=dict(arrowstyle='->', color=color, lw=2.0),
                     zorder=6)

    ax.scatter(tb_sites[:, 0], tb_sites[:, 1], s=140, c='#fddbc7', zorder=5,
               edgecolors='#b2182b', linewidths=1.0, marker='o')

    # --- Kondo coupling: vertical dashed lines from Tb to nearest Co neighbors ---
    # Like v1: connect each Tb to all Co sites within a short cutoff (surrounding kagome sites)
    jk_cutoff = 0.6
    for idx, ts in enumerate(tb_sites_2d):
        tp = proj(ts[0], ts[1], 0)
        for cs in co_sites_2d:
            d = np.linalg.norm(cs - ts)
            if d < jk_cutoff:
                cp = proj(cs[0], cs[1], layer_sep)
                ax.plot([tp[0], cp[0]], [tp[1], cp[1]], '--', color='#666666',
                        lw=0.8, alpha=0.5, zorder=3)

    # Label one Kondo coupling line
    label_idx = len(tb_sites_2d) // 2 + n_cells // 2 + 1
    if label_idx < len(tb_sites_2d):
        ts = tb_sites_2d[label_idx]
        tp = proj(ts[0], ts[1], 0)
        # Find one nearest Co for label placement
        min_d = float('inf')
        closest_co = None
        for cs in co_sites_2d:
            d = np.linalg.norm(cs - ts)
            if d < min_d:
                min_d = d
                closest_co = cs
        if min_d < jk_cutoff:
            cp = proj(closest_co[0], closest_co[1], layer_sep)
            mid = ((tp[0] + cp[0]) / 2 + 0.15, (tp[1] + cp[1]) / 2)
            ax.annotate(r'$J_K$', xy=mid, fontsize=12, color='#444444', fontweight='bold')

    # --- Draw Co layer (top, in front) ---
    bond_len = 0.5 + 0.02
    for i_s, si in enumerate(co_sites_2d):
        for j_s, sj in enumerate(co_sites_2d):
            if j_s > i_s:
                d = np.linalg.norm(si - sj)
                if d < bond_len:
                    p1 = proj(si[0], si[1], layer_sep)
                    p2 = proj(sj[0], sj[1], layer_sep)
                    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], '-', color='#2166ac',
                            lw=1.8, alpha=0.7, zorder=8)

    # Co sublattice coloring: 3 colors for 3 kagome basis sites
    basis_colors = ['#2166ac', '#4393c3', '#92c5de']  # dark, medium, light blue
    colors_co = [basis_colors[b] for i in range(n_cells) for j in range(n_cells) for b in range(3)]
    ax.scatter(co_sites[:, 0], co_sites[:, 1], s=100, c=colors_co, zorder=10,
               edgecolors='#053061', linewidths=0.8)

    # Hopping label on one bond
    ax.annotate(r'$t$', xy=(proj(0.25, 0.0, layer_sep)[0],
                             proj(0.25, 0.0, layer_sep)[1] + 0.12),
                fontsize=13, color='#2166ac', fontweight='bold', zorder=11)

    # Ising label
    ax.annotate(r'$J$', xy=(proj(0.5, 0.0, 0)[0],
                             proj(0.5, 0.0, 0)[1] - 0.25),
                fontsize=13, color='#b2182b', fontweight='bold', zorder=11)

    # --- Layer labels ---
    co_center = np.mean(co_sites, axis=0)
    tb_center = np.mean(tb_sites, axis=0)
    ax.annotate('Co kagome (itinerant)', xy=(co_center[0] + 1.5, co_center[1] + 0.3),
                fontsize=11, color='#2166ac', fontweight='bold',
                ha='center', zorder=12)
    ax.annotate('Tb (localized spins)', xy=(tb_center[0] + 1.5, tb_center[1] - 0.3),
                fontsize=11, color='#b2182b', fontweight='bold',
                ha='center', zorder=12)

    ax.set_xlim(-0.5, 5.0)
    ax.set_ylim(-0.8, 4.3)  # reduced from 5.0 to trim top whitespace
    ax.set_aspect('equal')
    ax.axis('off')

    # --- Hamiltonian at bottom ---
    ax.text(0.5, -0.06,
            r'$H = -t\!\sum_{\langle ij\rangle\sigma} c_{i\sigma}^\dagger c_{j\sigma}'
            r'\;+\; J_K\!\sum_i S_i^z s_i^z'
            r'\;-\; J\!\sum_{\langle ij\rangle} S_i^z S_j^z$',
            transform=ax.transAxes, ha='center', fontsize=14,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#f7f7f7', edgecolor='#999999'))

    # --- Legend ---
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#2166ac',
                    markeredgecolor='#053061', markersize=9, label='Co (sublattice A)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#4393c3',
                    markeredgecolor='#053061', markersize=9, label='Co (sublattice B)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#92c5de',
                    markeredgecolor='#053061', markersize=9, label='Co (sublattice C)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#fddbc7',
                    markeredgecolor='#b2182b', markersize=11, label='Tb (localized)'),
        plt.Line2D([0], [0], linestyle='-', color='#2166ac', lw=2, label=r'Hopping $t$'),
        plt.Line2D([0], [0], linestyle='-', color='#b2182b', lw=1.5, alpha=0.5,
                    label=r'Ising $J$'),
        plt.Line2D([0], [0], linestyle='--', color='#666666', lw=1.2,
                    label=r'Kondo $J_K$'),
        plt.Line2D([0], [0], marker='$↑$', color='#2166ac', markersize=12,
                    linestyle='none', label='Spin up'),
        plt.Line2D([0], [0], marker='$↓$', color='#d4a017', markersize=12,
                    linestyle='none', label='Spin down'),
    ]
    ax.legend(handles=legend_elements, loc='center left', fontsize=9,
              framealpha=0.95, edgecolor='#999999', bbox_to_anchor=(-0.05, 0.5))
    # bbox_to_anchor=(-0.05, 0.5)  # for fine-tuning legend position

    plt.tight_layout()
    outpath = Path('pics/kagome_kondo_ising_v2.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


# ---------------------------------------------------------------------------
# MC data plotting infrastructure
# ---------------------------------------------------------------------------

BINARY_DIR = Path(__file__).parent / 'examples' / 'build'


def parse_mc_output(stdout: str) -> dict:
    """Parse C++ MC output: JSON header + tab-separated data + optional time series."""
    result = {'params': {}, 'observables': {}, 'time_series': {}}
    lines = stdout.strip().split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('# {'):
            result['params'] = json.loads(line[2:])
        elif line.startswith('# algorithm:'):
            result['algorithm'] = line.split(':', 1)[1].strip()
        elif line.startswith('# seed:'):
            result['seed'] = int(line.split(':', 1)[1].strip())
        elif line.startswith('# time_series_begin:'):
            ts_name = line.split(':', 1)[1].strip()
            ts_values = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('# time_series_end'):
                try:
                    ts_values.append(float(lines[i].strip()))
                except ValueError:
                    pass
                i += 1
            result['time_series'][ts_name] = np.array(ts_values)
        elif line.startswith('#') or not line or line.startswith('observable'):
            pass
        else:
            parts = line.split('\t')
            if len(parts) >= 2:
                result['observables'][parts[0]] = {
                    'mean': float(parts[1]),
                    'variance': float(parts[2]) if len(parts) > 2 else 0.0,
                    'mean2': float(parts[3]) if len(parts) > 3 else 0.0,
                }
        i += 1
    return result


def run_mc_binary(algorithm: str, L: int, J: float, T: float,
                  sweeps: int = 5000, therm: int = 2000,
                  Jp: float = 0.0, seed: int = 0,
                  all_up: bool = False, ts: bool = False,
                  auto_therm: bool = False) -> dict:
    """Run a C++ MC binary and return parsed results."""
    binary = BINARY_DIR / algorithm
    if not binary.exists():
        raise FileNotFoundError(f'Binary not found: {binary}')
    cmd = [str(binary), '--L', str(L), '--J', str(J), '--T', str(T),
           '--sweeps', str(sweeps), '--therm', str(therm)]
    if Jp != 0.0:
        cmd += ['--Jp', str(Jp)]
    if seed != 0:
        cmd += ['--seed', str(seed)]
    if all_up:
        cmd += ['--all-up']
    if ts:
        cmd += ['--ts']
    if auto_therm:
        cmd += ['--auto-therm']
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise RuntimeError(f'{algorithm} failed: {proc.stderr}')
    result = parse_mc_output(proc.stdout)
    if proc.stderr:
        result['stderr'] = proc.stderr
    return result


def temperature_sweep(algorithm: str, T_list, L: int, J: float = 1.0,
                      sweeps: int = 5000, therm: int = 2000, **kwargs) -> list:
    """Run algorithm at each temperature, return list of result dicts."""
    results = []
    for T in T_list:
        r = run_mc_binary(algorithm, L, J, T, sweeps, therm, **kwargs)
        r['T'] = T
        results.append(r)
        print(f'  {algorithm} T={T:.2f} |m|={r["observables"].get("abs_magnetization", {}).get("mean", "N/A")}')
    return results


def onsager_exact_M(T_arr, J=1.0):
    """Onsager exact spontaneous magnetization for 2D Ising (NN only)."""
    T_c = 2.0 / np.log(1.0 + np.sqrt(2.0))
    M = np.zeros_like(T_arr)
    mask = T_arr < T_c
    with np.errstate(divide='ignore', invalid='ignore'):
        sinh_val = np.sinh(2.0 * J / T_arr[mask])
        ratio = 1.0 / sinh_val**4
    valid = ratio < 1.0
    M[mask] = np.where(valid, (1.0 - ratio)**(1.0/8.0), 0.0)
    return M


def find_tc_numerical(T, M):
    """Estimate Tc from MC data by finding T where -dM/dT is maximal."""
    idx = np.argsort(T)
    T_sorted, M_sorted = T[idx], M[idx]
    T_fine = np.linspace(T_sorted[0], T_sorted[-1], 1000)
    M_fine = np.interp(T_fine, T_sorted, M_sorted)
    dM_dT = np.gradient(M_fine, T_fine)
    tc_idx = np.argmin(dM_dT)
    return T_fine[tc_idx]


T_SWEEP_LIST = np.array([0.05, 0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8,
                         2.0, 2.1, 2.2, 2.25, 2.27, 2.29, 2.3, 2.35,
                         2.4, 2.5, 2.7, 3.0, 3.5, 4.0])


def plot_metropolis_mt():
    """Metropolis |m| vs T with Onsager exact overlay — verification plot."""
    print('Running Metropolis temperature sweep...')
    results = temperature_sweep('metropolis', T_SWEEP_LIST, L=16, J=1.0,
                                sweeps=10000, therm=5000)

    mc_T = np.array([r['T'] for r in results])
    mc_m = np.array([r['observables']['abs_magnetization']['mean'] for r in results])

    T_fine = np.linspace(0.0, 4.5, 300)
    M_exact = onsager_exact_M(T_fine, J=1.0)
    T_c_onsager = 2.0 / np.log(1.0 + np.sqrt(2.0))
    T_c_num = find_tc_numerical(mc_T, mc_m)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T_fine, M_exact, 'k-', lw=2, label='Onsager exact')
    ax.plot(mc_T, mc_m, 'ro', ms=6, label='Metropolis (L=16)', zorder=5)
    ax.axvline(T_c_onsager, color='gray', ls='--', lw=1, alpha=0.5)
    ax.axvline(T_c_num, color='red', ls=':', lw=1.5, alpha=0.8,
               label=rf'$T_c^{{\mathrm{{num}}}} = {T_c_num:.3f}$')
    ax.set_xlabel(r'$T$', fontsize=14)
    ax.set_ylabel(r'$\langle |m| \rangle$', fontsize=14)
    ax.set_title('Metropolis MC vs Onsager Exact Solution', fontsize=14)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_xlim(0.0, 4.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = Path('pics/metropolis_mt.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


def plot_thermalization(T=2.5, sweeps=100000, therm=0):
    """Dual-start convergence: plot ⟨m⟩ from all-up and random initial conditions."""
    print('Running thermalization convergence check...')
    #T = 2.5
    #sweeps = 10000
    #therm = 0  # no fixed thermalization — watch the full trajectory

    r_up = run_mc_binary('metropolis', L=16, J=1.0, T=T, sweeps=sweeps, therm=therm,
                         all_up=True, ts=True, seed=42)
    r_rand = run_mc_binary('metropolis', L=16, J=1.0, T=T, sweeps=sweeps, therm=therm,
                           all_up=False, ts=True, seed=123)

    ts_up = r_up['time_series'].get('magnetization', np.array([]))
    ts_rand = r_rand['time_series'].get('magnetization', np.array([]))

    if len(ts_up) == 0 or len(ts_rand) == 0:
        print('Warning: no time series data. Check --ts flag.')
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ts_up, 'b-', lw=0.8, alpha=0.7, label='All-up init')
    ax.plot(ts_rand, 'r-', lw=0.8, alpha=0.7, label='Random init')
    ax.axhline(y=0, color='gray', ls=':', lw=0.5)
    ax.set_xlabel('Sweep', fontsize=14)
    ax.set_ylabel(r'$\langle m \rangle$', fontsize=14)
    ax.set_title(f'Thermalization Convergence (T={T}, L=16)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = Path('pics/thermalization_convergence_2.5K.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


def plot_sw_mt():
    """Swendsen-Wang |m| vs T with Onsager exact overlay — verification plot."""
    print('Running Swendsen-Wang temperature sweep...')
    results = temperature_sweep('swendsen_wang', T_SWEEP_LIST, L=16, J=1.0,
                                sweeps=5000, therm=2000)

    mc_T = np.array([r['T'] for r in results])
    mc_m = np.array([r['observables']['abs_magnetization']['mean'] for r in results])

    T_fine = np.linspace(0.0, 4.5, 300)
    M_exact = onsager_exact_M(T_fine, J=1.0)
    T_c_onsager = 2.0 / np.log(1.0 + np.sqrt(2.0))
    T_c_num = find_tc_numerical(mc_T, mc_m)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T_fine, M_exact, 'k-', lw=2, label='Onsager exact')
    ax.plot(mc_T, mc_m, 's', color='#2196F3', ms=6, label='Swendsen-Wang (L=16)', zorder=5)
    ax.axvline(T_c_onsager, color='gray', ls='--', lw=1, alpha=0.5)
    ax.axvline(T_c_num, color='#2196F3', ls=':', lw=1.5, alpha=0.8,
               label=rf'$T_c^{{\mathrm{{num}}}} = {T_c_num:.3f}$')
    ax.set_xlabel(r'$T$', fontsize=14)
    ax.set_ylabel(r'$\langle |m| \rangle$', fontsize=14)
    ax.set_title('Swendsen-Wang MC vs Onsager Exact Solution', fontsize=14)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_xlim(0.0, 4.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = Path('pics/sw_mt.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


def plot_wolff_mt():
    """Wolff |m| vs T with Onsager exact overlay — verification plot."""
    print('Running Wolff temperature sweep...')
    results = temperature_sweep('wolff', T_SWEEP_LIST, L=16, J=1.0,
                                sweeps=5000, therm=2000)

    mc_T = np.array([r['T'] for r in results])
    mc_m = np.array([r['observables']['abs_magnetization']['mean'] for r in results])

    T_fine = np.linspace(0.0, 4.5, 300)
    M_exact = onsager_exact_M(T_fine, J=1.0)
    T_c_onsager = 2.0 / np.log(1.0 + np.sqrt(2.0))
    T_c_num = find_tc_numerical(mc_T, mc_m)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T_fine, M_exact, 'k-', lw=2, label='Onsager exact')
    ax.plot(mc_T, mc_m, '^', color='#4CAF50', ms=6, label='Wolff (L=16)', zorder=5)
    ax.axvline(T_c_onsager, color='gray', ls='--', lw=1, alpha=0.5)
    ax.axvline(T_c_num, color='#4CAF50', ls=':', lw=1.5, alpha=0.8,
               label=rf'$T_c^{{\mathrm{{num}}}} = {T_c_num:.3f}$')
    ax.set_xlabel(r'$T$', fontsize=14)
    ax.set_ylabel(r'$\langle |m| \rangle$', fontsize=14)
    ax.set_title('Wolff MC vs Onsager Exact Solution', fontsize=14)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_xlim(0.0, 4.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = Path('pics/wolff_mt.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'Saved: {outpath}')


def autocorrelation_time(ts):
    """Integrated autocorrelation time from magnetization time series.

    τ_int = 1/2 + Σ_{k=1}^M ρ(k), truncated at first negative ρ(k).
    """
    n = len(ts)
    mean = np.mean(ts)
    var = np.var(ts)
    if var == 0:
        return 1.0
    ts_ctr = ts - mean
    max_lag = min(n // 2, 500)
    tau = 0.5
    for lag in range(1, max_lag):
        rho = np.mean(ts_ctr[:n - lag] * ts_ctr[lag:]) / var
        if rho < 0:
            break
        tau += rho
    return tau


def compute_tau_int(algorithm, L, T, sweeps, therm, seed=42):
    """Run binary with --ts and return integrated autocorrelation time of |m|."""
    r = run_mc_binary(algorithm, L=L, J=1.0, T=T, sweeps=sweeps, therm=therm,
                      all_up=False, ts=True, seed=seed)
    ts_m = r['time_series'].get('magnetization', np.array([]))
    ts_abs_m = np.abs(ts_m)
    if len(ts_abs_m) == 0:
        print(f'  Warning: no time series for {algorithm} L={L}')
        return float('nan')
    tau = autocorrelation_time(ts_abs_m)
    return tau


def plot_dynamic_exponent():
    """Compute dynamic critical exponent z for all three algorithms.

    Runs at Tc for multiple L, computes τ_int(L), fits τ ∝ L^z on log-log scale.
    """
    T_c = 2.0 / np.log(1.0 + np.sqrt(2.0))
    print(f'\n=== Dynamic Critical Exponent z (T_c = {T_c:.4f}) ===')

    # L values and sweep counts per algorithm
    L_metro = np.array([8, 16, 32])
    L_cluster = np.array([8, 16, 32, 64])
    sweeps_metro, therm_metro = 50000, 10000
    sweeps_cluster, therm_cluster = 20000, 5000

    algorithms = {
        'metropolis':  ('Metropolis', L_metro, sweeps_metro, therm_metro, 'o', 'red'),
        'swendsen_wang': ('Swendsen-Wang', L_cluster, sweeps_cluster, therm_cluster, 's', '#2196F3'),
        'wolff':         ('Wolff', L_cluster, sweeps_cluster, therm_cluster, '^', '#4CAF50'),
    }

    results = {}
    for alg, (label, L_vals, sweeps, therm, marker, color) in algorithms.items():
        print(f'\n{alg}:')
        tau_vals = []
        for L in L_vals:
            tau = compute_tau_int(alg, L, T_c, sweeps, therm)
            tau_vals.append(tau)
            print(f'  L={L:3d}  τ_int={tau:.2f}')
        tau_vals = np.array(tau_vals)
        log_L = np.log(L_vals)
        log_tau = np.log(tau_vals)
        # Linear fit: log(τ) = z · log(L) + c
        z, c = np.polyfit(log_L, log_tau, 1)
        print(f'  → z = {z:.3f}')
        results[alg] = (label, L_vals, tau_vals, z, c, marker, color)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    for alg, (label, L_vals, tau_vals, z, c, marker, color) in results.items():
        L_fine = np.logspace(np.log10(L_vals[0] * 0.8), np.log10(L_vals[-1] * 1.2), 50)
        tau_fit = np.exp(c) * L_fine ** z
        ax.loglog(L_vals, tau_vals, marker, color=color, ms=8, zorder=5,
                  label=f'{label} ($z = {z:.2f}$)')
        ax.loglog(L_fine, tau_fit, '--', color=color, lw=1.2, alpha=0.6)

    ax.set_xlabel(r'$L$', fontsize=14)
    ax.set_ylabel(r'$\tau_{\mathrm{int}}$', fontsize=14)
    ax.set_title(r'Dynamic Critical Exponent: $\tau_{\mathrm{int}} \propto L^z$', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    # Reference lines for literature z values
    ax.text(0.05, 0.95, r'Literature: $z_{\mathrm{Metro}}\approx 2.17,\; z_{\mathrm{SW}}\approx 0.35,\; z_{\mathrm{Wolff}}\approx 0.25$',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    outpath = Path('pics/dynamic_exponent_z.png')
    outpath.parent.mkdir(exist_ok=True)
    fig.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'\nSaved: {outpath}')


if __name__ == '__main__':
    kagome_kondo_ising()
    kagome_kondo_ising_v2()
    plot_metropolis_mt()
    plot_thermalization()
    plot_sw_mt()
    plot_wolff_mt()
    plot_dynamic_exponent()
