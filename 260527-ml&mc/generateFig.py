"""Generate figures for FFO-CO + Frustrated Spin Systems slides."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

PIC_DIR = Path(__file__).parent / "pics"
PIC_DIR.mkdir(exist_ok=True)
DPI = 200


def kagome_afm():
    """2D kagome lattice with AFM coupling and frustration annotation."""
    fig, ax = plt.subplots(1, 1, figsize=(5, 4.5))

    # Kagome lattice: triangles pointing up and down
    a1 = np.array([1.0, 0.0])
    a2 = np.array([0.5, np.sqrt(3) / 2])

    # Basis vectors for kagome (3 sites per unit cell)
    b1 = a1 / 2
    b2 = a2 / 2
    b3 = (a1 + a2) / 2

    nx, ny = 4, 3
    sites = []
    for i in range(-1, nx + 1):
        for j in range(-1, ny + 1):
            origin = i * a1 + j * a2
            sites.append(origin + b1)
            sites.append(origin + b2)
            sites.append(origin + b3)

    sites = np.array(sites)

    # Clip to visible area
    mask = (sites[:, 0] > -0.5) & (sites[:, 0] < nx + 0.5) & \
           (sites[:, 1] > -0.5) & (sites[:, 1] < ny * np.sqrt(3) / 2 + 0.5)
    sites = sites[mask]

    # Draw bonds (nearest neighbors)
    drawn = set()
    for i, si in enumerate(sites):
        for j, sj in enumerate(sites):
            if j <= i:
                continue
            d = np.linalg.norm(si - sj)
            if abs(d - 0.5) < 0.05:  # NN distance = 0.5
                pair = (min(i, j), max(i, j))
                if pair not in drawn:
                    drawn.add(pair)
                    ax.plot([si[0], sj[0]], [si[1], sj[1]],
                            'k-', linewidth=0.8, alpha=0.6)

    # Draw sites as up arrows (spins) — all pointing up, AFM frustrated
    ax.scatter(sites[:, 0], sites[:, 1], s=80, c='#004983',
               zorder=5, edgecolors='k', linewidths=0.5)

    # Highlight one frustrated triangle
    tri_center = np.array([1.0, np.sqrt(3) / 6])
    tri_verts = [b1, b2, b3]  # vertices of one up triangle
    triangle = plt.Polygon(tri_verts, fill=False, edgecolor='red',
                           linewidth=2.5, linestyle='--', zorder=6)
    ax.add_patch(triangle)

    # Annotate frustration
    ax.annotate('Frustrated!\nCannot satisfy\nall AFM bonds',
                xy=(0.5, 0.15), fontsize=9, color='red',
                ha='center', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFE0E0',
                          edgecolor='red', alpha=0.9))

    ax.set_xlim(-0.3, 3.8)
    ax.set_ylim(-0.2, 2.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('2D Kagome AFM', fontsize=12, fontweight='bold', color='#004983')

    fig.tight_layout()
    fig.savefig(PIC_DIR / "kagome_afm.png", dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {PIC_DIR / 'kagome_afm.png'}")


def pyrochlore_frust():
    """3D pyrochlore: single tetrahedron with AFM frustration, semi-transparent faces."""
    fig = plt.figure(figsize=(5, 4.5))
    ax = fig.add_subplot(111, projection='3d')

    # Regular tetrahedron vertices (centered at origin)
    v = np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]
    ], dtype=float) * 0.5

    # Semi-transparent faces
    faces = [[v[0], v[1], v[2]],
             [v[0], v[1], v[3]],
             [v[0], v[2], v[3]],
             [v[1], v[2], v[3]]]
    face_col = Poly3DCollection(faces, alpha=0.10, facecolor='#4682B4',
                                edgecolor='#004983', linewidth=2.0)
    ax.add_collection3d(face_col)

    # AFM spin arrows: best arrangement is "2-in, 2-out" (spin ice rule)
    # Still leaves 2 out of 6 bonds frustrated
    spin_dirs = np.array([
        [0, 0, 0.55],       # v0: up (in)
        [0.35, 0.1, -0.35],  # v1: down-ish (out)
        [-0.25, 0.35, -0.35],  # v2: down-ish (out)
        [-0.1, -0.35, 0.45]   # v3: up-ish (in)
    ])

    for i in range(4):
        ax.quiver(v[i, 0], v[i, 1], v[i, 2],
                  spin_dirs[i, 0], spin_dirs[i, 1], spin_dirs[i, 2],
                  color='red', arrow_length_ratio=0.25, linewidth=2.2)
        ax.scatter3D([v[i, 0]], [v[i, 1]], [v[i, 2]], s=100,
                     c='#004983', edgecolors='k', linewidths=0.5, zorder=5)

    # Mark frustrated bonds (dashed red)
    # Bonds between same-direction spins are frustrated
    frustrated_pairs = [(0, 3), (1, 2)]  # both "in" or both "out"
    for i, j in frustrated_pairs:
        ax.plot3D([v[i, 0], v[j, 0]], [v[i, 1], v[j, 1]], [v[i, 2], v[j, 2]],
                  'r--', linewidth=1.2, alpha=0.6)

    ax.text(0, 0, -0.85,
            'AFM Frustrated\nBest: "2-in, 2-out"\n2 of 6 bonds unsatisfied',
            fontsize=8, color='red', ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFE0E0',
                      edgecolor='red', alpha=0.85))

    ax.set_xlim(-0.7, 0.7)
    ax.set_ylim(-0.7, 0.7)
    ax.set_zlim(-0.95, 0.95)
    ax.set_axis_off()
    ax.set_title('Pyrochlore AFM Tetrahedron', fontsize=12,
                 fontweight='bold', color='#004983')
    ax.view_init(elev=25, azim=135)

    fig.tight_layout()
    fig.savefig(PIC_DIR / "pyrochlore_frust.png", dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {PIC_DIR / 'pyrochlore_frust.png'}")


def ea_random():
    """3D Edwards-Anderson model with random ±J couplings."""
    fig, ax = plt.subplots(figsize=(5, 4.5))

    L = 4
    np.random.seed(42)
    bonds_h = np.random.choice([-1, 1], size=(L, L))
    bonds_v = np.random.choice([-1, 1], size=(L, L))

    # Draw lattice sites
    for i in range(L):
        for j in range(L):
            ax.plot(i, j, 'o', markersize=10, color='#004983',
                    markeredgecolor='k', markeredgewidth=0.5, zorder=5)

    # Draw horizontal bonds with color
    for i in range(L):
        for j in range(L):
            jj = (j + 1) % L
            c = '#CC0000' if bonds_h[i, j] > 0 else '#0066CC'
            lw = 1.8 if bonds_h[i, j] > 0 else 1.2
            ls = '-' if bonds_h[i, j] > 0 else '--'
            label_h = '$+J$' if bonds_h[i, j] > 0 else '$-J$'
            ax.plot([j, jj], [i, i], color=c, linewidth=lw, linestyle=ls,
                    zorder=3, alpha=0.8)

    # Draw vertical bonds with color
    for i in range(L):
        for j in range(L):
            ii = (i + 1) % L
            c = '#CC0000' if bonds_v[i, j] > 0 else '#0066CC'
            lw = 1.8 if bonds_v[i, j] > 0 else 1.2
            ls = '-' if bonds_v[i, j] > 0 else '--'
            ax.plot([j, j], [i, ii], color=c, linewidth=lw, linestyle=ls,
                    zorder=3, alpha=0.8)

    # Legend
    ferro = mpatches.Patch(color='#CC0000', label='$J_{ij} = +J$ (Ferro)')
    afro = mpatches.Patch(color='#0066CC', label='$J_{ij} = -J$ (AF)')
    ax.legend(handles=[ferro, afro], loc='upper right', fontsize=8,
              framealpha=0.9)

    ax.set_xlim(-0.5, L - 0.5)
    ax.set_ylim(-0.5, L - 0.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('2D EA Model (Random $\\pm J$)', fontsize=12,
                 fontweight='bold', color='#004983')

    fig.tight_layout()
    fig.savefig(PIC_DIR / "ea_random.png", dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {PIC_DIR / 'ea_random.png'}")


def rugged_ffo_co():
    """1D rugged landscape: FFO-CO discovers progressively deeper minima."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)

    x = np.linspace(0, 10, 1000)

    # Landscape with three progressively deeper wells
    def E(x):
        return (0.08 * (x - 5)**2
                - 1.5 * np.exp(-3 * (x - 2.5)**2)
                - 2.2 * np.exp(-3 * (x - 5.0)**2)
                - 3.0 * np.exp(-2 * (x - 7.5)**2)
                + 0.15 * np.sin(5 * x))

    y = E(x)

    # --- Panel (a): Small h → local descent to shallowest well ---
    ax = axes[0]
    ax.plot(x, y, 'k-', linewidth=1.5)
    ax.fill_between(x, y, y.min() - 0.5, alpha=0.1, color='#004983')

    x_start = 1.0
    x_local1 = 2.5
    E1 = E(x_local1)

    ax.plot(x_start, E(x_start), 'o', color='#004983', markersize=10, zorder=5,
            label='Start $x_0$')
    ax.plot(x_local1, E1, 's', color='red', markersize=10, zorder=5,
            label='Local min $x_1$')
    ax.annotate('', xy=(x_local1, E1 + 0.05),
                xytext=(x_start, E(x_start) - 0.05),
                arrowprops=dict(arrowstyle='->', color='#004983', lw=2))

    # Current best level set
    ax.axhline(y=E1, color='red', linestyle='--', alpha=0.6, linewidth=1.2)
    ax.text(0.3, E1 + 0.12, '$L_{x_1}$ (best so far)', fontsize=9, color='red')

    ax.set_title('(a) Step 1: Small $h$, local descent', fontsize=10,
                 fontweight='bold', color='#004983')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$E(x)$')
    ax.set_xlim(0, 10)

    # --- Panel (b): Increase h → CO discovers deeper well ---
    ax = axes[1]
    ax.plot(x, y, 'k-', linewidth=1.5)
    ax.fill_between(x, y, y.min() - 0.5, alpha=0.1, color='#004983')

    # Previous best (faded)
    ax.axhline(y=E1, color='red', linestyle='--', alpha=0.3, linewidth=1)
    ax.text(0.3, E1 + 0.12, '$L_{x_1}$', fontsize=9, color='red', alpha=0.5)
    ax.plot(x_local1, E1, 's', color='red', markersize=8, zorder=5, alpha=0.5)

    # CO probe finds deeper well at x≈5.0
    x_probe2 = 5.0
    E2 = E(x_probe2)
    ax.annotate('', xy=(x_probe2, E2), xytext=(x_local1, E1),
                arrowprops=dict(arrowstyle='<->', color='green', lw=2.5))
    ax.text(3.5, (E1 + E2) / 2 + 0.2, 'Large $h$\nCO: $x_2$ is better!',
            fontsize=9, color='green', ha='center', fontweight='bold')

    ax.plot(x_probe2, E2, '*', color='green', markersize=15, zorder=5,
            label='$x_2$ (deeper well)')

    # New lower level set — bound drops
    ax.axhline(y=E2, color='green', linestyle='--', alpha=0.6, linewidth=1.2)
    ax.text(0.3, E2 + 0.12, '$L_{x_2}$ (bound drops!)', fontsize=9, color='green')

    ax.set_title('(b) Step 2: Increase $h$, CO finds deeper well', fontsize=10,
                 fontweight='bold', color='#004983')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlabel('$x$')
    ax.set_xlim(0, 10)

    # --- Panel (c): Optimize, probe again → find deepest well ---
    ax = axes[2]
    ax.plot(x, y, 'k-', linewidth=1.5)
    ax.fill_between(x, y, y.min() - 0.5, alpha=0.1, color='#004983')

    # Previous level sets (faded)
    ax.axhline(y=E1, color='red', linestyle='--', alpha=0.15, linewidth=1)
    ax.axhline(y=E2, color='green', linestyle='--', alpha=0.3, linewidth=1)
    ax.plot(x_probe2, E2, 's', color='green', markersize=8, zorder=5, alpha=0.5)

    # CO probe finds even deeper well at x≈7.5
    x_probe3 = 7.5
    E3 = E(x_probe3)
    ax.annotate('', xy=(x_probe3, E3), xytext=(x_probe2, E2),
                arrowprops=dict(arrowstyle='<->', color='purple', lw=2.5))
    ax.text(6.3, (E2 + E3) / 2 + 0.2, 'CO: $x_3$\nis still\nbetter!',
            fontsize=9, color='purple', ha='center', fontweight='bold')

    ax.plot(x_probe3, E3, '*', color='purple', markersize=15, zorder=5,
            label='$x_3$ (deepest well)')

    # Final level set — bound drops again
    ax.axhline(y=E3, color='purple', linestyle='--', alpha=0.6, linewidth=1.2)
    ax.text(0.3, E3 + 0.12, '$L_{x_3}$ (bound drops!)', fontsize=9, color='purple')

    # Show overall bound progression
    ax.annotate('Best-so-far\nkeeps dropping',
                xy=(9.2, E3), xytext=(9.2, E1),
                arrowprops=dict(arrowstyle='->', color='#004983', lw=2.5),
                fontsize=9, color='#004983', fontweight='bold', ha='center')

    ax.set_title('(c) Step 3: Repeat, bound drops again', fontsize=10,
                 fontweight='bold', color='#004983')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_xlabel('$x$')
    ax.set_xlim(0, 10)

    fig.suptitle('FFO-CO Multi-Scale Strategy: Progressively Deeper Minima',
                 fontsize=13, fontweight='bold', color='#002F5F', y=1.02)
    fig.tight_layout()
    fig.savefig(PIC_DIR / "rugged_ffo_co.png", dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {PIC_DIR / 'rugged_ffo_co.png'}")


if __name__ == "__main__":
    kagome_afm()
    pyrochlore_frust()
    ea_random()
    rugged_ffo_co()
    print("All figures generated.")
