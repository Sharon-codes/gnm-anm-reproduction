"""Paper-derived, unit-spring Gaussian and anisotropic protein networks.

Covariances omit the global thermal/spring prefactor. They describe native-state
harmonic fluctuations, not molecular dynamics trajectories or folding.
"""
import csv
from pathlib import Path
import numpy as np
from scipy.linalg import eigh


def read_pdb_ca(path, chains):
    """First model, ATOM C-alpha only; one conformer per residue.

    Prefer blank altloc, then highest occupancy, then lexical altloc. HETATM
    calcium ions named CA are deliberately excluded. Keep insertion codes.
    """
    selected = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if not (line.startswith('ATOM  ') and line[12:16].strip() == 'CA'
                and line[21] in chains):
            continue
        row = dict(chain=line[21], resseq=int(line[22:26]), icode=line[26].strip(),
                   resname=line[17:20].strip(), altloc=line[16].strip(),
                   occupancy=float(line[54:60]), b=float(line[60:66]),
                   x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]))
        key = (row['chain'], row['resseq'], row['icode'])
        rank = (row['altloc'] != '', -row['occupancy'], row['altloc'])
        old = selected.get(key)
        if old is None or rank < old[0]:
            selected[key] = (rank, row)
    if not selected:
        raise ValueError(f'No C-alpha atoms for {path}, chains={chains!r}')
    return [item[1] for item in selected.values()]


def read_ca_csv(path):
    with open(path, newline='') as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row['resseq'] = int(row['resseq'])
        for key in ('x', 'y', 'z', 'b', 'occupancy'):
            row[key] = float(row[key])
    return rows


def coordinates(rows):
    return np.array([[r[k] for k in ('x', 'y', 'z')] for r in rows])


def contacts(xyz, cutoff):
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 2:
        raise ValueError('Coordinates must be an N by 3 array, N >= 2')
    if not np.isfinite(xyz).all() or not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('Finite coordinates and positive finite cutoff required')
    delta = xyz[:, None, :] - xyz[None, :, :]
    d2 = np.einsum('ijk,ijk->ij', delta, delta)
    if np.any((d2 == 0) & ~np.eye(len(xyz), dtype=bool)):
        raise ValueError('Coincident nodes are not a valid elastic network')
    adjacency = (d2 <= cutoff ** 2) & ~np.eye(len(xyz), dtype=bool)
    return xyz, delta, d2, adjacency


def kirchhoff(xyz, cutoff=7.0):
    _, _, _, adjacency = contacts(xyz, cutoff)
    return np.diag(adjacency.sum(axis=1)) - adjacency.astype(float)


def hessian(xyz, cutoff=13.0):
    xyz, delta, d2, adjacency = contacts(xyz, cutoff)
    out = np.zeros((3 * len(xyz), 3 * len(xyz)))
    for i, j in zip(*np.where(np.triu(adjacency, 1))):
        block = np.outer(delta[i, j], delta[i, j]) / d2[i, j]
        si, sj = slice(3*i, 3*i+3), slice(3*j, 3*j+3)
        out[si, si] += block
        out[sj, sj] += block
        out[si, sj] -= block
        out[sj, si] -= block
    return out


def modes_and_covariance(matrix, atol=1e-8):
    """Exclude all eigenvalues <= atol; report their number, never hide extras."""
    values, vectors = eigh(matrix, check_finite=True)
    if values.min() < -atol:
        raise ValueError('Network has a significantly negative eigenvalue')
    positive = values > atol
    covariance = (vectors[:, positive] / values[positive]) @ vectors[:, positive].T
    return values, vectors, covariance, int((~positive).sum())


def fluctuations(covariance, anisotropic=False):
    diagonal = np.diag(covariance)
    return diagonal.reshape(-1, 3).sum(axis=1) if anisotropic else diagonal


def scale_b_factors(unit_fluctuations, measured):
    """Match discrete areas, sum(B), with a single multiplicative constant.

    Uses observed B factors for scaling. This is descriptive calibration on the
    same structure, not held-out prediction. Pearson r is scale invariant.
    """
    unit_b = (8 * np.pi**2 / 3) * np.asarray(unit_fluctuations)
    measured = np.asarray(measured)
    if unit_b.sum() <= 0 or measured.sum() <= 0:
        raise ValueError('Positive total fluctuations and B factors required')
    multiplier = measured.sum() / unit_b.sum()
    return unit_b * multiplier, float(multiplier)
