"""Additional paper analyses, with explicit boundary and reference conventions."""
import numpy as np
from enm import contacts


def cross_correlation(covariance, anisotropic=False):
    c = np.asarray(covariance)
    if anisotropic:
        n = len(c) // 3
        c = np.einsum('iaja->ij', c.reshape(n, 3, n, 3))
    sd = np.sqrt(np.maximum(np.diag(c), 0))
    return np.divide(c, sd[:, None]*sd[None, :],
                     out=np.zeros_like(c), where=(sd[:, None]*sd[None, :]) > 0)


def chain_matrix(n, boundary='free'):
    """Free Rouse path; literal printed Eq.8; or ring interpretation (1997)."""
    if n < 3:
        raise ValueError('Use at least three nodes')
    k = 2*np.eye(n)-np.eye(n, k=1)-np.eye(n, k=-1)
    if boundary == 'free':
        k[0, 0] = k[-1, -1] = 1
    elif boundary == 'ring':
        k[0, -1] = k[-1, 0] = -1
    elif boundary != 'literal':
        raise ValueError(boundary)
    return k


def principal_frame(xyz):
    """Columns: centered SVD axes, decreasing spatial variance, right handed.

    Fix signs by the largest-magnitude loading; signs do not alter variances.
    The paper does not fully specify origin or axis sign conventions.
    """
    center = np.mean(xyz, axis=0)
    _, singular, vt = np.linalg.svd(np.asarray(xyz)-center, full_matrices=False)
    axes = vt.T
    for j in range(3):
        if axes[np.argmax(np.abs(axes[:, j])), j] < 0:
            axes[:, j] *= -1
    if np.linalg.det(axes) < 0:
        axes[:, -1] *= -1
    return center, axes, singular


def directional_variances(covariance, axes):
    n = len(covariance)//3
    blocks = covariance.reshape(n, 3, n, 3)[np.arange(n), :, np.arange(n), :]
    return np.einsum('ab,nac,cb->nb', axes, blocks, axes)


def equilibrium_matrix(xyz, cutoff):
    """3N x M force/extension operator B; the unit-spring Hessian is B B.T."""
    xyz, delta, d2, adjacency = contacts(xyz, cutoff)
    pairs = np.array(np.where(np.triu(adjacency, 1))).T
    b = np.zeros((3*len(xyz), len(pairs)))
    for column, (i, j) in enumerate(pairs):
        direction = delta[i, j]/np.sqrt(d2[i, j])
        b[3*i:3*i+3, column] = direction
        b[3*j:3*j+3, column] = -direction
    return b, pairs
