"""Offline reproduction benchmark. Optional independent ProDy validation."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from enm import (coordinates, read_ca_csv, kirchhoff, hessian,
                 modes_and_covariance, fluctuations, scale_b_factors)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def evaluate(case, model, cutoff, tolerance, validate=False):
    rows = read_ca_csv(ROOT / 'data/ca' / (case['id'] + '.csv'))
    xyz = coordinates(rows)
    matrix = kirchhoff(xyz, cutoff) if model == 'GNM' else hessian(xyz, cutoff)
    start = time.perf_counter()
    values, vectors, covariance, zeros = modes_and_covariance(matrix, tolerance)
    seconds = time.perf_counter() - start
    raw_f = fluctuations(covariance, model == 'ANM')
    mask = np.array([r['chain'] == case['evaluation_chain'] for r in rows])
    measured = np.array([r['b'] for r in rows])[mask]
    predicted, multiplier = scale_b_factors(raw_f[mask], measured)
    unit_b = (8*np.pi**2/3) * raw_f[mask]
    result = dict(case=case['id'], model=model, cutoff_A=cutoff,
                  network_ca=len(rows), evaluated_ca=int(mask.sum()),
                  published_network_count=case['published_network_count'],
                  published_evaluation_count=case['published_count'],
                  zero_modes=zeros, expected_rigid_zero_modes=1 if model == 'GNM' else 6,
                  pearson_r=float(np.corrcoef(predicted, measured)[0, 1]),
                  rmse_B_A2=float(np.sqrt(np.mean((predicted-measured)**2))),
                  same_structure_scale_multiplier=multiplier,
                  covariance_solve_wall_seconds=seconds)
    if model == 'GNM':
        reported = case['published_gstar']
        # Under the 1997 paper's equations 3--5, Gamma=g*L, B=const*diag(L+)/g*.
        result.update(published_gstar=reported,
                      fitted_gstar_equations_3_to_5=1/multiplier,
                      fitted_gstar_relative_difference=(1/multiplier-reported)/reported,
                      rmse_using_published_gstar_A2=float(np.sqrt(np.mean((unit_b/reported-measured)**2))),
                      fitted_gstar_trapezoid=float(np.trapezoid(unit_b)/np.trapezoid(measured)))
    validation = None
    if validate:
        import prody
        reference = prody.GNM() if model == 'GNM' else prody.ANM()
        if model == 'GNM':
            reference.buildKirchhoff(xyz, cutoff=cutoff, gamma=1.0)
            reference_matrix = reference.getKirchhoff()
        else:
            reference.buildHessian(xyz, cutoff=cutoff, gamma=1.0)
            reference_matrix = reference.getHessian()
        reference.calcModes(n_modes='all', zeros=False)
        reference_cov = reference.getCovariance()
        scale = max(float(np.max(np.abs(reference_cov))), 1e-12)
        reference_rmsf = prody.calcSqFlucts(reference[:])
        native_relative_cov_error = float(np.max(np.abs(covariance-reference_cov))/scale)
        # ProDy 2.6.1 uses absolute ZERO=1e-6. Match this convention for the
        # implementation comparison, and also retain our native discrepancy.
        # At the deliberately underconstrained 7 A ANM network they differ.
        from prody.utilities.eigtools import ZERO
        _, _, matched_cov, matched_zeros = modes_and_covariance(matrix, float(ZERO))
        matched_f = fluctuations(matched_cov, model == 'ANM')
        relative_cov_error = float(np.max(np.abs(matched_cov-reference_cov))/scale)
        relative_f_error = float(np.max(np.abs(matched_f-reference_rmsf))/max(np.max(np.abs(reference_rmsf)), 1e-12))
        validation = dict(case=case['id'], model=model, cutoff_A=cutoff,
                          prody_version=prody.__version__,
                          matrix_max_absolute_error=float(np.max(np.abs(matrix-reference_matrix))),
                          native_zero_atol=tolerance, reference_zero_atol=float(ZERO),
                          native_zero_count=zeros, matched_zero_count=matched_zeros,
                          native_covariance_max_relative_error=native_relative_cov_error,
                          covariance_max_relative_error=relative_cov_error,
                          fluctuations_max_relative_error=relative_f_error,
                          passed=bool(np.allclose(matrix, reference_matrix, rtol=1e-12, atol=1e-12)
                                      and relative_cov_error < 1e-8 and relative_f_error < 1e-8))
    evaluated_rows = [r for r, use in zip(rows, mask) if use]
    profile = [dict(ordinal=i+1, chain=r['chain'], resseq=r['resseq'], icode=r['icode'],
                    resname=r['resname'], observed_B_A2=float(obs), scaled_model_B_A2=float(pred),
                    unit_prefactor_ms_fluctuation=float(f))
               for i, (r, obs, pred, f) in enumerate(zip(evaluated_rows, measured, predicted, raw_f[mask]))]
    return result, profile, validation, (values, vectors, covariance)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validate-prody', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    config = json.loads((ROOT / 'config.json').read_text())
    tolerance = config['zero_eigenvalue_absolute_tolerance']
    results, validations, audits = [], [], []
    fig, axes = plt.subplots(4, 3, figsize=(14, 13), constrained_layout=True)
    for ax, case in zip(axes.flat, config['gnm_cases']):
        result, profile, val, details = evaluate(case, 'GNM', config['gnm_cutoff_angstrom'], tolerance, args.validate_prody)
        results.append(result)
        if val: validations.append(val)
        write_csv(OUT / (case['id']+'_gnm_profile.csv'), profile)
        audits.append({k: result[k] for k in ('case','network_ca','evaluated_ca','published_network_count','published_evaluation_count')})
        ax.plot([r['ordinal'] for r in profile], [r['observed_B_A2'] for r in profile], color='#777777', lw=1, label='Crystal B')
        ax.plot([r['ordinal'] for r in profile], [r['scaled_model_B_A2'] for r in profile], color='#157a82', lw=1, label='GNM; scale fitted')
        ax.set(title=f"{case['id']}  r={result['pearson_r']:.3f}", xlabel='Observed C-alpha index', ylabel=r'B ($\AA^2$)')
        print(json.dumps(result), flush=True)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('GNM: original 12 structures; current resolved C-alpha atoms\n7 A cutoff; each curve scaled using its own experimental B factors', fontsize=14)
    fig.savefig(OUT / 'gnm_profiles.png', dpi=160)
    fig.savefig(OUT / 'gnm_profiles.svg')
    plt.close(fig)
    write_csv(OUT / 'gnm_benchmark.csv', results)
    anm_results = []
    case = config['anm_case']
    primary = None
    for cutoff in config['anm_diagnostic_cutoffs_angstrom']:
        result, profile, val, details = evaluate(case, 'ANM', cutoff, tolerance, args.validate_prody)
        anm_results.append(result)
        if val: validations.append(val)
        if cutoff == config['anm_primary_cutoff_angstrom']:
            primary = result, profile, details
        print(json.dumps(result), flush=True)
    write_csv(OUT / 'anm_cutoff_diagnostics.csv', anm_results)
    result, profile, (values, vectors, covariance) = primary
    write_csv(OUT / '1aqb_A_anm_profile.csv', profile)
    audits.append({k: result[k] for k in ('case','network_ca','evaluated_ca','published_network_count','published_evaluation_count')})
    write_csv(OUT / 'data_audit.csv', audits)
    pred = np.array([r['scaled_model_B_A2'] for r in profile])
    peaks, _ = find_peaks(pred, prominence=8)
    write_csv(OUT / 'anm_profile_peaks.csv', [dict(residue=profile[i]['resseq'], scaled_B_A2=float(pred[i])) for i in peaks])
    # Seven clearly visible solid-curve local maxima, manually transcribed from
    # Figure 5's original raster. Calibration is fixed; no fit to our results.
    # Approximate figure data are not original author numerical arrays.
    figure_pixels = [(426,445),(537,206),(658,160),(768,394),(869,372),(1087,190),(1238,349)]
    figure_comparison = []
    for px, py in figure_pixels:
        residue = 40 * (px-180.5)/(470.5-180.5)
        bfig = 75 * (741-py)/(741-22)
        nearby = [i for i in peaks if abs(profile[i]['resseq']-residue) <= 3]
        if len(nearby) != 1:
            raise ValueError('Ambiguous figure landmark; do not automatically pick a best fit')
        i = nearby[0]
        figure_comparison.append(dict(figure='ANM Figure 5', source_pixel_x=px, source_pixel_y=py,
                                     figure_residue_approx=residue, figure_B_A2_approx=bfig,
                                     figure_residue_uncertainty=1.0, figure_B_A2_uncertainty=1.0,
                                     reproduced_peak_residue=profile[i]['resseq'],
                                     reproduced_peak_B_A2=float(pred[i]),
                                     B_difference_A2=float(pred[i]-bfig)))
    write_csv(OUT / 'anm_published_figure5_landmarks.csv', figure_comparison)
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    for ax in axes:
        ax.plot([r['resseq'] for r in profile], [r['observed_B_A2'] for r in profile], color='#777777', lw=1.2, label='1AQB crystal B')
        ax.plot([r['resseq'] for r in profile], pred, color='#a03c65', lw=1.2, label='Independent ANM; 13 A; fitted mean')
        ax.errorbar([r['figure_residue_approx'] for r in figure_comparison],
                [r['figure_B_A2_approx'] for r in figure_comparison],
                xerr=1, yerr=1, fmt='o', markersize=4, color='#252525', label='Published theory peak landmarks (approx.)')
        ax.set(ylabel=r'B ($\AA^2$)')
    axes[0].set_title(f"ANM: r={result['pearson_r']:.3f}; full range shows N-terminal discrepancy")
    axes[0].legend(fontsize=7, loc='upper right')
    axes[1].set(ylim=(0,90), xlabel='PDB residue number', title='Zoom to 0-90 A squared for comparison of internal peaks')
    fig.savefig(OUT / 'anm_published_comparison.png', dpi=160)
    fig.savefig(OUT / 'anm_published_comparison.svg')
    plt.close(fig)
    positive = np.flatnonzero(values > tolerance)
    mode_summary = []
    for number, j in enumerate(positive[:5], 1):
        amp = np.linalg.norm(vectors[:, j].reshape(-1,3), axis=1)
        mode_summary.append(dict(nonzero_mode=number, eigenvalue=float(values[j]),
                                 maximum_amplitude_residue=profile[int(amp.argmax())]['resseq'],
                                 top_five_residues=[profile[int(i)]['resseq'] for i in np.argsort(amp)[-5:][::-1]]))
    (OUT / 'anm_slow_mode_audit.json').write_text(json.dumps(mode_summary, indent=2)+'\n')
    if validations:
        write_csv(OUT / 'prody_validation.csv', validations)
    provenance = dict(run_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
                      platform=platform.platform(), command=' '.join(sys.argv),
                      packages={k: importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')},
                      zero_eigenvalue_atol=tolerance, prody_validation_requested=args.validate_prody,
                      source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in [ROOT/'enm.py', ROOT/'benchmark.py', ROOT/'config.json']},
                      input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in sorted((ROOT/'data/ca').glob('*.csv'))})
    (OUT / 'run_manifest.json').write_text(json.dumps(provenance, indent=2)+'\n')
    if any(not r['passed'] for r in validations):
        raise SystemExit('Independent ProDy validation failed; inspect results/prody_validation.csv')
    print('Finished: numerical reproduction performed; exact agreement with paper is not assumed.')


if __name__ == '__main__':
    main()
