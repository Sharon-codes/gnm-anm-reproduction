"""Verify evidence integrity and completeness, not publication equivalence."""
import csv
import hashlib
import json
import struct
import zlib
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))


def check_png(path):
    data=path.read_bytes()
    if data[:8]!=b'\x89PNG\r\n\x1a\n':return False
    offset=8
    while offset<len(data):
        length=struct.unpack('>I',data[offset:offset+4])[0]
        kind=data[offset+4:offset+8];body=data[offset+8:offset+8+length]
        crc=struct.unpack('>I',data[offset+8+length:offset+12+length])[0]
        if zlib.crc32(kind+body)&0xffffffff!=crc:return False
        offset+=12+length
        if kind==b'IEND':return offset==len(data)
    return False


def main():
    checks=[]
    def record(name,passed,detail):
        checks.append(dict(check=name,passed=bool(passed),detail=detail))
    main_manifest=json.loads((RESULTS/'run_manifest.json').read_text())
    extended=json.loads((RESULTS/'extended/run_manifest.json').read_text())
    for label,manifest in [('main',main_manifest),('extended',extended)]:
        for name,expected in manifest['source_sha256'].items():
            record(label+' source '+name,digest(ROOT/name)==expected,expected)
        for name,expected in manifest['input_sha256'].items():
            path=ROOT/name if label=='extended' else ROOT/'data/ca'/name
            record(label+' input '+name,digest(path)==expected,expected)
    validation=rows(RESULTS/'prody_validation.csv')
    record('independent ProDy comparisons',len(validation)==18 and all(r['passed']=='True' for r in validation),'18 expected cases; matched zero thresholds')
    log=(RESULTS/'test_log.txt').read_text()
    record('scientific test log','Ran 13 tests' in log and '\nOK\n' in log,'13 checks; see test_log.txt')
    primary=rows(RESULTS/'1aqb_A_anm_profile.csv')
    with np.load(RESULTS/'extended/1aqb_anm_covariance.npz',allow_pickle=False) as archive:
        c=archive['unit_covariance'];m=float(archive['mean_B_calibration_multiplier'])
        record('full ANM covariance shape',c.shape==(525,525),str(c.shape))
        record('full ANM covariance finite/symmetric',np.isfinite(c).all() and np.allclose(c,c.T), 'float64 unit covariance; thermal scale stored separately')
        b=np.diag(c).reshape(-1,3).sum(axis=1)*(8*np.pi**2/3)*m
        record('saved covariance reproduces primary B factors',np.allclose(b,[float(r['scaled_model_B_A2']) for r in primary],atol=1e-10),'same 13 A primary network')
    required=['anm_figure2_displacements.png','anm_figure3_force_constants.png',
              'anm_figure4_spectra.png','anm_figure7_anisotropy.png',
              'anm_figure8_global_mode.png','anm_structural_panels.png',
              'anm_figure10a_conformations.png','anm_figure10b_experimental_complex.png',
              'gnm_figure4_correlations.png','gnm_figures5_6_connectivity.png']
    for name in required:
        p=RESULTS/'extended'/name
        record('PNG integrity '+name,p.exists() and check_png(p),'all PNG chunk CRCs and final IEND checked')
    spectrum=rows(RESULTS/'extended/anm_all_frequencies.csv')
    record('complete spectra',len(spectrum)==4*519,'519 internal modes at each of 10/13/20/70 A')
    vectors=rows(RESULTS/'extended/anm_first12_mode_vectors.csv')
    record('complete exported mode vectors',len(vectors)==12*175,'12 modes, 175 residues, 3 components')
    peaks=rows(RESULTS/'anm_published_figure5_landmarks.csv')
    record('corrected theory landmarks',len(peaks)==6 and not any(r['source_pixel_x']=='658' for r in peaks),'experimental peak at pixel (658,160) excluded')
    force=rows(RESULTS/'extended/anm_figure3_digitized_comparison.csv')
    record('spring-constant published points',len(force)==12,'5 GNM and 7 ANM approximate source points')
    # This table records scientific failures as failures, not integrity failures.
    status=[dict(comparison='GNM Table 1 normalization',status='unreconciled',evidence='results/gnm_benchmark.csv'),
            dict(comparison='GNM Figure 4 marked correlation regions',status='qualitative agreement, sparse location comparison',evidence='results/extended/gnm_figure4_NMA_region_comparison.csv'),
            dict(comparison='ANM Figure 3 GNM scale',status='within approximate digitization error, assuming 300 K',evidence='results/extended/anm_figure3_digitized_comparison.csv'),
            dict(comparison='ANM Figure 3 standard ANM scale',status='disagrees; threefold alternative agrees under stated assumptions',evidence='results/extended/anm_figure3_digitized_comparison.csv'),
            dict(comparison='ANM Figure 5 theory amplitudes',status='unreconciled; 5 of 6 landmarks differ beyond 1 A squared',evidence='results/anm_published_figure5_landmarks.csv'),
            dict(comparison='ANM Figure 10a mode ratio',status='unreconciled: 0.8207 versus approximately 1.2794',evidence='results/extended/anm_mode_amplitude_comparison.csv')]
    with (RESULTS/'published_result_status.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(status[0]));writer.writeheader();writer.writerows(status)
    outputs={str(p.relative_to(ROOT)):digest(p) for p in sorted(RESULTS.rglob('*'))
             if p.is_file() and p.name!='evidence_verification.json'}
    result=dict(passed=all(c['passed'] for c in checks),integrity_checks=checks,
                verified_output_sha256=outputs,
                publication_equivalence='Not established; see published_result_status.csv and docs/figure_coverage.md')
    (RESULTS/'evidence_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} evidence integrity checks passed")
    for c in checks:
        if not c['passed']:print('FAILED:',c['check'])
    if not result['passed']:raise SystemExit(1)


if __name__=='__main__':main()
