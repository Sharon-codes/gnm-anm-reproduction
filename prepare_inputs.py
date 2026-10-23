"""Re-fetch original PDBs with SHA256 checks and export selected CA tables.

Already bundled CA tables allow the numerical benchmark to run offline.
"""
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from enm import read_pdb_ca

ROOT = Path(__file__).resolve().parent


def main():
    config = json.loads((ROOT / 'config.json').read_text())
    provenance = json.loads((ROOT / 'data/provenance.json').read_text())
    raw = ROOT / 'data/raw'
    raw.mkdir(exist_ok=True)
    ca = ROOT / 'data/ca'
    ca.mkdir(exist_ok=True)
    for record in provenance:
        path = raw / (record['pdb_id'] + '.pdb')
        if not path.exists():
            with urlopen(record['url'], timeout=60) as response:
                path.write_bytes(response.read())
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != record['sha256']:
            raise ValueError(f'Source changed for {path.name}; expected {record["sha256"]}, got {digest}')
    for case in config['gnm_cases'] + [config['anm_case']]:
        rows = read_pdb_ca(raw / (case['pdb'] + '.pdb'), case['network_chains'])
        with (ca / (case['id'] + '.csv')).open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(case['id'], len(rows))


if __name__ == '__main__':
    main()
