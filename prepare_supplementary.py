"""Rebuild the small supplementary atom tables from hash-verified raw PDBs.

Run prepare_inputs.py first to fetch the sources. Offline analyses use the
bundled tables; this script is only needed to audit parsing or refresh inputs.
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def atoms(path, include):
    selected = {}
    for line in path.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line[:6] not in ('ATOM  ', 'HETATM') or not include(line):
            continue
        element = line[76:78].strip()
        if element in ('H','D'):
            continue
        row = dict(atom=line[12:16].strip(), resname=line[17:20].strip(),
                   chain=line[21], resseq=int(line[22:26]), icode=line[26].strip(),
                   element=element, altloc=line[16].strip(), occupancy=float(line[54:60]),
                   b=float(line[60:66]),x=float(line[30:38]),y=float(line[38:46]),z=float(line[46:54]))
        key = tuple(row[k] for k in ('chain','resseq','icode','atom'))
        rank = (bool(row['altloc']),-row['occupancy'],row['altloc'])
        if key not in selected or rank < selected[key][0]:
            selected[key] = (rank,row)
    return [v[1] for v in selected.values()]


def main():
    target = ROOT/'data/supplementary'; target.mkdir(exist_ok=True)
    records = json.loads((ROOT/'data/provenance.json').read_text())
    for pdb in ('1aqb','1qab'):
        path = ROOT/'data/raw'/f'{pdb}.pdb'
        reference = next(r for r in records if r['pdb_id'] == pdb)
        if hashlib.sha256(path.read_bytes()).hexdigest() != reference['sha256']:
            raise ValueError(f'Raw source changed: {pdb}')
    selections = {
        '1qab_protein_atoms.csv': ('1qab', lambda s:s.startswith('ATOM  ')),
        '1aqb_retinol_atoms.csv': ('1aqb', lambda s:s[17:20]=='RTL'),
        '1aqb_loop_atoms.csv': ('1aqb', lambda s:s.startswith('ATOM  ') and s[21]=='A'
                               and (63<=int(s[22:26])<=68 or 94<=int(s[22:26])<=99)),
    }
    for name,(pdb,include) in selections.items():
        rows = atoms(ROOT/'data/raw'/f'{pdb}.pdb',include)
        with (target/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        print(name,len(rows))


if __name__=='__main__':
    main()
