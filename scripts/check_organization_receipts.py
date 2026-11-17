"""Replay paid-organization artifact integrity without recertifying the lattice."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    names=['organization.json','organization-verification.json','organization-continuum-check.json',
           'organization-force-ablation.json','organization-certificate-controls.json','organization-review.json',
           'organization-doc-verification.json']
    checked=0
    for name in names:
        record=json.loads((ROOT/'results'/name).read_bytes())
        values=record.get('source_sha256',record.get('sha256',{}))
        for file,expected in values.items():
            assert sha(ROOT/file)==expected,(name,file)
            checked+=1
    data=json.loads((ROOT/'results/organization.json').read_bytes())
    assert data['verification_sha256']==sha(ROOT/'results/organization-verification.json')
    docs=json.loads((ROOT/'results/organization-doc-verification.json').read_bytes())
    assert docs['status']=='PASS' and docs['tests']['return_code']==0
    proof=json.loads((ROOT/'results/organization-verification.json').read_bytes())
    assert proof['certificates'][-1]['spectator_minus_best_fixed_interval'][0]>0
    assert all(row['passed_dkw_band'] for row in proof['paths']['checks'])
    record={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
        'hashes_replayed':checked,'receipts':names,'tests_passed':docs['tests']['count'],
        'certificate_grid':proof['certificates'][-1]['n'],'continuum_error_enclosed':False,
        'scope':'Source/artifact integrity and retained acceptance replay, not a new scientific proof review.'}
    (ROOT/'results/organization-integrity.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
