"""Replay retained provenance checks without rerunning the Monte Carlo study."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'results/spatial-repair-delivery-manifest.json'
INTEGRITY=ROOT/'results/spatial-repair-integrity.json'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):return json.loads((ROOT/name).read_bytes())


def check_group(values):
    for relative,expected in values.items():
        assert digest(ROOT/relative)==expected, f'Hash mismatch: {relative}'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest',action='store_true',help='Seal current delivery after validating its receipts.')
    args=parser.parse_args()
    data=read('results/spatial-repair.json')
    check_group(data['source_sha256'])
    encounter=read('results/encounter-feasibility.json')
    check_group(encounter['source_sha256'])
    encounter_check=read('results/encounter-verification.json')
    assert encounter_check['return_code']==0
    assert encounter_check['data_sha256']==digest(ROOT/'results/encounter-feasibility.json')
    verified=read('results/spatial-repair-verification.json')
    assert verified['passed'] is True
    assert verified['input_sha256']==digest(ROOT/'results/spatial-repair.json')
    check_group(verified['source_sha256']);check_group(verified['figure_sha256'])
    ui=read('results/spatial-repair-ui-validation.json')
    assert ui['status']=='PASS'
    check_group(ui['source_and_artifact_sha256'])
    review=read('results/spatial-repair-review-check.json')
    assert review['all_cases_within_familywise_dkw_band'] is True
    assert review['input_sha256']==digest(ROOT/review['input_path'])
    assert review['formal_verification_sha256']==digest(ROOT/review['formal_verification_path'])
    check_group(review['source_sha256'])
    if args.write_manifest:
        paths={ROOT/'README.md',ROOT/'pyproject.toml',ROOT/'scripts/build_research_docs.py',
               ROOT/'scripts/verify_research_docs.py',ROOT/'docs/SPATIAL-REPAIR.html'}
        for pattern in ['src/purposeful/*.py','tests/test_*.py','scripts/*spatial*',
                        'docs/*spatial*','research/spatial-repair*','results/spatial-repair*',
                        'results/figures/spatial-repair*','results/document-previews/spatial-repair*',
                        'docs/assets/katex/**/*','scripts/*stochastic*','docs/*stochastic*',
                        'docs/STOCHASTIC-CALCULUS-PLAN.html','research/stochastic*',
                        'results/encounter*','results/stochastic*','results/figures/encounter*',
                        'results/document-previews/stochastic*']:
            paths.update(ROOT.glob(pattern))
        paths={p for p in paths if p.is_file() and p not in {MANIFEST,INTEGRITY}}
        items=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)}
               for p in sorted(paths)]
        MANIFEST.write_text(json.dumps({'generated_utc':datetime.now(timezone.utc).isoformat(),
            'scope':'Spatial repair delivery, scientific source dependencies and offline equation assets.',
            'exclusions':['Project-private MEMORY.md','.git','caches'],
            'artifacts':items},indent=2)+'\n',encoding='utf-8')
    manifest=json.loads(MANIFEST.read_bytes())
    check_group({item['path']:item['sha256'] for item in manifest['artifacts']})
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
        'manifest_sha256':digest(MANIFEST),'artifacts_checked':len(manifest['artifacts']),
        'formal_tests_passed':verified['summary']['test_count'],
        'formal_continuous_paths':verified['summary']['monte_carlo_paths_total'],
        'supplemental_fresh_seed_paths':review['paths_total'],
        'scope':'Integrity and retained-acceptance replay; not a new simulation or independent proof review.'}
    INTEGRITY.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
