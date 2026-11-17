"""Build a portable dossier archive, excluding caches and third-party papers."""
from pathlib import Path
import hashlib
import json
import zipfile
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def included(path):
    relative = path.relative_to(ROOT)
    return (path.is_file() and path.name not in {'ARTIFACT_MANIFEST.json', 'MEMORY.md'}
            and not any(part in {'.git', '.venv', '__pycache__'} or part.endswith('.egg-info') for part in relative.parts)
            and path.suffix != '.pyc'
            and not (relative.parts[:2] == ('research', 'sources') and path.suffix in {'.pdf', '.txt'}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tracked-only', action='store_true',
                        help='Package only staged/tracked files, excluding concurrent unpublished work.')
    args = parser.parse_args()
    if args.tracked_only:
        names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode('utf-8').split('\0')
        paths = sorted(ROOT / name for name in names if name and included(ROOT / name))
    else:
        paths = sorted(path for path in ROOT.rglob('*') if included(path))
    manifest = [{'path': path.relative_to(ROOT).as_posix(), 'bytes': path.stat().st_size,
                 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths]
    manifest_path = ROOT / 'ARTIFACT_MANIFEST.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    archive_path = ROOT.parent / 'tsimerman-research-dossier.zip'
    with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in [*paths, manifest_path]:
            archive.write(path, Path('tsimerman-research-dossier') / path.relative_to(ROOT))
        archive.writestr('START-HERE.txt',
                         'Extract the full archive before opening the HTML.\n'
                         'Start: tsimerman-research-dossier/docs/DEVELOPMENT-JOURNEY.html\n'
                         'Science: tsimerman-research-dossier/docs/SCIENTIFIC-DOSSIER.html\n'
                         'Paper: tsimerman-research-dossier/results/research-manuscript.pdf\n'
                         'Candidate research results; priority and independent review remain open.\n')
    digest = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    archive_path.with_suffix('.zip.sha256').write_text(f'{digest}  {archive_path.name}\n', encoding='utf-8')
    with zipfile.ZipFile(archive_path) as archive:
        assert archive.testzip() is None
        for entry in manifest:
            data = archive.read('tsimerman-research-dossier/' + entry['path'])
            assert hashlib.sha256(data).hexdigest() == entry['sha256']
    print(json.dumps({'archive': str(archive_path), 'files_in_manifest': len(paths),
                      'bytes': archive_path.stat().st_size, 'sha256': digest}, indent=2))


if __name__ == '__main__':
    main()
