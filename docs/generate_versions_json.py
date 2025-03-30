#!/usr/bin/env python
import argparse
import json
import os

# Mode: local or production
parser = argparse.ArgumentParser(
    description='Generate versions.json for pyremap documentation.')
parser.add_argument(
    '--local',
    action='store_true',
    help='Generate versions.json for local build.'
)
args = parser.parse_args()
local = args.local
output_dir = '_build/html/shared'
version_dir = '_build/html' if local else 'gh-pages'

entries = []

versions = sorted(os.listdir(version_dir))
if 'main' in versions:
    versions.insert(0, versions.pop(versions.index('main')))

for name in versions:
    path = os.path.join(version_dir, name)
    if os.path.isdir(path) and name not in ('shared',):
        entries.append({
            'version': name,
            'url': f'../{name}/' if local else f'/{name}/'
        })

os.makedirs(output_dir, exist_ok=True)
with open(os.path.join(output_dir, 'versions.json'), 'w') as f:
    json.dump(entries, f, indent=2)
