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
base_dir = '_build/html' if local else 'gh-pages'
shared_dir = os.path.join(base_dir, 'shared')

entries = []

versions = sorted(os.listdir(base_dir), reverse=True)
if 'main' in versions:
    versions.insert(0, versions.pop(versions.index('main')))

for name in versions:
    path = os.path.join(base_dir, name)
    if os.path.isdir(path) and name not in ('shared', '.git'):
        entries.append({
            'version': name,
            'url': f'../{name}/' if local else f'/pyremap/{name}/'
        })

os.makedirs(shared_dir, exist_ok=True)
with open(os.path.join(shared_dir, 'versions.json'), 'w') as f:
    json.dump(entries, f, indent=2)
