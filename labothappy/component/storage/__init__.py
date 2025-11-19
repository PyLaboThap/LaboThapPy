
import os
import sys

def find_project_root(starting_dir, markers=('connector', 'component')):
    current_dir = starting_dir
    while True:
        if all(os.path.isdir(os.path.join(current_dir, marker)) for marker in markers):
            return current_dir
        parent = os.path.abspath(os.path.join(current_dir, os.pardir))
        if parent == current_dir:
            return None  # Reached root of filesystem
        current_dir = parent

# Normalize paths for reliable comparison
def normalize_path(p):
    return os.path.normcase(os.path.abspath(p))

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = find_project_root(current_dir)

if project_root:
    normalized_root = normalize_path(project_root)
    normalized_sys_path = list(map(normalize_path, sys.path))
    if normalized_root not in normalized_sys_path:
        sys.path.insert(0, project_root)
else:
    raise RuntimeError("Project root not found — expected 'connector' and 'component' directories.")
