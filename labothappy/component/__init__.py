# import os
# import sys

# def find_project_root(starting_dir):
#     # Define markers that identify the project root
#     markers = ['connector', 'component']

#     current_dir = starting_dir
#     while True:
#         # Check if all markers exist in the current directory
#         if all(os.path.isdir(os.path.join(current_dir, marker)) for marker in markers):
#             return current_dir

#         # Move up one directory
#         parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
        
#         # If we have reached the root of the filesystem, stop searching
#         if parent_dir == current_dir:
#             return None
        
#         current_dir = parent_dir

# # Get the absolute path of the directory that contains the current script
# current_dir = os.path.dirname(os.path.abspath(__file__))

# # Find the project root directory
# project_root = find_project_root(current_dir)

# if project_root:
#     # Add the project root to sys.path if it's not already there
#     if project_root not in sys.path:
#         sys.path.insert(0, project_root)
# else:
#     raise RuntimeError("Project root not found. Make sure you have 'connector' and 'component' directories.")
# 

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
