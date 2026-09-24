# LaboThApPy: An open-source tool for advanced thermodynamic cycle simulation, designed for researchers and engineers.
#
# Copyright (C) <2025> <Université catholique de Louvain (UCLouvain), Belgique
#                       Université de Liège (ULiège), Belgique
#                       Université de Mons (UMONS), Belgique>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# List of the contributors to the development of LaboThApPy: see AUTHORS file.
# Description and complete License: see NOTICE & LICENSE files.

import sys
import os

def find_project_root(current_dir, marker_files=None):
    """
    Dynamically finds the root directory of the project based on the presence of marker files.
    Marker files can be a file like '.git', 'requirements.txt', or any specific folder structure.
    
    Parameters:
        current_dir (str): The starting directory (usually the current directory).
        marker_files (list): A list of files or directories that help identify the project root.
                             Example: ['.git', 'requirements.txt']
    
    Returns:
        str: The root directory path if found, else None.
    """
    if marker_files is None:
        marker_files = ['.git', 'requirements.txt']  # Default marker files

    current_dir = os.path.abspath(current_dir)

    while current_dir != os.path.dirname(current_dir):  # Stop when at the root of the file system
        if any(os.path.exists(os.path.join(current_dir, marker)) for marker in marker_files):
            return current_dir
        current_dir = os.path.dirname(current_dir)  # Move up one level

    return None

# Find the project root dynamically
project_root = find_project_root(os.path.dirname(__file__))

if project_root:
    sys.path.append(project_root)  # Add the project root to the sys.path
else:
    print("Project root not found!")