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