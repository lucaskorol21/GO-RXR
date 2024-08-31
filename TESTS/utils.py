
# This file contains utility functions for the tests.
import os

def find_root_dir(target_dir_name):
    """
    Recursively finds the root directory of the project by looking for a directory
    with the specified name ('target_dir_name') starting from the current working directory.

    Args:
    - target_dir_name (str): The name of the directory to find.

    Returns:
    - str: The absolute path to the root directory if found, otherwise None.
    """
    current_dir = os.getcwd()

    while True:
        if os.path.basename(current_dir) == target_dir_name:
            return current_dir
        parent_dir = os.path.dirname(current_dir)
        # Check if we have reached the root of the filesystem
        if parent_dir == current_dir:
            print(f"Root directory '{target_dir_name}' not found")
            return None
        current_dir = parent_dir






if __name__ == "__main__":
    
    pass