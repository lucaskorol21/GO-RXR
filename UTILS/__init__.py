# __init__.py file in the UTILS module

import os
import sys

# Get the current working directory as the root directory
CURRENT_DIR = os.getcwd()

ROOT_DIR = os.path.abspath(CURRENT_DIR)

while os.path.basename(ROOT_DIR) != 'GO-RXR':
    ROOT_DIR = os.path.dirname(ROOT_DIR)
    # Check if we have reached the root of the filesystem
    if ROOT_DIR == os.path.dirname(ROOT_DIR):
        print("Root directory not found")
        break

# Append the root directory to the Python path
sys.path.append(ROOT_DIR)