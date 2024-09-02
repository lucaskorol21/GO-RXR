# __init__.py file in the UTILS module

import os
import sys

# Get the root directory of the project
ROOT_DIR = os.getenv("ROOT_DIR", os.getcwd())
while os.path.basename(ROOT_DIR) != 'GO-RXR':
    ROOT_DIR = os.path.dirname(ROOT_DIR)
    if ROOT_DIR == os.path.dirname(ROOT_DIR):
        print("Root directory not found")
        break

DATA_DIR = os.path.join(ROOT_DIR, 'DATA')
# Append the root directory to the Python path
sys.path.append(ROOT_DIR)