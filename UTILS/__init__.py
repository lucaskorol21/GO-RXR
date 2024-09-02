# __init__.py file in the UTILS module

import os
import sys

# Get the root directory of the project
ROOT_DIR = os.getenv("ROOT_DIR", os.getcwd())

# Debugging output to see the initial ROOT_DIR value
print(f"Initial ROOT_DIR: {ROOT_DIR}")

# Adjust ROOT_DIR if necessary
while os.path.basename(ROOT_DIR) != 'GO-RXR':
    ROOT_DIR = os.path.dirname(ROOT_DIR)
    if ROOT_DIR == os.path.dirname(ROOT_DIR):
        print("Root directory not found")
        break

# Debugging output to see the final ROOT_DIR value
print(f"Final ROOT_DIR: {ROOT_DIR}")

DATA_DIR = os.path.join(ROOT_DIR, 'DATA')

# Debugging output to see the DATA_DIR value
print(f"DATA_DIR: {DATA_DIR}")

# Append the root directory to the Python path
sys.path.append(ROOT_DIR)