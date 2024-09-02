# __init__.py file in the UTILS module

import os
import sys

# Get the current working directory as the root directory
CURRENT_PATH = os.getcwd()
ROOT_PATH = os.path.abspath(CURRENT_PATH)

# # Check if ROOT_PATH ends with '/GO-RXR/GO-RXR' and remove the extra 'GO-RXR' if it does
# if ROOT_PATH.endswith('/GO-RXR/GO-RXR'):
#     ROOT_PATH = os.path.dirname(ROOT_PATH)

while os.path.basename(ROOT_PATH) != 'GO-RXR':
    ROOT_PATH = os.path.dirname(ROOT_PATH)
    # Check if we have reached the root of the filesystem
    if ROOT_PATH == os.path.dirname(ROOT_PATH):
        print("Root directory not found")
        break

DATA_PATH = os.path.join(ROOT_PATH, 'DATA')
TESTS_PATH = os.path.join(ROOT_PATH, 'TESTS')

# Append the root directory to the Python path
sys.path.append(ROOT_PATH)