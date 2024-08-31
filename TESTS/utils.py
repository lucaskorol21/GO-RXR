
# This file contains utility functions for the tests.
import os

def get_test_data_path(filename):
    """
    Returns the full path to a test data file, adjusting for different environments.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    # Adjust the base directory if needed
    test_data_dir = os.path.join(base_dir, 'test_data')
    
    return os.path.join(test_data_dir, filename)






if __name__ == "__main__":
    
    pass