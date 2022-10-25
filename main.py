from data_structure import *
from material_structure import *
from GlobalOptimization import *
import matplotlib.pyplot as plt
import numpy as np
from numba import *
from scipy import interpolate
from scipy.fft import fft, fftfreq,fftshift, rfft, irfft
from scipy.interpolate import UnivariateSpline
from scipy import signal
import time

class Logger():
    stdout = sys.stdout
    messages = []

    def start(self):
        sys.stdout = self
    def stop(self):
        sys.stdout = self.stdout
    def write(self, text):
        self.messages.append(text)


if __name__ == '__main__':
    log = Logger()
    log.start()
    print("sys")
    print("bye")
    log.stop()

    print(log.messages)
    print('hello')