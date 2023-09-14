<h1 align="center">
  <a href="https://imgur.com/3ZJbg84"><img src="https://i.imgur.com/3ZJbg84.jpg" title="source: imgur.com" /></a>
</h1>

<h1 align="center">
  <br>
  Global Optimization of Resonant X-ray Reflectometry
  <br>
</h1>

<h4 align="center">A data analysis application for material scientists.</h4>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#install">how-to-use</a> •
  <a href="#download">Download</a> •
  <a href="#credits">Credits</a> •
  <a href="#license">License</a> •
  <a href="#publications">Publications</a> 
</p>

<a href="https://imgur.com/kPS6E3D"><img src="https://i.imgur.com/kPS6E3D.jpg" title="source: imgur.com" /></a>

## Key Features

* Graphical User Interface
* Sample definition as a compound-profile
* Adaptive Layer Segmentation
* Internal database of form factors
  - allows for selection of form factors in project file
* Magnetism capabilities
* Compatibility with ReMagX
  - able to load in dataset from ReMagX '.all' file type
* Script functionality with built in functions
* Customizable cost function in data fitting
  - chi-square, L1-norm, and L2-norm regularization
  - create unique boundaries for select data scans with weights associated to them
  - shape parameterization using total variation
* Data smoothing
  - this is specifically used for the shape parameterization
* Progress workspace

## install
1. Download [Python 3.7 (332-bit)](https://www.python.org/downloads/release/python-370/)
2. Download [Microsoft C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
3. Clone GO-RXR from the main branch in this repository
4. Install the python libraries by running the setup python file:

```bash
$ python setup.py install
```
5. If the setup file does not work, then the following libraries can be installed:
    * Cython 0.29.24
    * numpy 1.21.4
    * PyQt5 5.15.7
    * h5py 2.9.0
    * matplotlib 3.4.3
    * numba 0.55.2
    * pyqt5 0.12.4
    * scipy 1.7.1
7. Install Python reflectivity by running
```bash
$ python setup_reflectivity.py install
```
9. If you recieve an error check to make sure that Python 3.7 (32-bit), Cython 0.29.24, and numpy are being used by your project environment.
10. Updates to Pythonreflectivity are required before recent version of Python and Cython can be used in the setup of GO-RXR.

## Download
The stand-alone application for [GO-RXR](https://research-groups.usask.ca/qmax/people.php) will be made available in October 2023.

## Credits

This software uses the following open source packages:
- [Pythonreflectivity](https://github.com/malaclypseII/PyXMRTool.git)

Conrtibution made by:
 - Dr. Robert J. Green
 - Dr. Raymond Spiteri
 - [QMax Group](https://research-groups.usask.ca/qmax/)
 - [Numerical Simulations Research Lab](https://simlab.usask.ca/)

GO-RXR would have not been possible without the University of Saskatchewan and the funding provided by the U of S Physics and Engineering Physics Department, the NSERC-CREATE to INSPIRE fellowship, and the NSERC CGS M.

## License
Copyright © 2023 QMaX and Numerical Simulations Lab

GO-RXR has been developed at in the Department of Physics and Engineering Physics at the University of Saskatchewan and is copyrighted by the QMaX and Numerical Simulation Lab. All rights are reserved by the authors, QMaX Research Group, and the Numerical Simulations Research Lab. 

GO-RXR is a a free software: you can redistribute it and/or modify it under the terms an conditions of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

ANY WORK USING THE SOFTWARE OR ANY RESULTS OBTAINED WITH THE HELP OF THIS SOFTWARE HAS TO CITE GO-RXR AND THE AUTHORS PROPERLY.

**Disclaimer of Warranty**

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

**Disclaimer of Liability**

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

END OF TERMS AND CONDITIONS

## Publications

There are currently no papers published using GO-RXR, but there are several papers in progress using GO-RXR.

## Bugs
* Resolved (5/10/2023): Data fitting will run the script regardless if script option is selected. It is suggested to comment out every line of the script for now until this issues is resolved.
* Issues with providing an element variation identifier with the same name as a previously defined element. If possible always provide a different element variation identifier name that is different than any of the element names provided in any layer.

---

If issues feel free to contact me at lsk601@usask.ca.


