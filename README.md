<h1 align="center">
  <a href="https://imgur.com/iM3mg3t"><img src="https://i.imgur.com/iM3mg3tm.jpg" title="source: imgur.com" /></a>
</h1>

<h1 align="center">
  <br>
  Global Optimization of Resonant X-ray Reflectometry
  <br>
</h1>

<h4 align="center">A data analysis application for material scientists.</h4>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#how-to-use">how-to-use</a> •
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

## how-to-use

The current version of GO-RXR only works for python 3.7 as this is required by the PythonReflectivity code. I would highly suggest downloading the stand-alone applications instead.

To clone the github repostitory the following line of code can be input into the command prompt window:
```bash
# Clone this repository
$ git clone https://github.com/lucaskorol21/MaterialReflection.git
```
Personally, I prefer to clone the repository in the 'code' tab in github.

The dependencies and libraries that need to be included are:
 - Python 3.7
 - Cython 0.29.24
 - PyQt5 5.15.7
 - h5py 2.9.0
 - matplotlib 3.4.3
 - numba 0.55.2
 - numpy 1.21.4
 - scipy 1.7.1
 - steuptools 62.3.2

After setting up your python interpreter run the following code in the python terminal.

```bash
# Installing Pythonreflectivity cython file
$ python setup.py install
```
Once all of this is setup the python file GUI_GO.py can be run. This will open up the GO-RXR workspace and you can begin analyzing your data. It is important to note that GO-RXR can only load in project workspaces with the '.h5' extensions. The documentation can be found in the GO-RXR documentation in this repository, or by clicking the help tab in the application.


## Download
Currently, there are a few updates that still need to be done before a stand-alone application can be downloaded. The application for GO-RXR is planned to be launched in July of 2023 and will be made accessible throught the link below.

- [Download](https://research-groups.usask.ca/qmax/people.php)


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


This software is copyrighted by Lucas Korol. All rights reserved.


Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Publications

There are currently no papers published using GO-RXR, but there are several papers in progress using GO-RXR.

## Bugs
* Data fitting will run the script regardless if script option is selected. It is suggested to comment out every line of the script for now until this issues is resolved.
* Issues with providing an element variation identifier with the same name as a previously defined element. If possible always provide a different element variation idenitifer name that is different than any of the element names provided in any layer.

---

If issues feel free to contact me at lsk601@usask.ca.


