
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
  <a href="#license">License</a>
</p>

<a href="https://imgur.com/2ucv2rp"><img src="https://i.imgur.com/2ucv2rp.png" title="source: imgur.com" /></a>

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

Currently, there are a few updates that need to be done and a stand-alone application of GO-RXR is not currently available. The application is planned to be launched in July of 2023.


## Credits

This software uses the following open source packages:
- [Pythonreflectivity](https://github.com/malaclypseII/PyXMRTool.git)

## License

University of Saskatchewan


---

If issues feel free to contact me at lsk601@usask.ca.


