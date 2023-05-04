
<h1 align="center">
  <br>
  Global Optimization of Resonant X-ray Reflectometry
  <br>
</h1>

<h4 align="center">A data analysis application for material scientists.</h4>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#Running-in-Python">How To Use</a> •
  <a href="#download">Download</a> •
  <a href="#credits">Credits</a> •
  <a href="#related">Related</a> •
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

## Running-in-Python

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

- [Electron](http://electron.atom.io/)
- [Node.js](https://nodejs.org/)
- [Marked - a markdown parser](https://github.com/chjj/marked)
- [showdown](http://showdownjs.github.io/showdown/)
- [CodeMirror](http://codemirror.net/)
- Emojis are taken from [here](https://github.com/arvida/emoji-cheat-sheet.com)
- [highlight.js](https://highlightjs.org/)

## Related

[markdownify-web](https://github.com/amitmerchant1990/markdownify-web) - Web version of Markdownify

## Support

<a href="https://www.buymeacoffee.com/5Zn8Xh3l9" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/purple_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>

<p>Or</p> 

<a href="https://www.patreon.com/amitmerchant">
	<img src="https://c5.patreon.com/external/logo/become_a_patron_button@2x.png" width="160">
</a>

## You may also like...

- [Pomolectron](https://github.com/amitmerchant1990/pomolectron) - A pomodoro app
- [Correo](https://github.com/amitmerchant1990/correo) - A menubar/taskbar Gmail App for Windows and macOS

## License

University of Saskatchewan


---

> [amitmerchant.com](https://www.amitmerchant.com) &nbsp;&middot;&nbsp;
> GitHub [@amitmerchant1990](https://github.com/amitmerchant1990) &nbsp;&middot;&nbsp;
> Twitter [@amit_merchant](https://twitter.com/amit_merchant)

