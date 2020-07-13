# vscreenml
![stack Overflow](vscreenml.png)
## Description
This repository contains the implementation of a novel machine learning classifier trained on the Dataset of Congruent Inhibitors and Decoys (D-COID)).  
## Dependencies
### Rosetta (version 3)
Request license here https://els.comotion.uw.edu/express_license_technologies/rosetta (academic license is free). Once you obtain the license, you can download Rosetta here https://www.rosettacommons.org/software/license-and-download .

In order to use this tool, it's better to build Rosetta from scratch. Once Rosetta repositiory is downloaded, before compiling do the following:
1. from your terminal change into path/to/Rosetta/main/source/src/apps/pilot/
```
cd path/to/Rosetta/main/source/src/apps/pilot/
```
2. Inside path/to/Rosetta/main/source/src/apps/pilot/ create a new directory "vscreenml_rosetta_extension"
```
mkdir vscreenml_rosetta_extension
```
3.  copy vscreenml_rosetta_extension/vscrenml_interface_ddg.cc into path/to/Rosetta/main/source/src/apps/pilot/vscreenml_rosetta_extension
```
cp /path/to/vscreenml/vscreenml_rosetta_extension/vscrenml_interface_ddg.cc vscreenml_rosetta_extension
```
4. add the json formated text below to path/to/Rosetta/main/source/src/pilot_apps.src.settings.all
```
"pilot/vscreenml_rosetta_extension" : [ 
                "vscrenml_interface_ddg",
        ],
 ```      
5. Compile Rosetta by changing to the source directory and running:
```
scons bin mode=release
```

### OpenEye OMEGA and Szybki
OpenEye is not a free software, however free academic licence can be obtained here https://www.eyesopen.com/
### Openbabel
Already contained in the requirements.txt file
### ChemAxon calcx, structurecheck
ChemAxon tools are not free, however free academic licence can be obtained here https://chemaxon.com/academic-license
### MGLTOOLS
MGLTOOLS tools is free and can be obtained here http://mgltools.scripps.edu/downloads
### BINANA 
BINANA is freely available and can be download at https://sourceforge.net/projects/binana/


## Installation and Use
Clone this repository

Install all the required python packages
```
conda create --name myenv --file requirements.txt
```
Install all the above dependencies and update configuration/congif.py file paths pointing to this dependencies.

Once all the dependencies have been installed, vscreenml is ready to use.

Usage:

For example, in order go get the vscreenml score of "mini_complex_sample.pdb" present in "test" directory
```python
python vscreenmlscore.py test sample

```
