![vscreenml](vscreenml.png)

<!--- These are examples. See https://shields.io for others or to customize this set of shields. You might want to include dependencies, project status and licence info here --->
![GitHub repo size](https://img.shields.io/github/repo-size/karanicolaslab/vscreenml)
![GitHub contributors](https://img.shields.io/github/contributors/karanicolaslab/vscreenml)
![GitHub stars](https://img.shields.io/github/stars/karanicolaslab/vscreenml?style=social)
![GitHub forks](https://img.shields.io/github/forks/karanicolaslab/vscreenml?style=social)

vScreenML is a tool`that allows rescoring of virtual screening hits to prune out false positives.

For more deatils, check out PNAS paper at <a href="https://www.biorxiv.org/content/10.1101/2020.01.10.902411v1.full.pdf">PNAS paper</a>

## Prerequisites

Before you begin, ensure you have met the following requirements:
<!--- These are just example requirements. Add, duplicate or remove as required --->
* You have installed the following packages [OpenEye](https://www.eyesopen.com), [ChemAxon calcx](https://chemaxon.com/academic-license), [ChemAxon structurecheck](https://chemaxon.com/academic-license), [MGLTOOLS](http://mgltools.scripps.edu/downloads), [BINANA](https://sourceforge.net/projects/binana)
* Set up Rosetta to work with vscreenml
    * Request license here https://els.comotion.uw.edu/express_license_technologies/rosetta (academic license is free). Once you obtain the license, you can download Rosetta here https://www.rosettacommons.org/software/license-and-download and do the following
    * In order to use this tool, it's better to build Rosetta from scratch. Once Rosetta repositiory is downloaded, before compiling do the following:
        1. from your terminal change into path/to/Rosetta/main/source/src/apps/pilot/
            ```
           cd path/to/Rosetta/main/source/src/apps/pilot/
           ```
        2. Inside path/to/Rosetta/main/source/src/apps/pilot/ create a new directory "vscreenml_rosetta_extension"
           ```
           mkdir vscreenml_rosetta_extension
           ```
        3. copy vscreenml_rosetta_extension/vscrenml_interface_ddg.cc into path/to/Rosetta/main/source/src/apps/pilot/vscreenml_rosetta_extension
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
          ```python
          scons bin mode=release
          ```
## Installing vScreenML
* Clone this repository

* Install all the required python packages
    ```
    conda create --name myenv --file requirements.txt
    ```
* Install all the above dependencies and update configuration/congif.py file paths pointing to this dependencies.

Once all the dependencies have been installed, vscreenml is ready to use.

## Using vScreenML

To use vScreenML, follow these steps:
For example, in order go get the vscreenml score of "mini_complex_sample.pdb" present in "test" directory
    ```
    python vscreenmlscore.py test sample
    ```
Add run commands and examples you think users will find useful. Provide an options reference for bonus points!

## Contributing to vScreenML
<!--- If your README is long or you have some specific process or steps you want contributors to follow, consider creating a separate CONTRIBUTING.md file--->
To contribute to vScreenML, follow these steps:

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Make your changes and commit them: `git commit -m '<commit_message>'`
4. Push to the original branch: `git push origin vscreenml/<location>`
5. Create the pull request.

Alternatively see the GitHub documentation on [creating a pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request).

## Contact

If you want to contact me you can reach me at <dryusufadeshina@gmail.com>.

## License
<!--- If you're not sure which open license to use see https://choosealicense.com/--->

This project uses the following license: [MIT License](<link>).

