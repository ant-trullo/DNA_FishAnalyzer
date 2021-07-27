# DNA_FishAnalyzer
This software is written in Python™ 3.8.

- Packages needed for DNA_FishAnalyzer_v0.0 to work:

	numpy         1.20.3;
	pyqtgraph     0.11.1;
	skimage       0.18.1;
	sklearn       0.22.1;
	PyQt5         5.13.1;
	xlsxwriter    1.4.3;
	czifile       2019.7.2;
	cython        0.29.23.

The number refers to the version we have used. DNA_FishAnalyzer_v0.0 can work anyway 
with newer version of the same packages, unless there are changes in syntax.
        Of course all the dependecies of these packages must be fulfilled: you will be 
        required to install matplotlib, PIL and some others. Depending on the installation 
        technique, the Python installer can take care of this directly, or you will have to do this 
        manually. 

DNA_FishAnalyzer_v0.0 was developed on the operative system Linux Ubuntu 20.04 64-bit 
        and tested on Linux Ubuntu 20.04 64-bit, Windows 10 Pro, Sierra OS.

 
- Run DNA_FishAnalyzer_v0.0:
	Before running the software, you have to compile the NucsPileUpUtility.pyx file:
	for this you need to have a C++ compiler (which is gcc for Linux and Mac) and can be
	installed on Windows following these instructions (https://wiki.python.org/moin/WindowsCompilers).
	Windows compiler can be installed even installing the VisualStudio Editor and loading 
	the C++ compiler in it. Once you are ready, donwload the files in the repository and
	put them all in a folder, then open a terminal 
	(or the cmd command prompt for windows users), go into the folder and type:
	
	'python3 setup.py build_ext --inplace'
	than enter. (This operation must be done only once, by the time you get the compiled file
	you don't need to redo it. A compiled file for Linux python 3.8.10 is already in the repository).
	If nothing goes wrong, type again
	
	'python3 SpotsDistanceGUI_v0_0.py'
	and press enter. The graphical user interface will pop up and let you work.
   
           
For any question or issue send an email at:
    antonello.trullo@gmail.com            
  
