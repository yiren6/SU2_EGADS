 ESP: The Engineering Sketch Pad
 Rev 1.21 -- July 2022

1. Distribution layout
	README.md	-	this file
	EngSketchPad	-	the $ESP_ROOT directory
	Examples	-	a directory with SU2 + CAPS example
	OCC741lin64	-	the 64 bit OpenCascade release for linux
	tetgen		-	volume mesh generator tetgen

2. ESP installation
	Please refer to the ESP installation guide at $ESP_ROOT/Readme.txt.
	$ESP_ROOT : https://acdl.mit.edu/ESP/
	
3. Run examples:
	The Python examples were included in the Examples folder. Please direct to the example subfolder and run python3 to execute the driver script. You need to ensure all the dependencies were included in $PATH to import python library successfully.

4. tetgen isntallation
	Please refer to the tetgen installation guide at 
	https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1
		
5. Environment Variables
	The following environment variables needs to be defined in ~/.bashrc or other bash profile files and sourced prior to running any scripts.
	 
	ESP/CAPS:
	ESP_ARCH
	ESP_ROOT
	CASROOT
	CASARCH
	CASREV
	^ these environment variables will be generated during ESP installation 
	PYTHONINC
	PYTHONLIB
	^ these enviroment variables, you need to set during EGADS installation
	PYTHONPATH
	^ make sure you add SU2 and EGADS PATH to system $PATH and $PYTHONPATH
	CAPS_PATH
	PW_HOME
	^ you need these set to use PointwiseAIM
	TETGEN
	^ you need this pointing to tetgen folder to use tetgenAIM
	
	

	 

 
