#
export ESP_ARCH=LINUX64
export ESP_ROOT=/home/aarons/Desktop/ESP1.21/EngSketchPad
export CASROOT=/home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1
export CASARCH=Linux
export CASREV=7.4
export PATH=/home/aarons/Desktop/ESP1.21/EngSketchPad/bin:$PATH
export PYTHONINC="/usr/include/python3.8"
export PYTHONLIB="/usr/lib/x86_64-linux-gnu/libpython3.8.so.1.0"
export PYTHONPATH=/home/aarons/Desktop/ESP1.21/EngSketchPad/pyESP:$PYTHONPATH
export EFCOMP=gfortran
#export AFLR=/home/aarons/Projects/AFLR
export AFLR_ARCH=Linux-x86-64
if [ -z "$LD_LIBRARY_PATH" ]; then
    export LD_LIBRARY_PATH=/home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1/Linux/lib:/home/aarons/Desktop/ESP1.21/EngSketchPad/lib
else
    export LD_LIBRARY_PATH=/home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1/Linux/lib:/home/aarons/Desktop/ESP1.21/EngSketchPad/lib:$LD_LIBRARY_PATH
fi
export UDUNITS2_XML_PATH=/home/aarons/Desktop/ESP1.21/EngSketchPad/src/CAPS/udunits/udunits2.xml
export CAPS_GLYPH=/home/aarons/Desktop/ESP1.21/EngSketchPad/src/CAPS/aim/pointwise/glyph
#export TETGEN=/home/aarons/Projects/TetGen/tetgen1.6.0
#export AWAVE=/home/aarons/Projects/awave/awavemod.f
export SLUGS_START="firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/SLUGS/Slugs.html &"
export ESP_START="firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/ESP/ESP.html &"
export WV_START="firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/wvClient/wv.html &"
