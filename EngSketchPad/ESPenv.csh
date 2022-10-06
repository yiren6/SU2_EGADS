#
setenv ESP_ARCH LINUX64
setenv ESP_ROOT /home/aarons/Desktop/ESP1.21/EngSketchPad
setenv CASROOT /home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1
setenv CASARCH Linux
setenv CASREV 7.4
setenv PATH /home/aarons/Desktop/ESP1.21/EngSketchPad/bin:$PATH
setenv PYTHONINC "/usr/include/python3.8"
setenv PYTHONLIB "/usr/lib/x86_64-linux-gnu/libpython3.8.so.1.0"
setenv PYTHONPATH /home/aarons/Desktop/ESP1.21/EngSketchPad/pyESP:$PYTHONPATH
setenv EFCOMP gfortran
#setenv AFLR /home/aarons/Projects/AFLR
setenv AFLR_ARCH Linux-x86-64
if ( $?LD_LIBRARY_PATH == 0 ) then
    setenv LD_LIBRARY_PATH /home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1/Linux/lib:/home/aarons/Desktop/ESP1.21/EngSketchPad/lib
else
    if ( "$LD_LIBRARY_PATH" == "" ) then
        setenv LD_LIBRARY_PATH /home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1/Linux/lib:/home/aarons/Desktop/ESP1.21/EngSketchPad/lib
    else
        setenv LD_LIBRARY_PATH /home/aarons/Desktop/ESP1.21/OCC741lin64/OpenCASCADE-7.4.1/Linux/lib:/home/aarons/Desktop/ESP1.21/EngSketchPad/lib:$LD_LIBRARY_PATH
    endif
endif
setenv UDUNITS2_XML_PATH /home/aarons/Desktop/ESP1.21/EngSketchPad/src/CAPS/udunits/udunits2.xml
setenv CAPS_GLYPH /home/aarons/Desktop/ESP1.21/EngSketchPad/src/CAPS/aim/pointwise/glyph
#setenv TETGEN /home/aarons/Projects/TetGen/tetgen1.6.0
#setenv AWAVE /home/aarons/Projects/awave/awavemod.f
setenv SLUGS_START "firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/SLUGS/Slugs.html &"
setenv ESP_START "firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/ESP/ESP.html &"
setenv WV_START "firefox /home/aarons/Desktop/ESP1.21/EngSketchPad/wvClient/wv.html &"
