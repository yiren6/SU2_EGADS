!
!      EGADS: Electronic Geometry Aircraft Design System
!
!             FORTRAN 90 Include
!
!      Copyright 2011-2022, Massachusetts Institute of Technology
!      Licensed under The GNU Lesser General Public License, version 2.1
!      See http://www.opensource.org/licenses/lgpl-2.1.php
!
!
!       OBJECT CLASSES
!
        INTEGER, PARAMETER :: CONTXT=0, TRANSFORM=1, TESSELLATION=2
        INTEGER, PARAMETER :: REFERENCE=5
        INTEGER, PARAMETER :: PCURVE=10, CURVE=11, SURFACE=12
        INTEGER, PARAMETER :: NODE=20, EDGE=21, LOOP=22, FACE=23
        INTEGER, PARAMETER :: SHELL=24, BODY=25, MODEL=26

!       MEMBER TYPES
!
!         PCURVES & CURVES
        INTEGER, PARAMETER :: LINE=1, CIRCLE=2, ELLIPSE=3, PARABOLA=4
        INTEGER, PARAMETER :: HYPERBOLA=5, TRIMMED=6, BEZIER=7, BSPLINE=8
        INTEGER, PARAMETER :: OFFSET=9
!
!         SURFACES
        INTEGER, PARAMETER :: PLANE=1, SPHERICAL=2, CYLINDRICAL=3, REVOLUTION=4
        INTEGER, PARAMETER :: TOROIDAL=5, CONICAL=10, EXTRUSION=11
!
!         TOPOLOGY
        INTEGER, PARAMETER :: SREVERSE=-1, NOMTYPE=0, SFORWARD=1, ONENODE=1
        INTEGER, PARAMETER :: TWONODE=2, OPEN=3, CLOSED=4, DEGENERATE=5
        INTEGER, PARAMETER :: WIREBODY=6, FACEBODY=7, SHEETBODY=8, SOLIDBODY=9
!
!	ATTRIBUTE TYPES
!
        INTEGER, PARAMETER :: ATTRINT=1, ATTRREAL=2, ATTRSTRING=3, ATTRCSYS=12
!
!       SOLID BOOLEAN OPERATIONS
!
        INTEGER, PARAMETER :: SUBTRACTION=1, INTERSECTION=2, FUSION=3
!
!	SOLID BODY TYPES
!
        INTEGER, PARAMETER :: BOX=1, SPHERE=2, CONE=3, CYLINDER=4, TORUS=5
!
!	NULL
!
        INTEGER(KIND=8), PARAMETER :: NULL=0