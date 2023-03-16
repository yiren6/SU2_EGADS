"""
Sample drive script for running a cross-flow heat exchanger test case.
The geometry of the fluid domain is ./channel.csm, which is a EGADS geometry
script that uses a constructive geometry approach to build 3D objectives
from premitives.

You need to install SU2, pysu2, pysu2ad, EGADS, CAPS with pyCAPS support
to satisfy the minimum environment requirements. Depending on the mesher
used, you may also need Pointwise or Tetgen installed. Please refer to the
EGADS and SU2 installation guide for installation help.

(c) Yiren Shen 2023
"""
# Import modules
# Import pyCAPS module

import pyCAPS
import os
import argparse
import platform
import time
import pysu2
import pysu2ad
import csv
import shutil
from utils import pysu2_egadsUtils
from utils import mpi
from parallel_computation import parallel_computation as su2Run

## Main script
# Import SU2 Python interface module


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description="Cross flow heat exchanger sensitivity Example",
                                 prog="channel_test",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Setup the available commandline options
parser.add_argument("-workDir", default=["." + os.sep], nargs=1, type=str, help="Set working/run directory")
parser.add_argument("-numberProc", default=1, nargs=1, type=float, help="Number of processors")
parser.add_argument("-outLevel", default=1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

project_name = "channel_test_case"
workDir = os.path.join(str(args.workDir[0]), project_name)  # log the current work dir

# Load CSM file
geometryScript = os.path.join(".", "channel.csm")  # point toward your current csm geom file
# build pyCAPS instance
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
# change design parameter
myProblem.geometry.despmtr.crosCircDiam = 0.25

###########Pointwise surface/volume mesh
"""
Pointwise was used using pyCAPS pointwise API. Alternatively
you can also use the latter section, the TETGEN API to do volume 
mesh with Tetgen. Usually Tetgen is not recommended to use for 
RANS problem, as the volume mesh has little boundary layer controls.
"""
# Load pointwise aim
pointwise = myProblem.analysis.create(aim="pointwiseAIM",
                                      name="pointwise")

# Global Min/Max number of points on edges
# you can manually set these parameters to change connector dim
# pointwise.input.Connector_Initial_Dim = 26
# pointwise.input.Connector_Min_Dim     = 22
# pointwise.input.Connector_Max_Dim     = 30
# pointwise.input.Connector_Turn_Angle = 10.0
# pointwise.input.Connector_Prox_Growth_Rate = 0.9

pointwise.input.Domain_Algorithm = "AdvancingFrontOrtho"
# pointwise.input.Domain_Adapt = True
pointwise.input.Domain_Max_Edge = 1.5  # relative to capsMeshLength in csm file for each domain

pointwise.input.Block_Full_Layers = 1
pointwise.input.Block_Max_Layers  = 30
pointwise.input.Block_Growth_Rate = 1.2
# Set mesh sizing parmeters, only Wing2 is viscous
viscousBC = {"boundaryLayerSpacing": 0.001,
             "boundaryLayerGrowthRate": 1.2,
             "tess_Params": [0.001, 0.01, 15]}
# tess_param physical meaning is [cell_size_relative_to_meshLength in csm file,
#           maximum_deviation_from_surface_relative_to_meshLength in csm file, maximum_turning_angle in degree]
pointwise.input.Mesh_Sizing = {"crossbox": viscousBC}

# Run AIM pre-analysis
pointwise.preAnalysis()

####### Run pointwise ####################
CAPS_GLYPH = os.environ["CAPS_GLYPH"]
for i in range(1):
    try:
        if platform.system() == "Windows":
            PW_HOME = os.environ["PW_HOME"]
            pointwise.system(
                PW_HOME + "\\win64\\bin\\tclsh.exe " + CAPS_GLYPH + "\\GeomToMesh.glf caps.egads capsUserDefaults.glf")
        else:
            pointwise.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
    except pyCAPS.CAPSError:
        time.sleep(10)  # wait and try again
        continue

    time.sleep(1)  # let the harddrive breathe
    if os.path.isfile(os.path.join(pointwise.analysisDir, "caps.GeomToMesh.gma")) and \
            (os.path.isfile(os.path.join(pointwise.analysisDir, "caps.GeomToMesh.ugrid")) or \
             os.path.isfile(os.path.join(pointwise.analysisDir, "caps.GeomToMesh.lb8.ugrid"))): break
    time.sleep(10)  # wait and try again
##########################################

# Run AIM post-analysis
pointwise.postAnalysis()

#### TETGEN :: uncomment this block if you are using tetgen as mesher
# mySurfMesh = myProblem.analysis.create(aim = "egadsTessAIM",
#                                        name = "tess")
# myMesh = myProblem.analysis.create(aim = "tetgenAIM",
#                                    name = "myMesh")
# # Set project name so a mesh file is generated
# mySurfMesh.input.Proj_Name = "egadsTessMesh"
#
# # Set new EGADS body tessellation parameters
# mySurfMesh.input.Tess_Params = [0.5, 0.1, 20.0]
#
# # Set output grid format since a project name is being supplied - Tecplot file
# mySurfMesh.input.Mesh_Format = "Tecplot"
#
# # Link surface mesh from EGADS to TetGen
# myMesh.input["Surface_Mesh"].link(mySurfMesh.output["Surface_Mesh"])
#
# # Preserve surface mesh while meshing
# myMesh.input.Preserve_Surf_Mesh = True


############ Load SU2 aim

su2 = myProblem.analysis.create(aim="su2AIM", name="su2")

su2.input["Mesh"].link(myProblem.analysis["pointwise"].output["Volume_Mesh"])

# for tetgen :: uncomment the following line for tetgen generated volume mesh
# su2.input["Mesh"].link(myMesh.output["Volume_Mesh"])

# Set SU2 Version
su2.input.SU2_Version = "Blackbird"
# Set project name
su2.input.Proj_Name = project_name
# Set AoA number
su2.input.Alpha = 0.0
# Set Mach number
su2.input.Mach = 0.01

# Set equation type
su2.input.Equation_Type = "Compressible"
su2.input.Physical_Problem = "INC_NAVIER_STOKES"
su2.input.Turbulence_Model = "NONE"
su2.input.Re = 1000
# Set the output file formats
su2.input.Output_Format = "Paraview"

# Set number of iterations
su2.input.Num_Iter = 800

# Set boundary conditions
BC1 = {"bcType": "Viscous", "wallTemperature": 600}
symBC = {"bcType": "Symmetry"}
backPressureBC = {"bcType": "SubsonicOutflow",
                  "staticPressure": 0.0}
inflow = {"bcType": "SubsonicInflow",
          "totalPressure": 0.2,
          "totalTemperature": 1900,
          "uVelocity": 1.0,
          "vVelocity": 0.0,
          "wVelocity": 0.0}
# set boundary conditions, the domain name corresponds to surface name
#   set in channel.csm
su2.input.Boundary_Condition = {"crossbox": BC1,
                                "tops": symBC,
                                "sides": symBC,
                                "inlet": inflow,
                                "outlet": backPressureBC}
su2.input.Reference_Area = 1
su2.input.Freestream_Velocity = 0.2
su2.input.Design_Functional = {"Drag": {
    "function": "TOTAL_HEATFLUX"}}  # "Drag" is a place holder here, not relevent to the objective function, which is total heat_flux
su2.input.Freestream_Pressure = 101300.0
su2.input.Freestream_Temperature = 1900
su2.input.Reference_Dimensionalization = "DIMENSIONAL"
su2.input.MultiGrid_Level = 0
# Declare design variables
su2.input.Design_Variable = {"crosCircDiam": "0.25"}
# Specify the boundaries used to compute forces
su2.input.Surface_Monitor = ["crossbox"]
su2.input.Convective_Flux = "JST"
# additional cfg file inputs, as an array
su2.input.Input_String = ["MARKER_DESIGNING= BC_1",
                          "INC_DENSITY_MODEL= CONSTANT",
                          "INC_ENERGY_EQUATION = YES",
                          "INC_DENSITY_INIT= 1",
                          "INC_VELOCITY_INIT= (0.2, 0, 0 )",
                          "INC_TEMPERATURE_INIT= 1900",
                          "INC_NONDIM= DIMENSIONAL",
                          "FLUID_MODEL= INC_IDEAL_GAS",
                          "VISCOSITY_MODEL= CONSTANT_VISCOSITY ",
                          "MU_CONSTANT = 2E-5",
                          "SPECIFIC_HEAT_CP= 1004.7036",
                          "INC_INLET_TYPE= VELOCITY_INLET",
                          "INC_OUTLET_TYPE= PRESSURE_OUTLET",
                          "INC_INLET_DAMPING = 0.05",
                          "INC_OUTLET_DAMPING = 0.05",
                          "MARKER_ANALYZE= ( outlet )",
                          "MARKER_ANALYZE_AVERAGE = AREA",
                          "NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES",
                          "CFL_ADAPT= YES ",
                          "CFL_ADAPT_PARAM= ( 0.1, 2.0, 1.0, 1E10 )",
                          "OUTPUT_WRT_FREQ= 100",
                          "LINEAR_SOLVER= FGMRES",
                          "LINEAR_SOLVER_PREC= ILU",
                          "LINEAR_SOLVER_ERROR= 1E-10",
                          "LINEAR_SOLVER_ITER=20",
                          "MUSCL_FLOW= YES",
                          "CONV_RESIDUAL_MINVAL= -10",
                          "CONV_STARTITER= 10",
                          "CONV_CAUCHY_ELEMS= 100",
                          "CONV_CAUCHY_EPS= 1E-6",
                          "VOLUME_OUTPUT= SOLUTION PRIMITIVE SOURCE RESIDUAL SENSITIVITY LOOKUP TIMESTEP",
                          "OBJECTIVE_FUNCTION= TOTAL_HEATFLUX"]

# Run AIM pre-analysis
su2.preAnalysis()

####### Run SU2 ######################
print("\n\nRunning SU2......")

current_dir = os.getcwd()  # Get our current working directory
os.chdir(su2.analysisDir)  # Move into test directory

# run su2 using pysu2 API
su2solver = pysu2.CSinglezoneDriver(project_name + ".cfg", 1, mpi.COMM)
su2solver.Preprocess(0)
su2solver.Run()
su2solver.Postprocess()
su2solver.Monitor(0)
su2solver.Output(0)

# obtain objective function val
objective_val = su2solver.GetObjective()
cell_normal = su2solver.GetMarkerVertexNormals(0)
cell_area = pysu2_egadsUtils.euclidean_norm(cell_normal)

print("Object function total heatflux is ", objective_val)

#######################################

# Run AIM post-analysis
su2.postAnalysis()
current_dir = os.getcwd()
# copy the restart file to solution file
old_file = f"{current_dir}/restart_flow_{project_name}.dat"
new_file = f"{current_dir}/solution_adj_cd.dat"
shutil.copy(old_file, new_file)
# copy the history file of direct sol
current_dir = os.getcwd()
old_file = current_dir + "/history.csv"
new_file = current_dir + "/history_direct.csv"
shutil.copy(old_file, new_file)

############ Adjoint calculation

# Request geometric design sensitivities for the Design_Functional
su2.input.Design_Sensitivity = True
# Declare design variables
su2.input.Design_Functional = {"Drag": {"function": "TOTAL_HEATFLUX"}}
su2.input.Design_Variable = {"crosCircDiam": "0.25"}
su2.input.Convective_Flux = "JST"
su2.input.Surface_Monitor = ["crossbox"]
su2.input.Num_Iter = 1400
su2.input.Input_String = ["MARKER_DESIGNING= BC_1",
                          "INC_DENSITY_MODEL= CONSTANT",
                          "INC_ENERGY_EQUATION = YES",
                          "INC_DENSITY_INIT= 1",
                          "INC_VELOCITY_INIT= (0.2, 0, 0 )",
                          "INC_TEMPERATURE_INIT= 1900",
                          "INC_NONDIM= DIMENSIONAL",
                          "FLUID_MODEL= INC_IDEAL_GAS",
                          "VISCOSITY_MODEL= CONSTANT_VISCOSITY ",
                          "MU_CONSTANT = 2E-5",
                          "SPECIFIC_HEAT_CP= 1004.7036",
                          "INC_INLET_TYPE= VELOCITY_INLET",
                          "INC_OUTLET_TYPE= PRESSURE_OUTLET",
                          "INC_INLET_DAMPING = 0.05",
                          "INC_OUTLET_DAMPING = 0.05",
                          "MARKER_ANALYZE= ( outlet )",
                          "MARKER_ANALYZE_AVERAGE = AREA",
                          "NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES",
                          "CFL_ADAPT= YES ",
                          "CFL_ADAPT_PARAM= ( 0.1, 2.0, 1.0, 1E10 )",
                          "OUTPUT_WRT_FREQ= 100",
                          "LINEAR_SOLVER= FGMRES",
                          "LINEAR_SOLVER_PREC= ILU",
                          "LINEAR_SOLVER_ERROR= 1E-10",
                          "LINEAR_SOLVER_ITER=20",
                          "MUSCL_FLOW= YES",
                          "CONV_RESIDUAL_MINVAL= -10",
                          "CONV_STARTITER= 10",
                          "CONV_CAUCHY_ELEMS= 100",
                          "CONV_CAUCHY_EPS= 1E-6",
                          "VOLUME_OUTPUT= SOLUTION PRIMITIVE SOURCE RESIDUAL SENSITIVITY LOOKUP TIMESTEP",
                          "SINGLEZONE_DRIVER= YES",
                          "DISCADJ_LIN_SOLVER= FGMRES",
                          "DISCADJ_LIN_PREC= ILU",
                          "SURFACE_ADJ_FILENAME= surface_adjoint"]

# Run AIM pre-analysis
su2.preAnalysis()

############ Adjoint calculation
print("\n\nRunning SU2 Adjoint......")
current_dir = os.getcwd()  # Get our current working directory
os.chdir(su2.analysisDir)  # Move into test directory

# Run Adjoint
su2ad = pysu2ad.CDiscAdjSinglezoneDriver(f"{project_name}_adjoint.cfg", 1, mpi.COMM)

su2ad.Preprocess(0)
su2ad.Run()
su2ad.Postprocess()
su2ad.Monitor(0)
su2ad.Output(0)

# extract residual from partial I over partial xv
surface_selected = "BC_1"
result_surface_coordinate_sensitivity = pysu2_egadsUtils.read_surface_sensitivity("SURFACE_ADJOINT.csv",
                                                                                  surface_selected, is_euler=False)
result_surface_coordinate_sensitivity = pysu2_egadsUtils.multiple_elementwise(result_surface_coordinate_sensitivity,
                                                                              cell_area, surface_selected)

marker_index = pysu2_egadsUtils.read_marker_index("SURFACE_ADJOINT.csv")
marker_index = {"BC_1": marker_index}
# ensure objective val is float
objective_val = float(objective_val)

# with updated surface sensitivity routine
# write surface sensitivity file into txt file
pysu2_egadsUtils.write_marker_sensitivity(result_surface_coordinate_sensitivity, marker_index, [surface_selected],
                                          "HeatFlux", objective_val, f"{project_name}.sens")
su2solver.Postprocessing()
su2ad.Postprocessing()

# Run AIM post-analysis
su2.postAnalysis()

############# calculate parametric sensitivity
#   this was done using the dynout function in SU2AIM in pyCAPS
dHeatfluxdHeight = su2.dynout["HeatFlux"].deriv("crosCircDiam")

print("dHeatFlux/dcrosCircDiam", dHeatfluxdHeight)
