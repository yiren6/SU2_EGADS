# Import pyCAPS module
import mpi4py.MPI
import numpy as np
import pyCAPS
import pyEGADS.egads
from pyEGADS import egads
# Import modules
import os
import argparse
import platform
import time
import pysu2
import pysu2ad
import csv
import shutil
import scipy

# subfunctions for marker vertices extraction and sensitivity writing

def unique(repList):
 # returns sorted unique elements list from a input list
 # input: unsorted list, list
 # output: sorted list, list
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in repList:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    unique_list.sort()
    return unique_list

def findBoundaryIndices(su2filepath):
    # Function Read and extract boundary condition and index from SU2 file data
    # input: relative path to su2 mesh file, string
    # output: marker vertices index for all surface marker, dict

    su2file_contents = {}
    # read file definition
    with open(su2filepath) as su2file:
        su2file=su2file.readlines()
        print("%d rows in SU2 grid file \n" % len(su2file))
        lines = int(0)
        ii = 1
        for lines in range(len(su2file)):
            current_line= su2file[lines]
            #print(current_line)
            # extract dimension data
            if "NDIME" in current_line:
                su2file_contents['NDIME'] = [int(i) for i in current_line.split() if i.isdigit()]
                su2file_contents['NDIME_rowidx'] = lines
            elif "NELEM" in current_line:
                su2file_contents['NELEM'] = [int(i) for i in current_line.split() if i.isdigit()]
                su2file_contents['NELEM_rowidx'] = lines
            elif "NPOIN" in current_line:
                su2file_contents['NPOIN'] = [int(i) for i in current_line.split() if i.isdigit()]
                su2file_contents['NPOIN_rowidx'] = lines
            elif "NMARK" in current_line:
                su2file_contents['NMARK'] = [int(i) for i in current_line.split() if i.isdigit()]
                su2file_contents['NMARK_rowidx'] =lines
            # loop and extract all marker info
            elif "MARKER_TAG" in current_line:
                su2file_contents['MARKER_'+str(ii)] = current_line.split('=')[1].strip() #remove \n
                su2file_contents['MARKER_'+str(ii)+'_ELEMS'] = int(su2file[lines+1].split('=')[1].strip())
                su2file_contents['MARKER_'+str(ii)+'_rowidx'] = lines
                ii += 1
    # print .su2 file summary
    print("File summary: \n")
    print(su2file_contents)

    # extract vertices index for all marker

    numMarker = su2file_contents['NMARK'][0]
    markerVerticesIdx = {}
    ii = 0
    for ii in range(numMarker):
        row_start = su2file_contents['MARKER_'+str(ii+1)+'_rowidx']+2 # data starts at 2 rows below
        cur_marker_numElem = su2file_contents['MARKER_'+str(ii+1)+'_ELEMS']
        # extract all row for current marker
        cur_marker_vertices_idx = []
        for curRow in range(row_start,row_start+cur_marker_numElem):
            curline = list(map(int,(su2file[curRow].split(' ')[1:-1]) ))
            #print(curline)
            cur_marker_vertices_idx = cur_marker_vertices_idx + curline
        # append current marker list to data struct
        markerVerticesIdx[su2file_contents['MARKER_'+str(ii+1)]] = unique(cur_marker_vertices_idx)

    return markerVerticesIdx

def writeMarkerSensitivity(sensitivityData, markerIndex, selectedMarker,ObjectiveVal, saveFileName):
    # Function write the surface sensitivity text file to current workspace
    # Input:
    # sensitivityData: the surface sensitivity data obtained from pysu2ad.getObjectiveCoordinateSensitivity, dict
    # marker_index: the vertices index for all surface marker, dict
    # selectedMarker: marker names selected for export, list
    # saveFileName: name of the export text file, string
    # Output: surface sensitivity text file

    with open(saveFileName,'a') as outfile:
        cur_marker_idx = 0 # current Marker Index
        cur_line_idx = 0   # current text file line Index
        accum_line  = 0   # dummy variable for total line wrote

        if cur_line_idx == 0:
            curline = "1 0 \n" # nfunctions nDesignVariables
            outfile.writelines(curline)
            cur_line_idx += 1
        if cur_line_idx == 1:
            curline = "%s \n" # function key
            outfile.writelines(curline%("Drag"))
            cur_line_idx += 1
            print(cur_line_idx)
        if cur_line_idx == 2:
            curline = "%s \n" #objective value
            outfile.writelines(curline%(ObjectiveVal))
            cur_line_idx += 1
            print(cur_line_idx)
        if cur_line_idx == 3:
            numTotalVertices = sum((len(selectedMarker) for selectedMarker in sensitivityData.values()))
            print(numTotalVertices)
            curline = "%s \n" # nnodes
            outfile.writelines(curline%(numTotalVertices ) )
            cur_line_idx += 1


        # write line by line for all surface markers
        for ii in range(len(selectedMarker)):
            curData = sensitivityData[selectedMarker[cur_marker_idx]]

            for jj in range(len(curData)):

                outfile.writelines("%s %s %s %s\n"%(markerIndex[selectedMarker[cur_marker_idx]][jj]+1,             \
                                                       sensitivityData[surface_selected[cur_marker_idx]][jj][0],    \
                                                       sensitivityData[surface_selected[cur_marker_idx]][jj][1],    \
                                                       sensitivityData[surface_selected[cur_marker_idx]][jj][2]))
                cur_line_idx += 1

            accum_line += len(curData)
            cur_marker_idx += 1

    return 0

def readObjectiveFunctionValue (forceBreakDownFilePath, ObjectiveValueLabel):
    ObjectiveValue = []
    with open(forceBreakDownFilePath) as forceBreakDownFile:
        lines = forceBreakDownFile.readlines()
        for l in lines: #start searching for data beginning from surface summary line
            if len(l) <= 90:
                continue
            else:
                tempLineContVec = l.split()
                if tempLineContVec[1] == ObjectiveValueLabel:
                    ObjectiveValue = tempLineContVec[4]
                    break
    return ObjectiveValue

def readSurfaceSensitivity (surfaceSensitivityFileName, selectedMarker):
    # * surfaceSensitivityFileName is a string of
    #       the file name of adjoint objective surface sensitivity
    #       usually named as "SURFACE_ADJOINT.csv"
    # * selectedMarker is a string of the name of the selected marker
    # * function output a dict of selected marker with a list of all
    #       corresponding marker surface objective sensitivities in x,y,z
    surfaceSens = ()
    with open(surfaceSensitivityFileName) as surfaceAdjoint:
        lines = surfaceAdjoint.readlines()
        for l in lines[1:]: #start extracting data from 2nd to end
            if len(l) <= 90: # remove empty lines
                continue
            else:
                tempLineContVec = l.split(',')
                #ptID = tempLineContVec[0]
                sensitivity_x = tempLineContVec[9]
                sensitivity_y = tempLineContVec[10]
                sensitivity_z = tempLineContVec[11]
            currentTuple = (sensitivity_x, sensitivity_y, sensitivity_z)
            surfaceSens = (*surfaceSens, currentTuple)
    surfaceSensdict = {selectedMarker: list(surfaceSens)}
    return surfaceSensdict

## Main script
# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'box sensitivity unit test Example',
                                 prog = 'box_unit_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

projectName = "boxSens"
workDir = os.path.join(str(args.workDir[0]), projectName)


# Load CSM file
geometryScript = os.path.join(".","block.csm")

myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
# change design parameter
myProblem.geometry.despmtr.blkheight = 1.25

# Load pointwise aim
pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                      name = "pointwise")

# Global Min/Max number of points on edges
pointwise.input.Connector_Initial_Dim = 31
pointwise.input.Connector_Min_Dim     = 25
pointwise.input.Connector_Max_Dim     = 35
#pointwise.input.Connector_Turn_Angle = 10.0
#pointwise.input.Connector_Prox_Growth_Rate = 0.9

pointwise.input.Domain_Algorithm = "AdvancingFrontOrtho"
#pointwise.input.Domain_Adapt = True
pointwise.input.Domain_Max_Edge = 1.0 #relative to capsMeshLength in csm file for each domain

pointwise.input.Block_Full_Layers = 1
pointwise.input.Block_Max_Layers  = 7
pointwise.input.Block_Growth_Rate = 2
# Set mesh sizing parmeters, only Wing2 is viscous
viscousBC  = {"boundaryLayerSpacing" : 0.2,
              "boundaryLayerGrowthRate" : 2,
              "tess_Params": [0.12,0.01,15]}
pointwise.input.Mesh_Sizing = {"body": viscousBC}

# Run AIM pre-analysis
pointwise.preAnalysis()

####### Run pointwise ####################
CAPS_GLYPH = os.environ["CAPS_GLYPH"]
for i in range(1):
    try:
        if platform.system() == "Windows":
            PW_HOME = os.environ["PW_HOME"]
            pointwise.system(PW_HOME + "\\win64\\bin\\tclsh.exe " + CAPS_GLYPH + "\\GeomToMesh.glf caps.egads capsUserDefaults.glf")
        else:
            pointwise.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
    except pyCAPS.CAPSError:
        time.sleep(10) # wait and try again
        continue

    time.sleep(1) # let the harddrive breathe
    if os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.gma')) and \
      (os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.ugrid')) or \
       os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.lb8.ugrid'))): break
    time.sleep(10) # wait and try again
##########################################

# Run AIM post-analysis
pointwise.postAnalysis()


############ Load SU2 aim

su2 = myProblem.analysis.create(aim = "su2AIM",
                                name = "su2")

su2.input["Mesh"].link(myProblem.analysis["pointwise"].output["Volume_Mesh"])
# Set SU2 Version
su2.input.SU2_Version = "Blackbird"

# Set project name
su2.input.Proj_Name = projectName

# Set AoA number
su2.input.Alpha = 0.0

# Set Mach number
su2.input.Mach = 0.3


# Set equation type
su2.input.Equation_Type = "Compressible"
su2.input.Physical_Problem = "EULER"

# Set the output file formats
su2.input.Output_Format = "Paraview"

# Set number of iterations
su2.input.Num_Iter = 800

# Set boundary conditions
inviscidBC1 = {"bcType" : "Inviscid"}
symbc = {"bcType" : "Symmetry"}
su2.input.Boundary_Condition = {"body"   : inviscidBC1,
                                "sym"  : symbc,
                                "channel"  : symbc,
                                "inlet":  "farfield",
                                "outlet": "farfield"}
su2.input.Reference_Area = 1
su2.input.Freestream_Velocity = 102.1
su2.input.Design_Functional = {"Drag": {'function':'DRAG'}}
su2.input.Freestream_Pressure = 101300.0
su2.input.Reference_Dimensionalization = "FREESTREAM_PRESS_EQ_ONE"
# Declare design variablesAZ
su2.input.Design_Variable = {"blkheight":"1.25"}
su2.input.MultiGrid_Level = 0

# Specifcy the boundares used to compute forces
su2.input.Surface_Monitor = ["body"]
su2.input.Convective_Flux = "JST"
su2.input.Input_String = ["CFL_ADAPT= YES ","MARKER_DESIGNING= BC_1"]

#"MARKER_DEFORM_MESH= BC_1", "DV_KIND= NO_DEFORMATION", "DV_MARKER= BC_1"

# Run AIM pre-analysis
su2.preAnalysis()

####### Run SU2 ######################
print ("\n\nRunning SU2......")

currentDirectory = os.getcwd() # Get our current working directory
os.chdir(su2.analysisDir) # Move into test directory

comm=mpi4py.MPI.COMM_WORLD
su2solver = pysu2.CSinglezoneDriver(projectName + ".cfg", 1, comm)
su2solver.Preprocess(0)

su2solver.Run()

su2solver.Postprocess()
su2solver.Monitor(0)
su2solver.Output(0)


I  = su2solver.GetObjective()

print('Object function CD is ', I)
#su2solver.Postprocessing()

#######################################

# Run AIM post-analysis
su2.postAnalysis()
currentDirectory = os.getcwd()
oldfile = currentDirectory+"/restart_flow_boxSens.dat"
newfile = currentDirectory+"/solution_adj_cd.dat"
shutil.copy(oldfile, newfile)

############ Adjoint calculation

# Request geometric design sensitivities for the Design_Functional
su2.input.Design_Sensitivity = True
# Declare design variables
su2.input.Design_Functional = {"Drag": {'function':'DRAG'}}
su2.input.Design_Variable = {"blkheight":"1.25"}
su2.input.Convective_Flux = "JST"
su2.input.Surface_Monitor = ["body"]

su2.input.Input_String = ["CFL_ADAPT= NO ", "SINGLEZONE_DRIVER= YES",
                          "DISCADJ_LIN_SOLVER= FGMRES", "DISCADJ_LIN_PREC= ILU",
                            "LINEAR_SOLVER_ITER= 20",
                            "SURFACE_ADJ_FILENAME= surface_adjoint", "MARKER_DESIGNING= BC_1",
                          "LINEAR_SOLVER_ERROR= 1E-5"]


# Run AIM pre-analysis
su2.preAnalysis()

############ Adjoint calculation
print ("\n\nRunning SU2 Adjoint......")
currentDirectory = os.getcwd() # Get our current working directory
os.chdir(su2.analysisDir) # Move into test directory

# Run Adjoint
comm=mpi4py.MPI.COMM_WORLD  # define intracommunicator instance
su2ad = pysu2ad.CDiscAdjSinglezoneDriver(projectName + "_adjoint.cfg", 1, comm)
#su2ad = pysu2ad.CDiscAdjSinglezoneDriver(adjoint_cfg_file_name, 1, comm)


#injectResidual = (10**10)* np.ones((3189,5))

su2ad.Preprocess(0)
#su2ad.SetAdjointStates(injectResidual)
su2ad.Run()

su2ad.Postprocess()
su2ad.Monitor(0)
su2ad.Output(0)

# extract sensitivity for volume mesh
result_coordinate_sensitivity = su2ad.GetObjectiveCoordinatesSensitivities()
result_farfield_sensitivity = su2ad.GetObjectiveFarfieldVariablesSensitivities()


#result_coordinate_sensitivity = su2ad.GetMarkerObjectiveDisplacementsSensitivities(0)


# extract residual from partial I over partial xv
result_objective_state_sensitivities = su2ad.GetObjectiveStatesSensitivities()
result_residual_state_sensitivities = su2ad.GetResidualsStatesSensitivities()

print(result_objective_state_sensitivities)

# extract sensitivity for surface mesh
markerVerticesIndex = findBoundaryIndices("../pointwise/caps.GeomToMesh.su2")

# write volume sensitivity file
with open('Adj_Coordinate_Sensitivities.csv','w',newline='') as csvfile:
    my_writer = csv.writer(csvfile, delimiter= ' ')
    my_writer.writerow(result_coordinate_sensitivity)

# write surface sensitivity file
result_surface_coordinate_sensitivity = {}
surface_selected = ['BC_1']

# filter out surface sensitivity vertices
for ii in range(len(surface_selected)):
    result_surface_coordinate_sensitivity[surface_selected[ii]] = \
        [result_coordinate_sensitivity[jj] for jj in markerVerticesIndex[surface_selected[ii] ] ]

## obtain objective value from froce break down file

ObjectiveVal = readObjectiveFunctionValue ("../su2/forces_breakdown_boxSens.dat", "CD")
ObjectiveVal = float(ObjectiveVal)

# with updated surface sensitivity routine
result_surface_coordinate_sensitivity = readSurfaceSensitivity('SURFACE_ADJOINT.csv', surface_selected[0])
# write surface sensitivity file into txt file
writeMarkerSensitivity(result_surface_coordinate_sensitivity, markerVerticesIndex,surface_selected, ObjectiveVal, \
                       projectName + ".sens")

su2solver.Postprocessing()
su2ad.Postprocessing()

# Run AIM post-analysis
su2.postAnalysis()


# calculate adjoint
dDragdHeight = su2.dynout["Drag"].deriv("blkheight")

print('dDrag/dblkheight', dDragdHeight)




