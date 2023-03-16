# Utilities function to facilitate EGADS-SU2 I/O
# Use these function by importing pysu2_egadsUtil.py
# Yiren Shen (c) 2023

import numpy as np
import os

# subfunctions for marker vertices extraction and sensitivity writing

def unique(rep_list):
 # returns sorted unique elements list from a input list
 # input: unsorted list, list
 # output: sorted list, list
    # initialize a null list
    unique_list = []
    # traverse for all elements
    for x in rep_list:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    unique_list.sort()
    return unique_list

def findBoundaryIndices(su2_file_path):
    # Function Read and extract boundary condition and index from SU2 file data
    # input: relative path to su2 mesh file, string
    # output: marker vertices index for all surface marker, dict

    su2file_contents = {}
    # read file definition
    with open(su2_file_path) as su2_file:
        su2_file=su2_file.readlines()
        print("%d rows in SU2 grid file \n" % len(su2_file))
        lines = int(0)
        ii = 1
        for lines in range(len(su2_file)):
            current_line= su2_file[lines]

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
                su2file_contents['MARKER_'+str(ii)+'_ELEMS'] = int(su2_file[lines+1].split('=')[1].strip())
                su2file_contents['MARKER_'+str(ii)+'_rowidx'] = lines
                ii += 1
    # print .su2 file summary
    print("File summary: \n")
    print(su2file_contents)

    # extract vertices index for all marker

    num_marker = su2file_contents['NMARK'][0]
    marker_vertices_idx = {}
    ii = 0
    for ii in range(num_marker):
        row_start = su2file_contents['MARKER_'+str(ii+1)+'_rowidx'] + 2 # data starts at 2 rows below
        cur_marker_numElem = su2file_contents['MARKER_'+str(ii+1)+'_ELEMS']
        # extract all row for current marker
        cur_marker_vertices_idx = []
        for cur_row in range(row_start,row_start+cur_marker_numElem):
            cur_line = list(map(int,(su2_file[cur_row].split(' ')[1:-1]) ))
            cur_marker_vertices_idx = cur_marker_vertices_idx + cur_line
        # append current marker list to data struct
        marker_vertices_idx[su2file_contents['MARKER_'+str(ii+1)]] = unique(cur_marker_vertices_idx)

    return marker_vertices_idx

def writeMarkerSensitivity(sensitivity_data, marker_index, selected_marker, objective_name, \
                           objective_val, save_file_name):
    # Function write the surface sensitivity text file to current workspace
    # Input:
    # sensitivity_data: the surface sensitivity data obtained from pysu2ad.getObjectiveCoordinateSensitivity, dict
    # marker_index: the vertices index for all surface marker, dict
    # selected_marker: marker names selected for export, list
    # saveFileName: name of the export text file, string
    # Output: surface sensitivity text file

    with open(save_file_name,'a') as outfile:
        cur_marker_idx = 0 # current Marker Index
        cur_line_idx = 0   # current text file line Index
        accum_line  = 0   # dummy variable for total line wrote

        if cur_line_idx == 0:
            curline = "1 0 \n" # nfunctions nDesignVariables
            outfile.writelines(curline)
            cur_line_idx += 1
        if cur_line_idx == 1:
            curline = "%s \n" # function key
            outfile.writelines(curline%(objective_name))
            cur_line_idx += 1
            print(cur_line_idx)
        if cur_line_idx == 2:
            curline = "%s \n" #objective value
            outfile.writelines(curline%(objective_val))
            cur_line_idx += 1
            print(cur_line_idx)
        if cur_line_idx == 3:
            numTotalVertices = sum((len(selected_marker) for selected_marker in sensitivity_data.values()))
            print(numTotalVertices)
            curline = "%s \n" # nnodes
            outfile.writelines(curline%(numTotalVertices ) )
            cur_line_idx += 1


        # write line by line for all surface markers
        for ii in range(len(selected_marker)):
            curData = sensitivity_data[selected_marker[cur_marker_idx]]

            for jj in range(len(curData)):

                outfile.writelines("%s %s %s %s\n"%(marker_index[selected_marker[cur_marker_idx]][jj],             \
                                                       sensitivity_data[selected_marker[cur_marker_idx]][jj][0],    \
                                                       sensitivity_data[selected_marker[cur_marker_idx]][jj][1],    \
                                                       sensitivity_data[selected_marker[cur_marker_idx]][jj][2]))
                cur_line_idx += 1

            accum_line += len(curData)
            cur_marker_idx += 1

    return 0


def readObjectiveFunctionValue (force_breakDown_file_path, objective_val_label):
    ObjectiveValue = []
    with open(force_breakDown_file_path) as forceBreakDownFile:
        lines = forceBreakDownFile.readlines()
        for l in lines: #start searching for data beginning from surface summary line
            if len(l) <= 90:
                # remove empty line
                continue
            else:
                temp_line_cont_vec = l.split()
                if temp_line_cont_vec[1] == objective_val_label:
                    ObjectiveValue = temp_line_cont_vec[4]
                    break
    return ObjectiveValue

def readSurfaceSensitivity (surface_sens_file_name, selected_marker):
    # * surfaceSensitivityFileName is a string of
    #       the file name of adjoint objective surface sensitivity
    #       usually named as "SURFACE_ADJOINT.csv"
    # * selected_marker is a string of the name of the selected marker
    # * function output a dict of selected marker with a list of all
    #       corresponding marker surface objective sensitivities in x,y,z
    surfaceSens = ()
    with open(surface_sens_file_name) as suface_adj:
        lines = suface_adj.readlines()
        for l in lines[1:]: #start extracting data from 2nd to end
            if len(l) <= 90:
                # remove empty lines
                continue
            else:
                temp_line_conv_vec = l.split(',')
                #ptID = temp_line_conv_vec[0]
                sensitivity_x = temp_line_conv_vec[14]
                sensitivity_y = temp_line_conv_vec[15]
                sensitivity_z = temp_line_conv_vec[16]
            currentTuple = (sensitivity_x, sensitivity_y, sensitivity_z)
            surfaceSens = (*surfaceSens, currentTuple)
    surfaceSensdict = {selected_marker: list(surfaceSens)}
    return surfaceSensdict

def readmarker_index (surfaceSensitivityFileName):
    # * surfaceSensitivityFileName is a string of
    #       the file name of adjoint objective surface sensitivity
    #       usually named as "SURFACE_ADJOINT.csv"
    # * function output a list of index corresponding to first colum of the csv
    marker_index = ()
    with open(surfaceSensitivityFileName) as suface_adj:
        lines = suface_adj.readlines()
        for l in lines[1:]: #start extracting data from 2nd to end
            if len(l) <= 90: # remove empty lines
                continue
            else:
                temp_line_conv_vec = l.split(',')
                #ptID = temp_line_conv_vec[0]
                markerIdx = temp_line_conv_vec[0]
            marker_index = (*marker_index, markerIdx)

    return marker_index


def euclideanNorm(data):
    # return a tuple of nested tuple of 3d cell normal vector to obtain cell area
    result = []
    for vector in data:
        x, y, z = vector
        norm = np.sqrt(x**2 + y**2 + z**2)
        result.append(norm)
    return tuple(result)


def multipleElementwise(tuple1, tuple2, marker_index):
    # return a tuple of the multiplication of these two tuple
    result = []
    for x, y in zip(tuple1[marker_index], tuple2):
        tmp_result = []
        for str in x:
            float_num = float(str)
            tmp_result.append(float_num * y)
        result.append( tmp_result )
    return_dict = {marker_index: result}
    return return_dict
