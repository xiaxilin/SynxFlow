#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Xiaodong Ming, Xilin Xia

"""
InputModel
===========

The module is based on the following assumptions:

    * Input DEM is a regular DEM file
    * its map unit is meter
    * its cellsize is the same in both x and y direction
    * its reference position is on the lower left corner of the southwest cell
    * All the other grid-based input files must be consistent with the DEM file

To do:

    * generate input (including sub-folder mesh and field) and output folders
    * generate mesh file (DEM.txt) and field files
    * divide model domain into small sections if multiple GPU is used

-----------------

"""

# Structure:

#    class InputModel

#        - Initialize an object: __init__

#        - set model parameters: 

#            set_boundary_condition, set_parameter, set_rainfall, 

#            set_gauges_position, set_case_folder, set_runtime, set_device_no,

#            add_user_defined_parameter, set_num_of_sections

__author__ = "Xiaodong Ming"
import os
import shutil
import copy
import warnings
import numpy as np
import glob
from datetime import datetime
from netCDF4 import Dataset
from .Raster import Raster
from .Boundary import Boundary
from .Rainfall import Rainfall
from .Landcover import Landcover
from .Summary import Summary
from .spatial_analysis import sub2map
from . import indep_functions as indep_f
#%% definition of class InputModel
class InputModel:
    """To define input files for a flood model case

    Read data, process data, write input files, and save data of a model case.

    Attributes:
        num_of_sections: (scalar) number of GPUs to run the model
        shape: shape of the DEM array
        header: (dict) header of the DEM grid
        __attributes_default: (dict) default model attribute names and values
        attributes: (dict) model attribute names and values
        times:  (list/numpy array) of four values reprenting model run time in
                seconds: start, end, output interval, backup interval
        device_no: (int) the gpu device id(s) to run model
        param_per_landcover: dict, argument to set grid parameters using the
            Landcover object. Keys are grid parameter names. Each Value is a 
            dict with three keys: param_value, land_value, default_value. 
            Refer to Landcover for more details.
        Sections: a list of objects of child-class InputModelSub
        Boundary: A boundary object for boundary conditions
        DEM: a Raster object to provide DEM data [alias: Raster].
        Rainfall: a Rainfall object to process rainfall data
        Landcover: a Landcover object to provide landcover data for setting
            gridded parameters        
        Summary: a Summary object to record model information
    """

#%%****************************************************************************
#***************************initialize an object*******************************
    # default parameters
    __attributes_default = {'h0':0, 'hU0x':0, 'hU0y':0,
                          'precipitation':0,
                          'precipitation_mask':0,
                          'precipitation_source':np.array([[0, 0], [1, 0]]),
                          'manning':0.035,
                          'sewer_sink':0,
                          'cumulative_depth':0, 'hydraulic_conductivity':0,
                          'capillary_head':0, 'water_content_diff':0,
                          'gauges_pos':np.array([[0, 0], [1, 1]])}
    _file_tag_list = ['z', 'h', 'hU', 'precipitation',
                      'manning', 'sewer_sink',
                      'cumulative_depth', 'hydraulic_conductivity',
                      'capillary_head', 'water_content_diff',
                      'precipitation_mask', 'precipitation_source',
                      'boundary_condition', 'gauges_pos']
    __grid_files = ['z', 'h', 'hU', 'precipitation_mask',
                  'manning', 'sewer_sink', 'precipitation',
                  'cumulative_depth', 'hydraulic_conductivity',
                  'capillary_head', 'water_content_diff']
    gridded_parameter_keys = ['h0', 'hU0x', 'hU0y', 'manning', 'sewer_sink',
                           'cumulative_depth', 'hydraulic_conductivity',
                           'capillary_head', 'water_content_diff']
    def __init__(self, dem_data=None, num_of_sections=1, case_folder=None,
                 data_path=None):
        """Initialise the InputModel object

        Args:
            dem_data: (Raster object) or (str) provides file name of the DEM data
            data_folder: a path contain at least a DEM file named as 'DEM' with a
                suffix .gz|.asc|.tif. 'landcover' and 'rain_mask' can also be read if
                these files were given with one of the three suffix
        """
        self.attributes = InputModel.__attributes_default.copy()
        self.num_of_sections = num_of_sections
        self.birthday = datetime.now()
        if case_folder is None:
            case_folder = os.getcwd()
        self._case_folder = case_folder
        if data_path is not None:
            self.data_path = data_path
            self._setup_by_files(data_path)
        else:
            if type(dem_data) is str:
                self.DEM = Raster(dem_data) # create Raster object
            elif type(dem_data) is Raster:
                self.DEM = dem_data
        self.Raster = self.DEM
        self.shape = self.DEM.shape
        self.header = self.DEM.header
        if not hasattr(self, 'Rainfall'):
            self.set_rainfall()
        _ = self.set_runtime()
        # get row and col index of all cells on DEM grid
        self._get_cell_subs()  # add _valid_cell_subs and _outline_cell_subs
        # divide model domain to several sections if it is not a sub section
        # each section contains a "model_IO_class.InputModelSub" object
        if isinstance(self, InputModelSub):
            pass
        else:
            self.__divide_grid()
        self.set_case_folder() # set data_folders
        self.set_device_no() # set the device number
        self.set_boundary_condition(outline_boundary='fall')
        self.set_gauges_position()
        self._initialize_summary_obj()# initialize a Model Summary object
        
    def __str__(self):
        """To show object summary information when it is called in console
        """
        self.Summary.display()
        time_str = self.birthday.strftime('%Y-%m-%d %H:%M:%S')
        return  self.__class__.__name__+' object created on '+ time_str
#    __repr__ = __str__
#******************************************************************************
#%%**********************Get object attributes**************************************
    def get_case_folder(self):
        """Get the directory of the case.

        Return: 
            string: The directory of the case
        """
        return self._case_folder

    def get_data_folders(self):
        """Get the directories and subdirectories of the input and output data.
        
        Return:
            string: a dict with four keys ('input', 'output', 'mesh', 'field')
                providing their absolute paths
        """
        return self._data_folders
    
    def get_input_filenames(self):
        """Get the input file names.
        
        Return:
            list: a list of string giving all input filenames without 
        """
        return self._file_tag_list

#%%**********************Setup the object**************************************
    def set_boundary_condition(self, boundary_list=None,
                               outline_boundary='fall'):
        """To create a Boundary object for boundary conditions
        
        Args:
            boundary_list: (list of dicts), each dict contain keys (polyPoints,
                type, h, hU) to define a IO boundary's position, type, and
                Input-Output (IO) sources timeseries. Keys including:

                1.polyPoints is a numpy array giving X(1st col) and Y(2nd col)
                    coordinates of points to define the position of a boundary.
                    An empty polyPoints means outline boundary.

                2.type: 'open'(flow out flatly), 'rigid'(no outlet),
                        'fall'(water flow out like a fall)

                3.h: a two-col numpy array. The 1st col is time(s). The 2nd col
                     is water depth(m)

                4.hU: a two-col numpy array. The 1st col is time(s). The 2nd
                    col is discharge(m3/s) or a three-col numpy array, the 2nd
                    col and the 3rd col are velocities(m/s) in x and y
                    direction, respectively.
            outline_boundary: (str) 'open'|'rigid', default outline boundary is
                open and both h and hU are set as zero

        """
        if boundary_list is None and hasattr(self, 'Boundary'):
            boundary_list = self.Boundary.boundary_list
        bound_obj = Boundary(boundary_list, outline_boundary)
        valid_subs = self._valid_cell_subs
        outline_subs = self._outline_cell_subs
        if not isinstance(self, InputModelSub):
            dem_header = self.header
        # add the subsripts and id of boundary cells on the domain grid
            bound_obj._fetch_boundary_cells(valid_subs,
                                            outline_subs, dem_header)
            bound_obj._convert_flow2velocity(self.DEM)
        self.Boundary = bound_obj
        if hasattr(self, 'Sections'):
            bound_obj._divide_domain(self)       
        if hasattr(self, 'Summary'):
            self.Summary.set_boundary_attr(self.Boundary)

    def set_initial_condition(self, parameter_name, parameter_value):
        """Set initial condition for h0, hU0x, hU0y
        
        Args:
            parameter_name: (str) h0, hU0x, hU0y
            parameter_value: scalar or numpy array with the same size of DEM.
        """
        if parameter_name not in ['h0', 'hU0x', 'hU0y']:
            raise ValueError('Parameter is not recognized: '+parameter_name)
        if type(parameter_value) is np.ndarray:
            if parameter_value.shape != self.shape:
                raise ValueError('The array of the parameter '
                                 'value should have the same '
                                 'shape with the DEM array')
        elif np.isscalar(parameter_value) is False:
            raise ValueError('The parameter value must be either '
                             'a scalar or an numpy array')
        self.attributes[parameter_name] = parameter_value
        self.Summary.set_initial_attr(**{parameter_name:parameter_value})
    
    def set_grid_parameter(self, **kwargs):
        """ Set grid parameter with Landcover object as name=value

        Args:
            kwargs: Keyword Arguments Specified by a Dictionary
                keyword: name, from grid_parameter_keys

                value:  1. scalar, a uniform parameter value
                        2. array, gridded parameter value with the same shape
                            of the DEM array
                        3. dict, contains param_value, land_value,  and 
                            default_value=0
        """
        if not hasattr(self, 'param_per_landcover'):
            # save arguments to call Landcover.to_grid_parameter()
            self.param_per_landcover = {}
        for keyword, value in kwargs.items():
            if keyword not in self.gridded_parameter_keys:
                raise ValueError('Parameter is not recognized: '+keyword)
            if type(value) is np.ndarray:
                if value.shape != self.shape:
                    raise ValueError('The array of the parameter '+keyword+
                           ' should have the same shape with the DEM array')
                else:
                    self.attributes[keyword] = value
            elif np.isscalar(value):
                self.attributes[keyword] = value
            elif type(value) is dict:
                self.param_per_landcover[keyword] = value
                value_array = self.Landcover.to_grid_parameter(**value)
                self.attributes[keyword] = value_array
            else:
                raise ValueError(keyword+' must be a scalar, array or dict')
            self.Summary.set_params_attr(**{keyword:value})

    def set_rainfall(self, rain_mask=None, rain_source=None):
        """ Set rainfall mask and rainfall source

        Args:
            rain_mask: str [filename of a Raster endswith .gz/asc/tif]
                    numpy int array with the same shape with DEM array
                    a Raster object
            rain_source: str [filename of a csv file for rainfall source data]
                    numpy array the 1st column is time in seconds, 2nd to
                    the end columns are rainfall rates in m/s.
        """
        if not hasattr(self, 'Rainfall'):
            self.Rainfall = Rainfall(self.attributes['precipitation_mask'], 
                                     self.attributes['precipitation_source'],
                                     dem_ras=self.DEM)
        # set rain mask
        if (rain_mask is None) & (rain_source is None): 
            # use default or pre-defined value
            rain_mask = self.attributes['precipitation_mask']
            rain_source = self.attributes['precipitation_source']
            self.Rainfall = Rainfall(rain_mask, rain_source, dem_ras=self.DEM)
        if rain_mask is not None:
            self.Rainfall.set_mask(rain_mask, dem_ras=self.DEM)
        if rain_source is not None:
            self.Rainfall.set_source(rain_source)
        if hasattr(self, 'Summary'):
            self.Summary.rain_attr = self.Rainfall.attrs

    def set_gauges_position(self, gauges_pos=None):
        """Set coordinates of monitoring gauges

        Args:
            gauges_pos: (numpy array) the 1st column is X coordinates, 
                2nd column is the Y coordinates
        """
        if gauges_pos is None:
            gauges_pos = self.attributes['gauges_pos']
        if type(gauges_pos) is list:
            gauges_pos = np.array(gauges_pos)
        if gauges_pos.shape[1] != 2:
            raise ValueError('The gauges_pos arraymust have two columns')
        self.attributes['gauges_pos'] = gauges_pos
#        self.Summary.add_param_infor('gauges_pos', gauges_pos)
        # for multi_GPU, divide gauges based on the extent of each section
        if hasattr(self, 'Sections'):
            pos_X = gauges_pos[:,0]
            pos_Y = gauges_pos[:,1]
            for obj_section in self.Sections:
                extent = obj_section.DEM.extent
                ind_x = np.logical_and(pos_X >= extent[0], pos_X <= extent[1])
                ind_y = np.logical_and(pos_Y >= extent[2], pos_Y <= extent[3])
                ind = np.where(np.logical_and(ind_x, ind_y))
                ind = ind[0]
                obj_section.attributes['gauges_pos'] = gauges_pos[ind,:]
                obj_section.attributes['gauges_ind'] = ind
        if hasattr(self, 'Summary'):
            self.Summary.set_model_attr(gauges_pos=gauges_pos)

    def set_case_folder(self, new_folder=None, make_dir=False):
        """ Initialize, renew, or create case and data folders
        
        Args:
            new_folder: (str) renew case and data folder if it is given
            make_dir: True|False create folders if it is True
        """
        # to change case_folder
        if new_folder is None:
            new_folder = self._case_folder
        self._case_folder = new_folder
        self._data_folders = indep_f._create_io_folders(self._case_folder,
                                               make_dir)
        # for multiple GPUs
        if hasattr(self, 'Sections'):
            for obj in self.Sections:
                sub_case_folder = os.path.join(new_folder, str(obj.section_id))
                obj.set_case_folder(sub_case_folder)                        
        if hasattr(self, 'Summary'):
            self.Summary.set_model_attr(case_folder=self._case_folder)

    def set_runtime(self, runtime=None):
        """set runtime of the model

        Args:
            runtime: a list of four values representing start, end, output 
                interval and backup interval respectively
        """
        if runtime is None:
            runtime = [0, 3600, 3600, 3600]
        runtime = np.array(runtime)
        self.times = runtime
        runtime_str = ('{0}-start, {1}-end, {2}-output interval, '
                       '{3}-backup interval')
        runtime_str = runtime_str.format(*runtime)
        if hasattr(self, 'Summary'):
#            self.Summary.add_items('Runtime(s)', runtime_str)
            self.Summary.set_model_attr(run_time=self.times)
        return runtime_str

    def set_device_no(self, device_no=None):
        """set device no of the model

        Args:
            device_no: int or a list of int corresponding to the number of 
                sections 
        """
        if device_no is None:
            device_no = np.arange(self.num_of_sections)
        device_no = np.array(device_no)
        self.device_no = device_no
    
    def set_landcover(self, landcover_data):
        """ Set Landcover object with a Raster object or file

        Args:
            landcover_data: (string or Raster) A Raster file or its file name 
                of land cover data
        """
        self.Landcover = Landcover(landcover_data, self.DEM)

    def add_user_defined_parameter(self, param_name, param_value):
        """ Add a grid-based user-defined parameter to the model

        Args:
            param_name: (str) name the parameter and the input file name too
            param_value: (scalar) or (numpy arary) with the same shape of the 
                DEM array
        """
        if param_name not in InputModel.__grid_files:
            InputModel.__grid_files.append(param_name)
            InputModel._file_tag_list.append(param_name)
        self.attributes[param_name] = param_value
        print(param_name+ ' is added to the InputModel object')
    
    def set_num_of_sections(self, num_of_sections):
        """ set the number of divided sections to run a case

        It can transfer single-gpu to multi-gpu or the opposite way

        Args:
            num_of_sections: (int) number of domains
        """
        self.num_of_sections = num_of_sections
        self.set_device_no()
        if num_of_sections==1: # to single GPU
            if hasattr(self, 'Sections'):
                del self.Sections
        else: # to multiple GPU
            self.__divide_grid()
            outline_boundary = self.Boundary.outline_boundary
            self.set_boundary_condition(outline_boundary=outline_boundary)
            self.set_gauges_position()
            self.Boundary._divide_domain(self)
        self.set_case_folder(self._case_folder)
        self.birthday = datetime.now()
        self.Summary.set_model_attr(num_GPU=self.num_of_sections)
#        time_str = self.birthday.strftime('%Y-%m-%d %H:%M:%S')


#%%****************************************************************************
#************************Write input files*************************************
    def write_input_files(self, file_tag=None):
        """ Write input files

        To classify the input files and call functions needed to write each
            input files
        
        Args:
            file_tag: a string or list of string giving the name(s) of input
                files without suffix
        """
        self._make_data_dirs()
        grid_files = InputModel.__grid_files
        if file_tag is None or file_tag == 'all': # write all files
            write_list = self._file_tag_list
            if self.num_of_sections > 1:
                self.write_halo_file()
            self.write_mesh_file()
            self.write_runtime_file()
            self.write_device_file()
        elif type(file_tag) is str:
            write_list = [file_tag]
        elif type(file_tag) is list:
            write_list = file_tag
        else:
            print(self._file_tag_list)
            raise ValueError(('file_tag should be a string or a list of string'
                             'from the above list'))
        for one_file in write_list:
            if one_file in grid_files: # grid-based files
                self.write_grid_files(one_file)
            elif one_file == 'boundary_condition':
                self.write_boundary_conditions()
            elif one_file == 'precipitation_source':
                self.write_rainfall_source()
            elif one_file == 'gauges_pos':
                self.write_gauges_position()
            else:
                raise ValueError(one_file+' is not recognized')

    def write_grid_files(self, file_tag, is_single_gpu=False):
        """Write grid-based files

        Public version for both single and multiple GPUs

        Args:
            file_tag: the pure name of a grid-based file
        """
        self._make_data_dirs()
        grid_files = InputModel.__grid_files
        if file_tag not in grid_files:
            raise ValueError(file_tag+' is not a grid-based file')
        if is_single_gpu or self.num_of_sections == 1:
            # write as single GPU even the num of sections is more than one
            self._write_grid_files_netcdf(file_tag, is_multi_gpu=False)
        else:
            self._write_grid_files_netcdf(file_tag, is_multi_gpu=True)
        readme_filename = os.path.join(self._case_folder,'readme.txt')
        self.Summary.to_json(readme_filename)
        print(file_tag+' created')

    def write_boundary_conditions(self):
        """ Write boundary condtion files

        if there are multiple domains, write in the first folder
            and copy to others
        """
        self._make_data_dirs()
        if self.num_of_sections > 1:  # multiple-GPU
            field_dir = self.Sections[0]._data_folders['field']
            file_names_list = self.__write_boundary_conditions(field_dir)
            self.__copy_to_all_sections(file_names_list)
        else:  # single-GPU
            field_dir = self._data_folders['field']
            self.__write_boundary_conditions(field_dir)
        readme_filename = os.path.join(self._case_folder,'readme.txt')
        self.Summary.to_json(readme_filename)
        print('boundary condition files created')

    def write_rainfall_source(self):
        """Write rainfall source data

        rainfall mask can be written by function write_grid_files
        """
        self._make_data_dirs()
        if hasattr(self, 'Rainfall'):
            rain_source = self.Rainfall.get_source_array()
        else:
            rain_source = self.attributes['precipitation_source']
        case_folder = self._case_folder
        num_of_sections = self.num_of_sections
        indep_f.write_rain_source(rain_source, case_folder, num_of_sections)
        readme_filename = os.path.join(self._case_folder,'readme.txt')
        self.Summary.to_json(readme_filename)

    def write_gauges_position(self, gauges_pos=None):
        """ Write the gauges position file

        Public version for both single and multiple GPUs

        Args:
            gauge_pos: (2-coloumn numpy array) positions of gauges 
        """
        self._make_data_dirs()
        if gauges_pos is not None:
            self.set_gauges_position(np.array(gauges_pos))
        if self.num_of_sections > 1:  # multiple-GPU
            for obj_section in self.Sections:
                field_dir = obj_section._data_folders['field']
                obj_section.__write_gauge_ind(field_dir)
                obj_section.__write_gauge_pos(field_dir)
        else:  # single-GPU
            field_dir = self._data_folders['field']
            self.__write_gauge_pos(field_dir)
        readme_filename = os.path.join(self._case_folder,'readme.txt')
        self.Summary.to_json(readme_filename)
        print('gauges_pos.dat created')

    def write_halo_file(self):
        """ Write overlayed cell IDs
        """
        num_section = self.num_of_sections
        case_folder = self._case_folder
        file_name = os.path.join(case_folder, 'halo.dat')
        with open(file_name, 'w') as file2write:
            file2write.write("No. of Domains\n")
            file2write.write("%d\n" % num_section)
            for obj_section in self.Sections:
                file2write.write("#%d\n" % obj_section.section_id)
                overlayed_id = obj_section.overlayed_id
                for key in ['bottom_low', 'bottom_high',
                            'top_high', 'top_low']:
                    if key in overlayed_id.keys():
                        line_ids = overlayed_id[key]
                        line_ids = np.reshape(line_ids, (1, line_ids.size))
                        np.savetxt(file2write,
                                   line_ids, fmt='%d', delimiter=' ')
                    else:
                        file2write.write(' \n')
        print('halo.dat created')

    def write_mesh_file(self, is_single_gpu=False):
        """ Write mesh file DEM.txt, compatoble for both single and multiple
            GPU model
        """
        self._make_data_dirs()
        if is_single_gpu is True or self.num_of_sections == 1:
            file_name = os.path.join(self._data_folders['mesh'],
                                     'DEM.txt')
            self.DEM.write_asc(file_name)
        else:
            for obj_section in self.Sections:
                file_name = os.path.join(obj_section._data_folders['mesh'],
                                         'DEM.txt')
                obj_section.DEM.write_asc(file_name)
        readme_filename = os.path.join(self._case_folder,'readme.txt')
        self.Summary.to_json(readme_filename)
    
    def write_runtime_file(self, time_values=None):
        """ write times_setup.dat file
        """
        if time_values is None:
            time_values = self.times
        indep_f.write_times_setup(self._case_folder, self.num_of_sections,
                                  time_values)
        self.Summary.to_json(self._case_folder+'/readme.txt')
    
    def write_device_file(self, device_no=None):
        """Create device_setup.dat for choosing GPU number to run the model
        """
        if device_no is None:
            device_no = self.device_no
        indep_f.write_device_setup(self._case_folder, self.num_of_sections,
                                   device_no)

    def write_landslide_config(self, rheology_type = 1, rheology_params = [0.5, 0, 2000], gravity_correction_type = 0, curvature_on = False, filter_mass_flux = False):
        """
        Write the configuration values only to a file, based on provided parameters.
        :param filename: The name of the file to write to.
        :param rheology_type: Type of rheology (1, 2, or 3)
        :param rheology_params: Parameters specific to the rheology type
        :param gravity_correction_type: Type of gravity correction
        :param curvature_on: Boolean indicating if curvature is on
        :param filter_mass_flux: Boolean indicating if mass flux filtering is on
        """
        directory = os.path.join(self._case_folder, 'input')
        if not os.path.exists(directory):
            os.makedirs(directory)
        filename = os.path.join(directory,'setup.conf')
        with open(filename, 'w') as file:
            print("[Type of rheology]", file=file)
            print(rheology_type, file=file)
            print("[Parameter values]", file=file)
            if rheology_type == 1:
                _miu, _cohesion, _rho = rheology_params
                print(f"{_miu} {_cohesion} {_rho}", file=file)
            elif rheology_type == 2:
                _miu1, _miu2, _L, _beta = rheology_params
                print(f"{_miu1} {_miu2} {_L} {_beta}", file=file)
            elif rheology_type == 3:
                _miu1, _miu2, _U = rheology_params
                print(f"{_miu1} {_miu2} {_U}", file=file)

            print("[Gravity correction type]", file=file)
            print(gravity_correction_type, file=file)
            print("[Consider curvature effect or not]", file=file)
            print(str(curvature_on).lower(), file=file)
            print("[Filter mass flux or not]", file=file)
            print(str(filter_mass_flux).lower(), file=file)

    def write_debris_config(self, rhoW=1000, rhoS=2650, particleD=0.00308, phi=0.4, criticalSlope=1.0, alpha=1.0, beta=0.01):
        """
        Writes given parameters to a parameters.dat file in the specified case folder.
        
        Parameters:
        - case_folder (str): Path to the case folder where the parameters file will be written.
        - rhoW (float): Water density in kg/m3. Default is 1000.
        - rhoS (float): Solids density in kg/m3. Default is 2650.
        - particleD (float): Particle diameter in m. Default is 0.00308.
        - phi (float): Porosity. Default is 0.4.
        - criticalSlope (float): Critical slope. Default is 1.0.
        - alpha (float): Alpha parameter. Default is 1.0.
        - beta (float): Beta parameter. Default is 0.01.
        """
        directory = os.path.join(self._case_folder, 'input')
        if not os.path.exists(directory):
            os.makedirs(directory)
        filename = os.path.join(directory,'parameters.dat')
        
        with open(filename, 'w') as f:
            f.write(f"{rhoW} {rhoS} {particleD} {phi} {criticalSlope} {alpha} {beta}")

    def load_backup(self, time):
        """
        Prepare running simulations from backup time points, after calling this function, the simulation can start from that point by calling run

        Parameters:
        - time (float): Time since when the simulation can start again.
        """
        # Construct the search pattern for the given time
        time_str = str(time)
        out_folder = os.path.join(self._case_folder, 'output')
        search_pattern = os.path.join(out_folder,  f"*backup__{time_str}.nc")
        
        # Find all matching files
        matching_files = glob.glob(search_pattern)
        
        if not matching_files:
            print(f"No files found with pattern *backup__{time_str}.nc in {out_folder}")
            return
        
        in_folder = os.path.join(self._case_folder, 'input', 'field')
        # Ensure the input folder exists
        if not os.path.exists(in_folder):
            print("Input folder did not exist, make sure generating them first")
            os.makedirs(in_folder)
        
        # Copy each matching file to the input folder
        for file_path in matching_files:
            try:
                file_name = os.path.basename(file_path)
                new_file_name = file_name.split("_backup__")[0] + ".nc"
                new_file_path = os.path.join(in_folder, new_file_name)
                
                shutil.copy(file_path, new_file_path)
                print(f"Copied {file_path} to {new_file_path}")
            except Exception as e:
                print(f"Error copying {file_path}: {e}")

        new_runtime = np.array([time] + self.times[1:].tolist())
        self.write_runtime_file(new_runtime)


    def save_object(self, file_name):
        """ Save object as a pickle file
        """
        indep_f.save_as_dict(self, file_name)

#%%****************************************************************************
#******************************* Visualization ********************************
    def domain_show(self, title='Domain Map', show_split=False, **kwargs):
        """Show domain map of the object
        
        Args:
            title: string to set the figure titile, 'Domain Map' is the defualt
            show_split: logical, show split lines for multi-gpu division 
        
        Retrun:
            fig: figure handle
            ax: axis handle
        """
        obj_dem = copy.deepcopy(self.DEM)
        fig, ax = obj_dem.mapshow(title=title, cax_str='DEM(m)', **kwargs)
        cell_subs = self.Boundary.cell_subs
        legends = []
        if show_split:            
            if hasattr(self, 'Sections'):
                x_overlayed = []
                y_overlayed = []
                for obj_sub in self.Sections:
                    overlayed_subs = obj_sub.overlayed_cell_subs_global
                    rows = overlayed_subs[0]
                    cols = overlayed_subs[1]
                    X, Y = sub2map(rows, cols, self.DEM.header)
                    x_overlayed = np.append(x_overlayed, X)
                    y_overlayed = np.append(y_overlayed, Y)
                ax.plot(x_overlayed, y_overlayed, '.k')
                legends.append('Splitting cells')
            else:
                warnings.warn("show_split is only for multi-gpu model")
        num = 0
        for cell_sub in cell_subs:
            rows = cell_sub[0]
            cols = cell_sub[1]
            X, Y = sub2map(rows, cols, self.DEM.header)
            ax.plot(X, Y, '.')
            if num==0:
                legends.append('Outline cells')
            else:
                legends.append('Boundary '+str(num))
            num = num+1
        
        ax.legend(legends, edgecolor=None, facecolor=None, loc='best',
                  fontsize='x-small')
        return fig, ax
    
    def plot_rainfall_map(self, figname=None, method='sum', **kw):
        """plot rainfall map within model domain
        """
        fig, ax = self.Rainfall.plot_rainfall_map(method='sum', **kw)
        cell_subs = self._outline_cell_subs
        rows = cell_subs[0]
        cols = cell_subs[1]
        X, Y = sub2map(rows, cols, self.DEM.header)
        ax.plot(X, Y, '.k')
        return fig, ax 

    def plot_rainfall_curve(self, start_date=None, method='mean', **kw):
        """ Plot time series of average rainfall rate inside the model domain

        Args:
            start_date: a datetime object to give the initial datetime of rain
            method: 'mean'|'max','min','mean', method to calculate gridded
                rainfall over the model domain
        """
        fig, ax = self.Rainfall.plot_time_series(method, **kw)
        return fig, ax

#%%****************************************************************************
#*************************** Protected methods ********************************
    def _get_cell_subs(self, dem_array=None):
        """ To get valid_cell_subs and outline_cell_subs for the object

        To get the subscripts of each valid cell on grid
        """
        if dem_array is None:
            dem_array = self.DEM.array
        valid_id, outline_id = indep_f._get_cell_id_array(dem_array)
        subs = np.where(~np.isnan(valid_id))
        id_vector = valid_id[subs]
        # sort the subscripts according to cell id values
        sorted_vectors = np.c_[id_vector, subs[0], subs[1]]
        sorted_vectors = sorted_vectors[sorted_vectors[:, 0].argsort()]
        self._valid_cell_subs = (sorted_vectors[:, 1].astype('int32'),
                                 sorted_vectors[:, 2].astype('int32'))
        subs = np.where(outline_id == 0) # outline boundary cell
        outline_id_vect = outline_id[subs]
        sorted_array = np.c_[outline_id_vect, subs[0], subs[1]]
        self._outline_cell_subs = (sorted_array[:, 1].astype('int32'),
                                   sorted_array[:, 2].astype('int32'))
    def __divide_grid(self):
        """Divide DEM grid to sub grids

        Create objects based on sub-class InputModelSub
        """
        if isinstance(self, InputModelSub):
            return 0  # do not divide InputModelSub objects, return a number
        else:
            if self.num_of_sections == 1:
                return 1 # do not divide if num_of_sections is 1
        num_of_sections = self.num_of_sections
        dem_header = self.header
        # subscripts of the split row [0, 1,...] from bottom to top
        split_rows = indep_f._get_split_rows(self.DEM.array, num_of_sections)
        array_local, header_local = \
            indep_f._split_array_by_rows(self.DEM.array, dem_header,
                                         split_rows)
        # to receive InputModelSub objects for sections
        Sections = []
        section_sequence = np.arange(num_of_sections)
        header_global = dem_header
        for i in section_sequence:  # from bottom to top
            case_folder = os.path.join(self._case_folder, str(i))
            # create a sub object of InputModel
            sub_model = InputModelSub(array_local[i], header_local[i],
                                        case_folder, num_of_sections)
            # get valid_cell_subs on the global grid
            valid_cell_subs = sub_model._valid_cell_subs
            valid_subs_global = \
                 indep_f._cell_subs_convertor(valid_cell_subs, header_global,
                                      header_local[i], to_global=True)
            sub_model.valid_subs_global = valid_subs_global
            # record section sequence number
#            sub_model.section_id = i
            #get overlayed_id (top two rows and bottom two rows)
            top_h = np.where(valid_cell_subs[0] == 0)
            top_l = np.where(valid_cell_subs[0] == 1)
            bottom_h = np.where(
                valid_cell_subs[0] == valid_cell_subs[0].max()-1)
            bottom_l = np.where(valid_cell_subs[0] == valid_cell_subs[0].max())
            if i == 0: # the bottom section
                overlayed_id = {'top_high':top_h[0], 'top_low':top_l[0]}
            elif i == self.num_of_sections-1: # the top section
                overlayed_id = {'bottom_low':bottom_l[0],
                                'bottom_high':bottom_h[0]}
            else:
                overlayed_id = {'top_high':top_h[0], 'top_low':top_l[0],
                                'bottom_high':bottom_h[0],
                                'bottom_low':bottom_l[0]}
            sub_model.overlayed_id = overlayed_id
            all_ids = list(overlayed_id.values())
            all_ids = np.concatenate(all_ids).ravel()
            all_ids.sort()
            overlayed_cell_subs_global = (valid_subs_global[0][all_ids],
                                          valid_subs_global[1][all_ids])
            sub_model.overlayed_cell_subs_global = overlayed_cell_subs_global
            Sections.append(sub_model)
        # reset global var section_id of InputModelSub
        InputModelSub.section_id = 0
        self.Sections = Sections
        self._initialize_summary_obj()# get a Model Summary object

    def _get_vector_value(self, attribute_name, is_multi_gpu=True,
                          add_initial_water=True):
        """ Generate a single vector for values in each grid cell sorted based
                on cell IDs
        
        Args:
            attribute_name: attribute names based on a grid

        Return:
            output_vector: a vector or list vectors, giving values on valid 
                global cells or sub-domain cells
        """
        # get grid value
        dem_shape = self.shape
        grid_values = np.zeros(dem_shape)
        # set grid value for the entire domain
        if attribute_name == 'z':
            grid_values = self.DEM.array
        elif attribute_name == 'h':
            grid_values = grid_values+self.attributes['h0']
            # traversal each boundary to add initial water
            if add_initial_water:
                wet_cell_subs = self.Boundary.cell_subs_wet_io
                if wet_cell_subs is not None:
                    grid_values[wet_cell_subs] = grid_values[wet_cell_subs]+0.001
        elif attribute_name == 'hU':
            grid_values0 = grid_values+self.attributes['hU0x']
            grid_values1 = grid_values+self.attributes['hU0y']
            grid_values = [grid_values0, grid_values1]
        elif attribute_name == 'precipitation_mask':
            if hasattr(self, 'Rainfall'):
                grid_values = self.Rainfall.get_mask_array()
            else:
                grid_values = grid_values+self.attributes[attribute_name]
        else:
            if hasattr(self, 'param_per_landcover'):
                arg_dicts = self.param_per_landcover
                if attribute_name in arg_dicts.keys():
                    arg_dict = arg_dicts[attribute_name]
                    grid_values = self.Landcover.to_grid_parameter(**arg_dict)
                else:
                    grid_values = grid_values+self.attributes[attribute_name]
            else:
                grid_values = grid_values+self.attributes[attribute_name]

        def grid_to_vect(grid_values, cell_subs):
            """ Convert grid values to 1 or 2 col vector values
            """
            if type(grid_values) is list:
                vector_value0 = grid_values[0][cell_subs]
                vector_value1 = grid_values[1][cell_subs]
                vector_value = np.c_[vector_value0, vector_value1]
            else:
                vector_value = grid_values[cell_subs]
            return vector_value
        #
        if is_multi_gpu: # generate vector value for multiple GPU
            output_vector = []
            for obj_section in self.Sections:
                cell_subs = obj_section.valid_subs_global
                vector_value = grid_to_vect(grid_values, cell_subs)
                output_vector.append(vector_value)
        else:
            output_vector = grid_to_vect(grid_values, self._valid_cell_subs)
        return output_vector

    def _get_boundary_id_code_array(self, file_tag='z'):
        """To generate a 4-col array of boundary cell id (0) and code (1~3)
        """
        bound_obj = self.Boundary
        output_array_list = []
        for ind_num in np.arange(bound_obj.num_of_bound):
            if file_tag == 'h':
                bound_code = bound_obj.data_table.h_code[ind_num]
            elif file_tag == 'hU':
                bound_code = bound_obj.data_table.hU_code[ind_num]
            elif file_tag == 'C':
                bound_code = bound_obj.data_table.C_code[ind_num]
            else:
                bound_code = np.array([[2, 0, 0]]) # shape (1, 3)
            if bound_code.ndim < 2:
                bound_code = np.reshape(bound_code, (1, bound_code.size))
            cell_id = bound_obj.cell_id[ind_num]
            if cell_id.size > 0:
                bound_code_array = np.repeat(bound_code, cell_id.size, axis=0)
                id_code_array = np.c_[cell_id, bound_code_array]
                output_array_list.append(id_code_array)
        # add overlayed cells with [4, 0, 0]
        # if it is a sub section object, there should be attributes:
        # overlayed_id, and section_id
        if hasattr(self, 'overlayed_id'):
            cell_id = []
            if 'top_high' in self.overlayed_id:
                cell_id += self.overlayed_id['top_high'].tolist()
            if 'bottom_low' in self.overlayed_id:
                cell_id += self.overlayed_id['bottom_low'].tolist()
            cell_id = np.array(cell_id)
            bound_code = np.array([[4, 0, 0]]) # shape (1, 3)
            bound_code_array = np.repeat(bound_code, cell_id.size, axis=0)
            id_code_array = np.c_[cell_id, bound_code_array]
            output_array_list.append(id_code_array)
        output_array = np.concatenate(output_array_list, axis=0)
        # when unique the output array according to cell id
        # keep the last occurrence rather than the default first occurrence
        output_array = np.flipud(output_array) # make the IO boundaries first
        _, ind = np.unique(output_array[:, 0], return_index=True)
        output_array = output_array[ind]
        return output_array

    def _initialize_summary_obj(self):
        """ Initialize the model summary object
        """
        summary_obj = Summary(self)
        self.Summary = summary_obj

    def _write_grid_files_netcdf(self, file_tag, is_multi_gpu=True):
        """ Write input files consistent with the DEM grid in NetCDF format

        Private function called by public function write_grid_files.
        """
        if is_multi_gpu is True:  # write for multi-GPU, use child object
            vector_value_list = self._get_vector_value(file_tag, is_multi_gpu)
            for obj_section in self.Sections:
                vector_value = vector_value_list[obj_section.section_id]
                cell_id = np.arange(vector_value.shape[0])
                cells_vect = np.c_[cell_id, vector_value]
                file_name = os.path.join(obj_section._data_folders['field'], file_tag + '.nc')

                if file_tag == 'precipitation_mask':
                    bounds_vect = None
                else:
                    bounds_vect = obj_section._get_boundary_id_code_array(file_tag)

                # Call the updated function for writing NetCDF
                indep_f._write_two_arrays_netcdf(file_name, cells_vect, bounds_vect)
        else:  # single GPU, use global object
            file_name = os.path.join(self._data_folders['field'],
                                     file_tag+'.nc')
            vector_value = self._get_vector_value(file_tag, is_multi_gpu=False)
            cell_id = np.arange(vector_value.shape[0])
            cells_vect = np.c_[cell_id, vector_value]
            if file_tag == 'precipitation_mask':
                bounds_vect = None
            else:
                bounds_vect = self._get_boundary_id_code_array(file_tag)
            indep_f._write_two_arrays_netcdf(file_name, cells_vect, bounds_vect)
        return None

    def _write_grid_files(self, file_tag, is_multi_gpu=True):
        """ Write input files consistent with the DEM grid

        Private function called by public function write_grid_files
        file_name: includes ['h','hU','precipitation_mask',
                             'manning','sewer_sink',
                             'cumulative_depth', 'hydraulic_conductivity',
                             'capillary_head', 'water_content_diff']
        """
        if is_multi_gpu is True:  # write for multi-GPU, use child object
            vector_value_list = self._get_vector_value(file_tag, is_multi_gpu)
            for obj_section in self.Sections:
                vector_value = vector_value_list[obj_section.section_id]
                cell_id = np.arange(vector_value.shape[0])
                cells_vect = np.c_[cell_id, vector_value]
                file_name = os.path.join(obj_section._data_folders['field'],
                                         file_tag+'.dat')
                if file_tag == 'precipitation_mask':
                    bounds_vect = None
                else:
                    bounds_vect = \
                        obj_section._get_boundary_id_code_array(file_tag)
                indep_f._write_two_arrays(file_name, cells_vect, bounds_vect)
        else:  # single GPU, use global object
            file_name = os.path.join(self._data_folders['field'],
                                     file_tag+'.dat')
            vector_value = self._get_vector_value(file_tag, is_multi_gpu=False)
            cell_id = np.arange(vector_value.shape[0])
            cells_vect = np.c_[cell_id, vector_value]
            if file_tag == 'precipitation_mask':
                bounds_vect = None
            else:
                bounds_vect = self._get_boundary_id_code_array(file_tag)
            indep_f._write_two_arrays(file_name, cells_vect, bounds_vect)
        return None

    def _make_data_dirs(self):
        """ Create folders in current device
        """
        if hasattr(self, 'Sections'):
            for obj_section in self.Sections:
                indep_f._create_io_folders(obj_section.get_case_folder(),
                                           make_dir=True)
        else:
            indep_f._create_io_folders(self._case_folder, make_dir=True)
    
    def _dict2grid(self, mask_dict):
        """Convert mask_dict to a grid array with the same shape of DEM
        """
        num_values = mask_dict['value'].size
        grid_array = np.zeros(self.DEM.shape).astype(mask_dict['value'].dtype)
        for i in np.arange(num_values):
            grid_array[mask_dict['index'][i]] = mask_dict['value'][i]
        return grid_array
        
#------------------------------------------------------------------------------
#*************** Private methods only for the parent class ********************
#------------------------------------------------------------------------------
    def __write_boundary_conditions(self, field_dir, file_tag='all'):
        """ Write boundary condition source files,if hU is given as flow
        timeseries, convert flow to hUx and hUy.
        Private function to call by public function write_boundary_conditions
        file_tag: 'h', 'hU','C', 'all'
        h_BC_[N].dat, hU_BC_[N].dat, C_BC_[N].dat
        if hU is given as flow timeseries, convert flow to hUx and hUy
        """
        obj_boundary = self.Boundary
        file_names_list = []
        fmt_h = ['%g', '%g']
        fmt_hu = ['%g', '%g', '%g']
        fmt_c = ['%g', '%g']
        # write h_BC_[N].dat
        if file_tag in ['all', 'h']:
            h_sources = obj_boundary.data_table['hSources']
            ind_num = 0
            for i in np.arange(obj_boundary.num_of_bound):
                h_source = h_sources[i]
                if h_source is not None:
                    file_name = os.path.join(field_dir,
                                             'h_BC_'+str(ind_num)+'.dat')
                    np.savetxt(file_name, h_source, fmt=fmt_h, delimiter=' ')
                    ind_num = ind_num+1
                    file_names_list.append(file_name)
        # write hU_BC_[N].dat
        if file_tag in ['all', 'hU']:
            hU_sources = obj_boundary.data_table['hUSources']
            ind_num = 0
            for i in np.arange(obj_boundary.num_of_bound):
                hU_source = hU_sources[i]
                if hU_source is not None:
                    file_name = os.path.join(field_dir,
                                             'hU_BC_'+str(ind_num)+'.dat')
                    np.savetxt(file_name, hU_source, fmt=fmt_hu, delimiter=' ')
                    ind_num = ind_num+1
                    file_names_list.append(file_name)
        # write C_BC_[N].dat
        if file_tag in ['all', 'C']:
            C_sources = obj_boundary.data_table['CSources']
            ind_num = 0
            for i in np.arange(obj_boundary.num_of_bound):
                C_source = C_sources[i]
                if C_source is not None:
                    file_name = os.path.join(field_dir,
                                             'C_BC_'+str(ind_num)+'.dat')
                    np.savetxt(file_name, C_source, fmt=fmt_c, delimiter=' ')
                    ind_num = ind_num+1
                    file_names_list.append(file_name)
        return file_names_list

    def __write_gauge_pos(self, file_folder):
        """write monitoring gauges
        gauges_pos.dat
        file_folder: folder to write file
        gauges_pos: 2-col numpy array of X and Y coordinates
        """
        gauges_pos = self.attributes['gauges_pos']
        file_name = os.path.join(file_folder, 'gauges_pos.dat')
        fmt = ['%g %g']
        fmt = '\n'.join(fmt*gauges_pos.shape[0])
        gauges_pos_str = fmt % tuple(gauges_pos.ravel())
        with open(file_name, 'w') as file2write:
            file2write.write(gauges_pos_str)
        return file_name

    def __write_gauge_ind(self, file_folder):
        """write monitoring gauges index for mult-GPU sections
        gauges_ind.dat
        file_folder: folder to write file
        gauges_ind: 1-col numpy array of index values
        """
        gauges_ind = self.attributes['gauges_ind']
        file_name = os.path.join(file_folder, 'gauges_ind.dat')
        fmt = ['%g']
        fmt = '\n'.join(fmt*gauges_ind.shape[0])
        gauges_ind_str = fmt % tuple(gauges_ind.ravel())
        with open(file_name, 'w') as file2write:
            file2write.write(gauges_ind_str)
        return file_name

    def __copy_to_all_sections(self, file_names):
        """ Copy files that are the same in each sections
        file_names: (str) files written in the first seciton [0]
        boundary source files: h_BC_[N].dat, hU_BC_[N].dat
        rainfall source files: precipitation_source_all.dat
        gauges position file: gauges_pos.dat
        """
        if type(file_names) is not list:
            file_names = [file_names]
        for i in np.arange(1, self.num_of_sections):
            field_dir = self.Sections[i]._data_folders['field']
            for file in file_names:
                shutil.copy2(file, field_dir)
    
    def _setup_by_files(self, data_path):
        """Read files to setup model input object
        DEM, landcover, rain_mask, endswith '.gz', '.asc', or '.tif'
        rain_source.csv
        """
        data_path = self.data_path
        ras_file = indep_f._check_raster_exist(os.path.join(data_path, 'DEM'))
        if ras_file is not None:
            self.DEM = Raster(ras_file)
            print(ras_file+ ' read')
        ras_file = indep_f._check_raster_exist(os.path.join(data_path,
                                                            'landcover'))        
        if ras_file is not None:
            self.Landcover = Landcover(ras_file, dem_ras=self.DEM)
            print(ras_file+ ' read')
        
        mask_file = indep_f._check_raster_exist(os.path.join(data_path,
                                                         'rain_mask'))
        if mask_file is None:
            mask_file = 0
        source_file = os.path.join(data_path, 'rain_source.csv')
        if not os.path.isfile(source_file):
            source_file = np.array([[0, 0], [1, 0]])
        self.Rainfall = Rainfall(mask_file, source_file, dem_ras=self.DEM)

#%%****************************************************************************
#************************sub-class definition**********************************
class InputModelSub(InputModel):
    """object for each section, child class of InputModel

    Attributes:
        sectionNO: the serial number of each section
        _valid_cell_subs: (tuple, int) two numpy array indicating rows and cols
            of valid cells on the local grid
        valid_cell_subsOnGlobal: (tuple, int) two numpy array indicating rows
            and cols of valid cells on the global grid
        shared_cells_id: 2-row shared Cells id on a local grid
        case_folder: input folder of each section
        _outline_cell_subs: (tuple, int) two numpy array indicating rows and 
            cols of valid cells on a local grid
    """
    section_id = 0
    def __init__(self, dem_array, header, case_folder, num_of_sections):
        self.section_id = InputModelSub.section_id
        InputModelSub.section_id = self.section_id+1
        dem_data = Raster(array=dem_array, header=header)
        super().__init__(dem_data, num_of_sections, case_folder)

#%%****************************************************************************
#********************************Static method*********************************
def load_input_object(filename):
    """load object from a dictionary and return as an InputModel object
    Args:
        filename: a string giving the object file name
    Return: 
        An object of InputModel
    """
    obj_dict = indep_f.load_object(filename)
    
    if 'DEM' in obj_dict:
        dem_dict = obj_dict['DEM']
        obj_dem = Raster(array=dem_dict['array'],
                                header=dem_dict['header'])
        obj_dict.pop('DEM')
    else:
        raise ValueError(filename+' has no key: DEM')
    obj_in = InputModel(dem_data=obj_dem,
                         num_of_sections=obj_dict['num_of_sections'],
                         case_folder=obj_dict['_case_folder'])
    
    if 'Landcover' in obj_dict:
        ld_dict = obj_dict['Landcover']
        mask_header = ld_dict['mask_header']
        mask_dict = ld_dict['mask_dict']
        array_shape = (mask_header['nrows'], mask_header['ncols'])
        mask_array = indep_f._dict2grid(mask_dict, array_shape)
        ras_landcover = Raster(array=mask_array, header=mask_header)
        obj_in.set_landcover(ras_landcover)
        obj_dict.pop('Landcover')
    
    if 'Rainfall' in obj_dict:
        rain_dict = obj_dict['Rainfall']
        mask_header = rain_dict['mask_header']
        mask_dict = rain_dict['mask_dict']
        array_shape = (mask_header['nrows'], mask_header['ncols'])
        mask_array = indep_f._dict2grid(mask_dict, array_shape)
        rain_mask = Raster(array=mask_array, header=mask_header)
        rain_source = np.c_[rain_dict['time_s'], rain_dict['rain_rate']]
        obj_in.set_rainfall(rain_mask, rain_source)
        obj_dict.pop('Rainfall')
    
    bound_dict = obj_dict['Boundary']
    for key, value in bound_dict.items():
        obj_in.Boundary.__dict__[key] = value
    obj_dict.pop('Boundary')
    
    summ_dict = obj_dict['Summary']
    for key, value in summ_dict.items():
        obj_in.Summary.__dict__[key] = value
    obj_dict.pop('Summary')
    
    for key, value in obj_dict.items():
        obj_in.__dict__[key] = value
    
    return obj_in

def main():
    print('Class to setup input data')

if __name__=='__main__':
    main()