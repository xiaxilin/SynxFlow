#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Wed Apr  1 23:46:54 2020
# Author: Xiaodong Ming
"""
OutputModel
============

To do:
    Process output files

----------------------

"""
import os
import datetime
import warnings
import numpy as np
import pandas as pd
from . import spatial_analysis as sp
from .Raster import Raster
from .Summary import Summary
from .indep_functions import save_as_dict, _create_io_folders, load_object
class OutputModel:
    """To read and analyze otuput files from a flood model

    Attributes:
        case_folder: (str) the absolute path of the case folder
        input_folder: (str|list of strings) the absolute path of the 
            input folder(s)
        output_folder: (str|list of strings) the absolute path of the 
            output folder(s)
        number_of_sections: (int) the number of subdomains of the model
        header: (dict or list of dict) provide the header information
        header_list: a list of sub headers [only for multi-gpu model]
        ref_datetime: 
        times_simu: (dataframe) with variable 'times' (simulated time in 
                    seconds) and 'date_times' if ref_datetime is defined.
        gauge_values_all: (dict) 'h', 'eta', 'hU', array with values of all 
            gauge positions, not time column
        gauge_values: gauge timeseries summarized from a series of gauge points
        grid_results: a list of Raster objects for gridded results

    """  
    def __init__(self, input_obj=None, case_folder=None,
                 num_of_sections=1, header_file_tag=None):
        """Initialize the object with a InputModel object or a case folder and
            the number of sections

        Args:
            header_file_tag: the output file to read grid header, e.g. 'h_0'
        """
        # pass argument values
        if input_obj is None:
            self.case_folder = case_folder
            self.num_of_sections = num_of_sections
            self._set_IO_folders()
            self._set_grid_header(asc_file=header_file_tag)
        elif type(input_obj) is dict: # load from a dict
            for key, value in input_obj.items():
                if key=='Summary':
                    self.Summary = Summary(value)
                else:
                    self.__dict__[key] = value
        elif hasattr(input_obj, 'DEM'):
            # get information from the input object
            # case_folder, num_of_sections, header
            self.case_folder = input_obj.get_case_folder()
            self.num_of_sections = input_obj.num_of_sections
            self.header = input_obj.header
            self.dem_array = input_obj.DEM.array
            self.output_folder = input_obj._data_folders['output']
            self.input_folder = input_obj._data_folders['input']
            self.Summary = input_obj.Summary
            if input_obj.num_of_sections>1:
                header_list = []
                output_folder = []
                input_folder = []
                for sub_obj in input_obj.Sections:
                    header_list.append(sub_obj.Raster.header)
                    output_folder.append(sub_obj._data_folders['output'])
                    input_folder.append(sub_obj._data_folders['input'])
                self.header_list = header_list  
                self.output_folder = output_folder
                self.input_folder = input_folder
        else:
            raise IOError('The first argument (input_obj) must be '
                          'a InputModel object')

    def read_gauges_file(self, file_tag='h', compressed=False):
        """ Read gauges files for time seires of values at the monitored gauges

        Args:
            file_tag: h, hU, eta, corresponding to h_gauges.dat, hU_gauges.dat,
                and eta_gauges.dat, respectively

        Return:
            gauges_pos, the coordinates of gauges within the model domain

            times, time in seconds

            values, gauge values corresponding to the gauges position


        """
        if self.num_of_sections==1:
            output_folder = self.output_folder
            gauge_output_file = os.path.join(output_folder,
                                             file_tag+'_gauges.dat')
            gauge_pos_file = os.path.join(self.input_folder,
                                          'field', 'gauges_pos.dat')
            if compressed:
                gauge_output_file = gauge_output_file+'.gz'
                gauge_pos_file = gauge_pos_file+'.gz'
            times, values = _read_one_gauge_file(gauge_output_file)
            gauges_pos = np.loadtxt(gauge_pos_file, dtype='float64', ndmin=2)
        else: # multi-GPU
            gauges_value_all, gauges_pos = _combine_gauges_data_via_ind(
                    self.case_folder, self.num_of_sections, file_tag)
            if file_tag == 'hU':
                values_x = np.array(gauges_value_all[0].iloc[:,:-1])
                values_y = np.array(gauges_value_all[1].iloc[:,:-1])
                values = np.array([values_x, values_y])
                times = np.array(gauges_value_all[0]['times'])
            else:
                values = np.array(gauges_value_all.iloc[:,:-1])
                times = np.array(gauges_value_all['times'])
        self.times_simu = pd.DataFrame({'times':times})
        if hasattr(self, 'ref_datetime'):
            times_delta = times.astype('timedelta64[s]')
            date_times = np.datetime64(self.ref_datetime)+times_delta                    
            self.times_simu['date_times'] = date_times
        if not hasattr(self, 'gauge_values_all'):
            self.gauge_values_all = {}
        self.gauge_values_all[file_tag] = values
        self.gauges_pos = gauges_pos
        return gauges_pos, times, values
    
    def read_grid_file(self, file_tag='h_0', grid_type='netCDF'):
        """Read asc grid files from output

        Return:
            Raster: a raster object of the grid

        """
        if grid_type == 'netCDF':
            file_tag = file_tag+'.nc'
        elif grid_type == 'compressed':
            file_tag = file_tag+'.gz'
        elif grid_type == 'asc':
            file_tag = file_tag+'.asc'
        if self.num_of_sections==1:
            file_name = os.path.join(self.output_folder, file_tag)
            if file_name.lower().endswith('.nc'):
                grid_array, _, _ = sp.ncgridread(file_name)
            else:
                grid_array, _, _ = sp.arcgridread(file_name)
        else: # multi-GPU
            grid_array = self._combine_multi_gpu_grid_data(file_tag)
        grid_obj = Raster(array=grid_array, header=self.header)
        return grid_obj
    
    def add_gauge_results(self, var_name, gauge_name='All', gauge_ind=None,  
                          compressed=False):
        """ add simulated value to the object gauge by gauge

        Args:
            var_name: 'h', 'hU', 'eta'
            gauge_name: 'All' add all gauges, then gauge_ind not needed

        """
        # read all gauge data in all positions
        if not hasattr(self, 'gauge_values_all'):
            self.read_gauges_file(var_name, compressed)
        else:
            if var_name not in self.gauge_values_all.keys():
                self.read_gauges_file(var_name, compressed)
        if gauge_ind is not None:
            gauge_ind = np.array(gauge_ind, ndmin=2)
            # add position data for a gauge
            values = self.gauge_values_all[var_name]+0
            values_pd = self.times_simu.copy()
            if not hasattr(self, 'gauge_values'):
                self.gauge_values = {}
            if var_name == 'h': # calculation method is min
                one_gauge_v = values[:, gauge_ind]
                one_gauge_v = one_gauge_v.max(axis=1)
                values_pd['values'] = one_gauge_v
            elif var_name == 'hU':
                one_gauge_v = values[:, :, gauge_ind]
                one_gauge_v = one_gauge_v.sum(axis=2)*self.header['cellsize']
                values_pd['values_x'] = one_gauge_v[0]
                values_pd['values_y'] = one_gauge_v[1]
            elif var_name == 'eta':
                if gauge_ind.size > 1:
                    raise ValueError('gauge_ind for eta must be a scalar')
                else:
                    one_gauge_v = values[:, gauge_ind[0]]
                    values_pd['values'] = one_gauge_v
            if gauge_name in self.gauge_values.keys():
                gauge_dict = self.gauge_values[gauge_name]
                gauge_dict[var_name] = values_pd
            else:
                gauge_dict = {var_name:values_pd}
            self.gauge_values[gauge_name] = gauge_dict
    
    def add_grid_results(self, result_names, compressed=False):
        """Read and return Raster object to attribute 'grid_results'

        Args:
            result_names: string or list of string, gives the name of grid file

        """
        if not hasattr(self, 'grid_results'):
            self.grid_results = {}
        if type(result_names) is not list: # for a list of files
            result_names = [result_names]
        for file_tag in result_names:
            grid_obj = self.read_grid_file(file_tag, compressed)
            self.grid_results[file_tag] = grid_obj.array

    def set_ref_datetime(self, date_time, str_format='%Y-%m-%d %H:%M:%S'):
        """Set the refernce datetime of the simulation
        """
        if type(date_time) is str:
            self.ref_datetime = datetime.datetime.strptime(date_time,
                                                           str_format)
        elif type(date_time) is datetime.datetime:
            self.ref_datetime = date_time
        else:
            raise IOError('date_time must be a datetime object or a string')
    
    def _combine_multi_gpu_grid_data(self, asc_file_name):
        """Combine multi-gpu grid files into a single file

        asc_file_name: string endswith '.asc'
        """
        header_global = self.header
        header_list = self.header_list
        output_folder = self.output_folder
        grid_shape = (header_global['nrows'], header_global['ncols'])
        array_global = np.zeros(grid_shape)
        for header0, folder0 in zip(header_list, output_folder):
            ind_top, ind_bottom = _header2row_numbers(header0, header_global)
            file_name = os.path.join(folder0, asc_file_name)
            if file_name.lower().endswith('.nc'):
                array_local, _, _ = sp.ncgridread(file_name)
            else:
                array_local, _, _ = sp.arcgridread(file_name)
            array_global[ind_top:ind_bottom+1,:] = array_local
        return array_global
    
    def _set_IO_folders(self):
        """ Set input and output folder/folders

        case_folder, num_of_sections are required
        """
        case_folder = self.case_folder
        num_of_sections = self.num_of_sections
        if num_of_sections == 1: # single gpu
            data_folders = _create_io_folders(case_folder)
            output_folder = data_folders['output']
            input_folder = data_folders['input']
        else: #multi-gpu model
            output_folder = []
            input_folder = []
            for i in range(num_of_sections):
                section_case_folder = os.path.join(case_folder, str(i))
                data_folders = _create_io_folders(section_case_folder)
                output_folder.append(data_folders['output'])  
                input_folder.append(data_folders['input'])
        self.output_folder = output_folder
        self.input_folder = input_folder
    
    def _set_grid_header(self, asc_file=None):
        """set header for the grid of the output object

        num_of_sections, input_folder, output_folder are required

        asc_file: the file to get grid header, default is DEM.txt

        """
        
        num_of_sections = self.num_of_sections
        output_folder = self.output_folder
        input_folder = self.input_folder
        warning_message = ('The grid header is not set, '
                           'try _set_grid_header(asc_file) to set a '
                           'header if it is required in your next steps')
        if num_of_sections == 1:
            if asc_file is None:
                file_name = os.path.join(input_folder, 'mesh', 'DEM.txt')
            else:
                file_name = os.path.join(output_folder, asc_file)
            if os.path.exists(file_name):
                self.header = sp.arc_header_read(file_name)
            else:
                warnings.warn(warning_message)
        else: #multi-gpu model
            headers = []
            header_set = True
            for i in np.arange(num_of_sections):
                if asc_file is None:
                    file_name = os.path.join(input_folder[i],'mesh', 'DEM.txt')
                else:
                    file_name = os.path.join(output_folder[i], asc_file)
                if os.path.exists(file_name):
                    header = sp.arc_header_read(file_name)
                    headers.append(header)
                    header_set = True*header_set
                else:
                    warnings.warn(warning_message)
                    header_set = False
            if header_set:
                self.header_list = headers
                self.header = _header_local2global(headers)
    
    def save_object(self, file_name):
        """Save the object to a pickle file
        """
#        save_object(self, file_name, compression=True)
        save_as_dict(self, file_name)
    
#%% =======================Supporting functions===============================
def load_output_object(filename):
    """load object from a dictionary and return as an OutputModel object
    
    Args:
        filename: a string giving the object file name
    Return: 
        An object of OutputModel
    """
    obj_dict = load_object(filename)
    if type(obj_dict) is OutputModel:
        obj_out = obj_dict
    elif type(obj_dict) is dict:
        obj_out = OutputModel(obj_dict)
    else:
        raise ValueError(filename+' should store either a dict or OutputModel'
                         ' object')
    return obj_out
    
def _combine_gauges_data_via_ind(case_folder, num_section, file_tag):
    """Combine gauges outputs from multi-gpu models according to gauges

    position index.

    gauge_ind.dat

    gauge_pos.dat in each subdoamin input/field folder

    """
    # read index
    ind_list =[]
    ind_max = 0
    gauges_pos_all = []
    for i in range(num_section):
        gauge_ind_file = os.path.join(case_folder, str(i), 'input', 'field',
                                      'gauges_ind.dat')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ind_1 = np.loadtxt(gauge_ind_file, dtype='int')
        if ind_1.size>0:
            ind_max = max(ind_max, ind_1.max())
            file_name = os.path.join(case_folder, str(i), 'input', 'field',
                                     'gauges_pos.dat')
            gauge_xy = np.loadtxt(file_name, dtype='float64', ndmin=2)
            gauges_pos_all.append(gauge_xy)
        ind_list.append(ind_1)
    gauges_pos_all = np.concatenate(gauges_pos_all, axis=0)
    num_gauges = ind_max+1
    if file_tag == 'hU':
        gauges_value_x = pd.DataFrame(columns=np.arange(num_gauges))
        gauges_value_y = pd.DataFrame(columns=np.arange(num_gauges))
        gauges_value = [gauges_value_x, gauges_value_y]
    else:
        gauges_value = pd.DataFrame(columns=np.arange(num_gauges))
    for i in range(num_section):
        file_name = os.path.join(case_folder, str(i), 'output',
                                 file_tag+'_gauges.dat')
        if ind_list[i].size>0:
            t_value = np.loadtxt(file_name, dtype='float64', ndmin=2)
            values = t_value[:,1:]
            times = t_value[:, 0]
            if file_tag == 'hU':
                gauges_value_x['times'] = times
                gauges_value_y['times'] = times
                values_x = values[:, 0::2]
                values_y = values[:, 1::2]
                gauges_value_x.iloc[0:values.shape[0], ind_list[i]] = values_x
                gauges_value_y.iloc[0:values.shape[0], ind_list[i]] = values_y
            else:
                gauges_value['times'] = times
                gauges_value.iloc[0:values.shape[0], ind_list[i]] = values        
    return gauges_value, gauges_pos_all
    
def _combine_multi_gpu_gauges_data(header_list, case_folder, file_tag):
    """ Combine gauges outputs from multi-gpu models according to gauges

    position data.

    gauges_pos.dat for each domain must be available

    """
    gauges_array = []
    value_array = []
    for i in range(len(header_list)):
        domain_header = header_list[i]
        gauge_pos_file = os.path.join(case_folder, str(i), 'input', 'field',
                                      'gauges_pos.dat')
        gauge_xy = np.loadtxt(gauge_pos_file, dtype='float64', ndmin=2)
        if gauge_xy.size>=2:
            gauge_ind = _find_gauges_inside_domain(domain_header, gauge_xy)
            gauges_array.append(gauge_xy[gauge_ind,:])
            file_name = os.path.join(case_folder, str(i), 'output',
                                     file_tag+'_gauges.dat')
            times, values = _read_one_gauge_file(file_name, gauge_ind)
            value_array.append(values)
    gauges_array = np.concatenate(gauges_array, axis=0)
    gauges_array, ind = np.unique(gauges_array, axis=0, return_index=True)
    if values.ndim == 2:
        value_array = np.concatenate(value_array, axis=1)
        value_array = value_array[:, ind]
    else: # values.ndim == 3
        value_array = np.concatenate(value_array, axis=2)
        value_array = value_array[:, :, ind]
    return gauges_array, times, value_array

def _find_gauges_inside_domain(domain_header, gauge_xy):
    """ Find the gauges inside a domain

    domain_header: (dict) header of the domain DEM

    gauge_xy: xy coordinate of the gauges

    Return: (numpy array) index of gauges inside the model domain

    """
    extent = sp.header2extent(domain_header)
    ind_x = np.logical_and(gauge_xy[:, 0] > extent[0],
                           gauge_xy[:, 0] < extent[1])
    ind_y = np.logical_and(gauge_xy[:, 1] > extent[2],
                           gauge_xy[:, 1] < extent[3])
    gauge_ind = np.where(np.logical_and(ind_x, ind_y))
    gauge_ind = gauge_ind[0]
    return gauge_ind
        
def _read_one_gauge_file(file_name, gauge_ind=None):
    """ Read a gauge file and return time series and values with outside gauges
    removed

    Supporting function to read_gauges_file

    """
    t_value = np.loadtxt(file_name, dtype='float64', ndmin=2)
    times = t_value[:, 0]
    values = t_value[:, 1:]
    if 'hU_gauges.dat' in file_name:
        values = np.array([values[:, 0::2], values[:, 1::2]])
    if gauge_ind is not None:
        if values.ndim==2:
            values = values[:, gauge_ind]
        else: #ndim=3
            values = values[:, :, gauge_ind]
    return times, values

def _header2row_numbers(local_header, global_header):
    """Calculate local grid starting and ending rows in global grid

    Return:

        ind_top: the index of the top row

        ind_bottom: the index of the bottom row
        
    """
    y_bottom_gap = local_header['yllcorner']-global_header['yllcorner']
    row_gap = round(y_bottom_gap/local_header['cellsize'])
    ind_bottom = global_header['nrows']-1-row_gap
    ind_top = ind_bottom-local_header['nrows']+1
    ind_top = int(ind_top)
    ind_bottom = int(ind_bottom)
    return ind_top, ind_bottom
    
def _header_local2global(header_list):
    """Convert local headers to a global header
    """
    extent_array = []
    for header in header_list:
        extent = sp.header2extent(header)
        extent_array.append(extent)
    extent_array = np.asarray(extent_array)
    header_global = header_list[0].copy()
    y_bottom = extent_array[:,2].min()
    header_global['yllcorner'] = y_bottom
    y_top = extent_array[:,3].max()
    nrows = (y_top-y_bottom)/header_global['cellsize']
    header_global['nrows'] = int(round(nrows))
    return header_global
