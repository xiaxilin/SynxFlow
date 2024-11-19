// ======================================================================================
// Name                :    High-Performance Integrated Modelling System
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software. Legacy name: GeoClasses
// ======================================================================================
// Version             :    1.0.1 
// Author              :    Xilin Xia
// Create Time         :    2014/10/04
// Update Time         :    2020/04/26
// ======================================================================================
// LICENCE: GPLv3 
// ======================================================================================

/*!
  \file field_reader.cc
  \brief Source file for basic field file reader class

*/

#include "field_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;


namespace GC{

  fieldReader::fieldReader(const char* filename){
    readin_field(filename);
  }

  void fieldReader::readin_field_netcdf(const char* filename) {
      try {
          // Open the NetCDF file
          NcFile ncFile(filename, NcFile::read);

          // Reading number of elements
          size_t num_elements = ncFile.getDim("num_cells").getSize();

          // Read cell IDs and cell values
          NcVar cell_ids_var = ncFile.getVar("cell_ids");
          NcVar cell_values_var = ncFile.getVar("cell_values");

          if (cell_ids_var.isNull() || cell_values_var.isNull()) {
              std::cout << "Error: Required variables not found in the NetCDF file." << std::endl;
              return;
          }

          // Read data from variables
          std::vector<int> cell_ids(num_elements);
          cell_ids_var.getVar(cell_ids.data());

          // Determine the number of columns in the cell values
          NcDim xy_dim = ncFile.getDim("xy_dim");
          size_t num_columns = xy_dim.isNull() ? 1 : xy_dim.getSize();

          // Read cell values data
          std::vector<double> cell_values(num_elements * num_columns);
          cell_values_var.getVar(cell_values.data());

          // Convert and insert data into the map
          for (size_t i = 0; i < num_elements; i++) {
              Vector3 value;
              if (num_columns == 1) {
                  value = Vector3(cell_values[i], 0.0, 0.0);  // Single value, set (x, 0, 0)
              } else if (num_columns == 2) {
                  value = Vector3(cell_values[i * 2], cell_values[i * 2 + 1], 0.0);  // Two values, set (x, y, 0)
              } else {
                  std::cout << "Unexpected number of columns in cell values: " << num_columns << std::endl;
                  return;
              }

              data.push_back({ cell_ids[i], value });
          }

          // Check and read boundary data if it exists
          if (!ncFile.getDim("num_boundaries").isNull()) {
              size_t num_boundaries = ncFile.getDim("num_boundaries").getSize();
              NcVar boundary_ids_var = ncFile.getVar("boundary_ids");
              NcVar boundary_values_var = ncFile.getVar("boundary_values");

              if (boundary_ids_var.isNull() || boundary_values_var.isNull()) {
                  std::cout << "Error: Required boundary variables not found in the NetCDF file." << std::endl;
                  return;
              }

              // Read boundary data into vectors
              std::vector<int> boundary_ids(num_boundaries);
              std::vector<int> boundary_values(num_boundaries * 3);  // Assuming 3 columns for boundary data

              boundary_ids_var.getVar(boundary_ids.data());
              boundary_values_var.getVar(boundary_values.data());

              // Insert boundary data into the map
              for (size_t i = 0; i < num_boundaries; i++) {
                  int primary_type = boundary_values[i * 3];
                  int secondary_type = boundary_values[i * 3 + 1];
                  int source_id = boundary_values[i * 3 + 2];
                  boundary_type.push_back({ boundary_ids[i], ShortTripleFlag(primary_type, secondary_type, source_id) });
              }

          }

          std::cout << "Data successfully read from file:" << filename << std::endl;

      } catch (NcException& e) {
          std::cerr << "NetCDF error: " << e.what() << std::endl;
      } catch (std::exception& e) {
          std::cerr << "Standard exception: " << e.what() << std::endl;
      }
  }

  void fieldReader::readin_field(const char* filename){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      std::cout << "error: unable to open input file: " << filename << std::endl;
      return;
    }
    std::string line;
    std::string word;
    //reading elements
    getline(input, line);
    size_t n_elements;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> n_elements;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Vector3 value;
    Flag id_element;
    for (size_t i = 0; i < n_elements; i++){
      getline(input, line);
      stream.str(line);
      value = 0.0;
      stream >> id_element >> value;
      data.push_back({ id_element, value });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
    getline(input, line);
    getline(input, line);
    size_t n_boundaries;
    stream.str(line);
    stream >> n_boundaries;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Flag id_boundary, primary_type, secondary_type, source_id;
    for (size_t i = 0; i < n_boundaries; i++){
      getline(input, line);
      stream.str(line);
      primary_type = 0;
      secondary_type = 0;
      source_id = 0;
      stream >> id_boundary >> primary_type >> secondary_type >>source_id;
      boundary_type.push_back({ id_boundary, ShortTripleFlag(primary_type, secondary_type, source_id) });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
  }

  bool completeFieldReader::readin_region_mask_netcdf(const char* filename) {

      try {
          // Open the NetCDF file
          NcFile ncFile(filename, NcFile::read);

          // Read the number of elements from the dimension
          NcDim num_elements_dim = ncFile.getDim("num_cells");
          if (num_elements_dim.isNull()) {
              std::cerr << "Error: Dimension 'num_cells' not found in NetCDF file." << std::endl;
              return false;
          }

          size_t n_elements = num_elements_dim.getSize();

          // Read the IDs and values
          NcVar element_ids_var = ncFile.getVar("cell_ids");
          NcVar element_values_var = ncFile.getVar("cell_values");

          if (element_ids_var.isNull() || element_values_var.isNull()) {
              std::cerr << "Error: Variables 'cell_ids' or 'cell_values' not found in NetCDF file." << std::endl;
              return false;
          }

          // Read data into vectors
          std::vector<Flag> element_ids(n_elements);
          std::vector<Flag> element_values(n_elements);

          element_ids_var.getVar(element_ids.data());
          element_values_var.getVar(element_values.data());

          // Insert data into the region_mask map
          for (size_t i = 0; i < n_elements; i++) {
              region_mask.push_back({ element_ids[i], element_values[i] });
          }

          return true;

      } catch (NcException& e) {
          //std::cerr << "NetCDF error: " << e.what() << std::endl;
          return false;
      } catch (std::exception& e) {
          //std::cerr << "Standard exception: " << e.what() << std::endl;
          return false;
      }
  }

  bool completeFieldReader::readin_region_mask(const char* filename){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    std::string line;
    std::string word;
    //reading elements
    getline(input, line);
    size_t n_elements;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> n_elements;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Flag value;
    Flag id_element;
    for (size_t i = 0; i < n_elements; i++){
      getline(input, line);
      stream.str(line);
      value = 0;
      stream >> id_element >> value;
      region_mask.push_back({ id_element, value });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
    return true;
  }

  completeFieldReader::completeFieldReader(const char* path, const char* field_name){
    std::string filename = std::string(path) + std::string(field_name) + ".nc";
    readin_field_netcdf(filename.c_str());
    for (unsigned int i = 0;; ++i){
      std::string filename;
      std::ostringstream file_id;
      file_id << i;
      filename = std::string(path) + std::string(field_name) + "_BC_" + file_id.str() + ".dat";
      bool success = readin_boundary_source(filename.c_str(), i);
      if (!success){
        break;
      }
    }
    filename = std::string(path) + std::string(field_name) + "_mask.nc";
    if (readin_region_mask_netcdf(filename.c_str())){
      std::string filename = std::string(path) + std::string(field_name) + "_source_all.dat";
      std::ifstream input(filename);
      if(input){
        int num_of_regions;
        std::istringstream stream;
        std::string line;
        Scalar t;
        Scalar value;
        getline(input, line);
        stream.str(line);
        stream >> num_of_regions;
        stream.str(std::string());
        stream.clear();		//clear the istringstream
        data_time_series.resize(num_of_regions);
        data_source.resize(num_of_regions);
        while (!input.eof()){
          getline(input, line);
          stream.str(line);
          stream >> t;
          for (int i = 0; i < num_of_regions; i++){
            stream >> value;
            data_time_series[i].push_back(t);
            data_source[i].push_back(value);
          }
          stream.str(std::string());
          stream.clear();		//clear the istringstream
        }
      }
      else{
        for (unsigned int i = 0;; ++i){
          std::string filename;
          std::ostringstream file_id;
          file_id << i;
          filename = std::string(path) + std::string(field_name) + "_source_" + file_id.str() + ".dat";
          bool success = readin_data_source(filename.c_str(), i);
          if (!success){
            break;
          }
        }
      }
    }
  }

  bool completeFieldReader::readin_boundary_source(const char* filename, const unsigned int cnt){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    else{
      time_series.resize(cnt + 1);
      boundary_source.resize(cnt + 1);
      std::istringstream stream;
      std::string line;
      Scalar t;
      Vector3 value;
      while (!input.eof()){
        getline(input, line);
        stream.str(line);
        stream >> t >> value;
        time_series[cnt].push_back(t);
        boundary_source[cnt].push_back(value);
        stream.str(std::string());
        stream.clear();		//clear the istringstream
      }
      return true;
    }
  }

  bool completeFieldReader::readin_data_source(const char* filename, const unsigned int cnt){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    else{
      data_time_series.resize(cnt + 1);
      data_source.resize(cnt + 1);
      std::istringstream stream;
      std::string line;
      Scalar t;
      Vector3 value;
      while (!input.eof()){
        getline(input, line);
        stream.str(line);
        stream >> t >> value;
        data_time_series[cnt].push_back(t);
        data_source[cnt].push_back(value);
        stream.str(std::string());
        stream.clear();		//clear the istringstream
      }
      return true;
    }
  }

}
