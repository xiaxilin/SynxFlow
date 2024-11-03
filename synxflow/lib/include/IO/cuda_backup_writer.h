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
\file cuda_backup_writer.h
\brief Header file for cuda simple field writer class

*/

#ifndef CUDA_BACKUP_WRITER_H
#define CUDA_BACKUP_WRITER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include "Flag.h"
#include "Vector.h"
#include "Scalar.h"
#include "cuda_mapped_field.h"
#include <netcdf>
#include <map>

using namespace netCDF;
using namespace netCDF::exceptions;

namespace GC{


/*
  template <class T>
  void cuBackupWriter(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t, const char* directory = "output/"){
    std::ofstream fout;
    std::ostringstream file_id, out_time;
    out_time << t;
    file_id.fill('0');
    file_id.width(3);
    std::string name = std::string(directory) + std::string(filename) + "_" + out_time.str() + ".dat";
    fout.open(name.c_str());
    if (!fout){
      std::cout << "Unable to create output file." << std::endl;
    }
    phi.data.sync();
    phi.boundary.sync();
    phi.mesh->boundary2opposite_handles.sync();
    fout << "$Element Number" << std::endl;
    fout << phi.data.size() << std::endl;
    fout << "$Element_id  Value" << std::endl;
    fout.flags(std::ios::scientific);
    auto size = phi.data.size();
    auto data_array = phi.data.host_ptr();
    auto boundary_size = phi.boundary.size();
    auto boundary_array = phi.boundary.host_ptr();
    auto boundary_cell_handles = phi.mesh->boundary2opposite_handles.host_ptr();
    for (unsigned int i = 0; i < size; ++i){
      fout << i << "  " << std::setprecision(20) << data_array[i] << std::endl;
    }
    std::map<Flag, ShortTripleFlag> boundary_type;
    for (unsigned int i = 0; i < boundary_size; ++i){
      unsigned int cell_id = boundary_cell_handles[i].get_global_id();
      boundary_type.insert({ cell_id, boundary_array[i] });
    }
    fout << "$Boundary Numbers" << std::endl;
    fout << boundary_type.size() << std::endl;
    fout << "$Element_id Boundary_type" << std::endl;
    for (auto cellid_boundary : boundary_type){
      fout << cellid_boundary.first << "  " << cellid_boundary.second << std::endl;
    }
    fout.close();
  } */

  template <class T>
  void cuBackupWriter(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t, const char* directory = "output/") {
      try {
          std::ostringstream file_id, out_time;
          out_time << t;
          file_id.fill('0');
          file_id.width(3);
          std::string name = std::string(directory) + std::string(filename) + "_" + out_time.str() + ".nc";

          // Create a new NetCDF file
          NcFile dataFile(name, NcFile::replace);

          // Sync data to ensure they are up-to-date
          phi.data.sync();
          phi.boundary.sync();
          phi.mesh->boundary2opposite_handles.sync();

          auto size = phi.data.size();
          auto data_array = phi.data.host_ptr();
          auto boundary_size = phi.boundary.size();
          auto boundary_array = phi.boundary.host_ptr();
          auto boundary_cell_handles = phi.mesh->boundary2opposite_handles.host_ptr();

          // Define dimensions
          NcDim numElementsDim = dataFile.addDim("num_elements", size);
          NcDim numBoundariesDim = dataFile.addDim("num_boundaries", boundary_size);
          NcDim tripletDim = dataFile.addDim("triplet_components", 3); // For three int values


          // Define variables based on data type
          NcVar elementIdsVar = dataFile.addVar("cell_ids", ncUint, numElementsDim);
          if constexpr (std::is_same<T, Scalar>::value) {
              // Scalar data
              NcVar elementValuesVar = dataFile.addVar("cell_values", ncDouble, numElementsDim);
              elementValuesVar.putVar(data_array);
          } else if constexpr (std::is_same<T, Vector2>::value) {
              // Vector2 data
              NcDim componentDim = dataFile.addDim("components", 2);
              NcVar elementValuesVar = dataFile.addVar("cell_values", ncDouble, {numElementsDim, componentDim});

              // Create a 2D array to store the (x, y) components
              std::vector<std::vector<double>> vectorData(size, std::vector<double>(2));
              for (unsigned int i = 0; i < size; ++i) {
                  vectorData[i][0] = data_array[i].x; // x component
                  vectorData[i][1] = data_array[i].y; // y component
              }

              // Flatten the 2D array for writing to NetCDF
              std::vector<double> flattenedData;
              for (const auto& vec : vectorData) {
                  flattenedData.insert(flattenedData.end(), vec.begin(), vec.end());
              }
              elementValuesVar.putVar(flattenedData.data());
          }

          // Prepare and write boundary data
          std::map<Flag, ShortTripleFlag> boundary_type;
          for (unsigned int i = 0; i < boundary_size; ++i) {
              unsigned int cell_id = boundary_cell_handles[i].get_global_id();
              boundary_type.insert({cell_id, boundary_array[i]});
          }

          // Create arrays for boundary data
          std::vector<unsigned int> boundaryElementIds;
          std::vector<std::vector<int>> boundaryTypes(boundary_type.size(), std::vector<int>(3));

          int index = 0;
          for (const auto& cellid_boundary : boundary_type) {
              boundaryElementIds.push_back(cellid_boundary.first);
              boundaryTypes[index][0] = cellid_boundary.second.getx(); // First int value
              boundaryTypes[index][1] = cellid_boundary.second.gety(); // Second int value
              boundaryTypes[index][2] = cellid_boundary.second.getz(); // Third int value
              ++index;
          }

          // Flatten the boundary types array for NetCDF storage
          std::vector<int> flattenedBoundaryTypes;
          for (const auto& triplet : boundaryTypes) {
              flattenedBoundaryTypes.insert(flattenedBoundaryTypes.end(), triplet.begin(), triplet.end());
          }

          // Write boundary data to NetCDF file
          NcVar boundaryElementIdsVar = dataFile.addVar("boundary_ids", ncUint, numBoundariesDim);
          NcVar boundaryTypesVar = dataFile.addVar("boundary_values", ncInt, {numBoundariesDim, tripletDim});

          boundaryElementIdsVar.putVar(boundaryElementIds.data());
          boundaryTypesVar.putVar(flattenedBoundaryTypes.data());

          std::cout << "Data successfully written to NetCDF file: " << name << std::endl;
      } catch (NcException& e) {
          std::cerr << "NetCDF error: " << e.what() << std::endl;
      } catch (std::exception& e) {
          std::cerr << "Standard exception: " << e.what() << std::endl;
      }
  }

}

#endif
