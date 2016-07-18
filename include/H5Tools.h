
/*******************************************************************************************//** @file
 *  		
 * 	file: 		H5Tools.h
 * 	contents:  	Helper Functions to write containers to hdf5 file
 * 
 ****************************************************************************************************/


#pragma once

#include <string>
#include <grid.h>
#include <vector>
#include <boost/multi_array.hpp>
#include <H5Cpp.h>

typedef std::complex<double> dcomplex; 

// -- Convenient write functions to HDF5 group
void write( const F_Grid& fgrid, H5::Group& group, const std::string& dataset_name = std::string("fgrid") ); 	///<   	Write fermionic grid to given group of hdf5 file
void write( const Bos_Grid& bgrid, H5::Group& group, const std::string& dataset_name = std::string("bgrid") ); 	///<   	Write bosonic grid to given group of hdf5 file
void write( double scalar, H5::Group& group, const std::string& dataset_name ); 				///<   	Write a scalar value to given group of hdf5 file
void write( int integer, H5::Group& group, const std::string& dataset_name ); 					///<   	Write a integer value to given group of hdf5 file
void write( std::string text, H5::Group& group, const std::string& dataset_name ); 				///<   	Write a scalar value to given group of hdf5 file
void write( std::vector<double> vec_scalar, H5::Group& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file

// -- Convenient write functions for Boost Multi-Array and Grids
template< size_t ndims > void write( const boost::multi_array<double, ndims>& my_arr, H5::Group& group, const std::string& dataset_name );  ///< Writes double valued boost multi_array to given group of hdf5 file
template< size_t ndims > void write( const boost::multi_array<std::complex<double>, ndims>& my_arr, H5::Group& group, const std::string& dataset_name = std::string("") );  ///< Writes complex<double> valued boost multi_array to given group of hdf5 file
template< size_t ndims > void read( boost::multi_array<double, ndims>& my_arr, H5::Group& group, const std::string& dataset_name );  ///< Reads double valued boost multi_array to given group from hdf5 file
template< size_t ndims > void read( boost::multi_array<std::complex<double>, ndims>& my_arr, H5::Group& group, const std::string& dataset_name = std::string("") );  ///< Reads complex<double> valued boost multi_array from given group of hdf5 file

#include <H5Tools_impl.h> 	// contains implementations of template functions
