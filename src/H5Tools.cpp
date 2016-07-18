
/*******************************************************************************************//** @file
 *  		
 * 	file: 		H5Tools.cpp
 * 	contents:  	See H5Tools.h
 * 
 ****************************************************************************************************/


#include <H5Tools.h>

using namespace H5; 
using std::string; 
using std::vector; 


void write( const F_Grid& fgrid, Group& group, const string& dataset_name )
{
   
   hsize_t grid_dim[1] = { 2* fgrid.get_pos_freq_count() };
   DataSpace grid_dataspace( 1, grid_dim );

   DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );

   grid_dataset.write( fgrid.data(), PredType::NATIVE_DOUBLE );

}

void write( const Bos_Grid& bgrid, Group& group, const string& dataset_name )
{
   
   hsize_t grid_dim[1] = { 2* bgrid.get_pos_freq_count() + 1 };
   DataSpace grid_dataspace( 1, grid_dim );

   DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );

   grid_dataset.write( bgrid.data(), PredType::NATIVE_DOUBLE );

}

void write( double scalar, Group& group, const string& attr_name )
{

   DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
   Attribute attr = group.createAttribute( attr_name, PredType::NATIVE_DOUBLE, attr_dsp );
   attr.write( PredType::NATIVE_DOUBLE, &scalar );

}

void write( int integer, Group& group, const string& attr_name )
{

   DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
   Attribute attr = group.createAttribute( attr_name, PredType::NATIVE_INT, attr_dsp );
   attr.write( PredType::NATIVE_DOUBLE, &integer );

}

void write( string text, H5::Group& group, const string& dataset_name )
{
   
   DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
   StrType strdatatype( PredType::C_S1, H5T_VARIABLE ); // String type with variable length
   Attribute attr = group.createAttribute( dataset_name, strdatatype, attr_dsp );
   attr.write( strdatatype, &text );

}

void write( vector<double> vec_scalar, H5::Group& group, const string& dataset_name )
{

   hsize_t dims[1] = { vec_scalar.size() }; 
   DataSpace dataspace( 1, dims );
   DataSet dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, dataspace );
   dataset.write( vec_scalar.data(), PredType::NATIVE_DOUBLE );
   
}
