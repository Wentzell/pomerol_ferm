
/****************************************************************************************************
 *  		
 * 	file: 		grid.cpp
 * 	contents:   	See grid.h
 * 
 ****************************************************************************************************/


#include <grid.h>
#include <fstream>
#include <string>

using namespace std; 

/********************* F_Grid Class  ********************/

F_Grid::F_Grid( unsigned int pos_freq_count_, double step_size_ ) :
   pos_freq_count( pos_freq_count_ ), step_size( step_size_ )
{
   for( int i = -pos_freq_count; i < pos_freq_count ; ++i )
      grid_points.push_back( step_size * i + step_size / 2.0 );
}

void F_Grid::print( string fname )
{
   fstream data;
   data.open("dat/" + fname, ios::out | ios::trunc );
   for( double freq : grid_points )
      data << freq << endl;
   data.flush();
   data.close();
}

unsigned int F_Grid::get_pos_freq_count() const
{
   return pos_freq_count;
}

double F_Grid::get_step_size() const
{
   return step_size;
}

/********************* Bos_Grid Class  ********************/

Bos_Grid::Bos_Grid( unsigned int pos_freq_count_, double step_size_ ) :
   pos_freq_count( pos_freq_count_ ), step_size( step_size_ )
{
   for( int i = -pos_freq_count; i < pos_freq_count + 1 ; ++i )
      grid_points.push_back( step_size * i );
}

void Bos_Grid::print( string fname )
{
   fstream data;
   data.open("dat/" + fname, ios::out | ios::trunc );
   for( double freq : grid_points )
      data << freq << endl;
   data.flush();
   data.close();
}

unsigned int Bos_Grid::get_pos_freq_count() const
{
   return pos_freq_count;
}

double Bos_Grid::get_step_size() const
{
   return step_size;
}

/********************* K_Grid Class ********************/

K_Grid::K_Grid()
{
   push_back( pair<double, double>( 0.0, 	0.0 ));		// PATCH 0
   push_back( pair<double, double>( PI/2.0, 	PI/2.0 ));	// PATCH 1
   push_back( pair<double, double>( PI/2.0, 	-PI/2.0 ));	// PATCH 2
   push_back( pair<double, double>(-PI/2.0, 	-PI/2.0 ));	// PATCH 3
   push_back( pair<double, double>(-PI/2.0, 	PI/2.0 ));	// PATCH 4
   push_back( pair<double, double>( 0.0, 	PI ));		// PATCH 5
   push_back( pair<double, double>( PI, 	PI ));		// PATCH 6
   push_back( pair<double, double>( PI, 	0.0 ));		// PATCH 7
   patch_count = size();
}

unsigned int K_Grid::get_patch_count()
{
   return patch_count;
}

void K_Grid::print( string fname )
{
   fstream data;
   data.open("dat/" + fname, ios::out | ios::trunc );
   data << "# kx\t\tky" << endl;         
   for( pair<double, double> mom : *this )
      data << mom.first << "\t\t" << mom.second << endl;
   data.flush();
   data.close();
}

