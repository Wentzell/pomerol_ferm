
/******************************************************************************************//** @file
 *  		
 * 	file: 		grid.h
 * 	contents:	Classes for frequency and momentum grid  
 * 
 ****************************************************************************************************/


#pragma once

#include <vector>
#include <string>
#include <iostream>

const double PI = 3.141592653589793238462; 

class F_Grid
{	
   public: 	
      /********************* Constructor ********************/
      F_Grid( unsigned int pos_freq_count, double step_size );	///< Initialing constructor for equidistant frequency grid
      //F_Grid( F_Grid&& ) = default;			///< Default move constructor

      /***************************************************/
      void print( std::string fname = "f_grid.dat");	///< Prints frequency grid to file "f_grid.dat"

      inline const double& operator[]( int w )
      {
	 return grid_points[w + pos_freq_count]; 
      }
      
      inline const double * data() const
      {
	 return grid_points.data(); 
      }

      unsigned int get_pos_freq_count() const;
      double get_step_size() const;

   private:

      std::vector<double> grid_points;  		///< Values of grid
      const int 	 pos_freq_count;		///< Number of elements in grid
      const double 	 step_size;			///< Step size of equidistant grid 

};

class Bos_Grid
{	
   public:
      /********************* Constructor ********************/
      Bos_Grid( unsigned int pos_freq_count, double step_size );	///< Initialing constructor for equidistant frequency grid
      //Bos_Grid( Bos_Grid&& ) = default;			///< Default move constructor

      /***************************************************/
      void print( std::string fname = "bos_grid.dat");	///< Prints frequency grid to file "f_grid.dat"

      inline const double& operator[]( int w )
      {
	 return grid_points[w + pos_freq_count]; 
      }

      inline const double * data() const
      {
	 return grid_points.data(); 
      }

      unsigned int get_pos_freq_count() const;
      double get_step_size() const;

   private:

      std::vector<double> grid_points; 			///< Values of grid
      const int 	 pos_freq_count;		///< Number of elements in grid
      const double 	 step_size;			///< Step size of equidistant grid 

};

class K_Grid : public std::vector<std::pair<double, double>>
{	
   public:

      /********************* Constructor ********************/
      K_Grid();						///< Constructor for 8-patch momentum grid

      /***************************************************/
      void print( std::string fname = "k_grid.dat");	///< Prints frequency grid to file "f_grid.dat"
      unsigned int get_patch_count();

   private:

      int	patch_count;				///< Number of elements in grid
};

const int mirror_mom_vert_arr[8] = { 0, 4, 3, 2, 1, 5, 6, 7 }; ///< Array that specifies how to mirror single momentum index at vertical axis
const int mirror_mom_diag_arr[8] = { 0, 1, 4, 3, 2, 7, 6, 5 }; ///< Array that specifies how to mirror single momentum index at diagonal ( bottom left to top right ) axis
const int mirror_mom_pipi_arr[8] = { 6, 1, 2, 3, 4, 7, 0, 5 }; ///< Array that specifies how to calculate ( pi, pi ) - k for single momentum k
const int neg_k_arr[8] = {0, 3, 4, 1, 2, 5, 6, 7}; ///< Array that specifies how to change sign of single momentum index ( corresponds to 2 rotations )
const int rot_k_arr[8] = {0, 2, 3, 4, 1, 7, 6, 5}; ///< Array that specifies how to rotate single momentum index by 90 degrees clockwise

/**
 *	Array that specifies how to add two momenta in the language of the patch indeces
 */
const int add_k_arr[8][8] =   { { 0, 1, 2, 3, 4, 5, 6, 7 }, 
	   			{ 1, 6, 7, 0, 5, 2, 3, 4 }, 
				{ 2, 7, 6, 5, 0, 1, 4, 3 }, 
				{ 3, 0, 5, 6, 7, 4, 1, 2 }, 
   				{ 4, 5, 0, 7, 6, 3, 2, 1 }, 
				{ 5, 2, 1, 4, 3, 0, 7, 6 }, 
				{ 6, 3, 4, 1, 2, 7, 0, 5 }, 
				{ 7, 4, 3, 2, 1, 6, 5, 0 } } ;

/**
 *	Array that specifies how to subtract two momenta in the language of the patch indeces
 */
const int dif_k_arr[8][8] =  { { 0, 3, 4, 1, 2, 5, 6, 7 },
  				{ 1, 0, 5, 6, 7, 2, 3, 4 }, 
				{ 2, 5, 0, 7, 6, 1, 4, 3 }, 
  				{ 3, 6, 7, 0, 5, 4, 1, 2 }, 
  				{ 4, 7, 6, 5, 0, 3, 2, 1 }, 
  				{ 5, 4, 3, 2, 1, 0, 7, 6 }, 
  				{ 6, 1, 2, 3, 4, 7, 0, 5 }, 
  				{ 7, 2, 1, 4, 3, 6, 5, 0 } };  


inline int add_k( int k1, int k2 )
{
   return add_k_arr[k1][k2]; 
}

inline int dif_k( int k1, int k2 )
{
   return dif_k_arr[k1][k2]; 
}

inline int neg_k( int k )
{
   return neg_k_arr[k]; 
}

inline int rot_k( int k )
{
   return rot_k_arr[k]; 
}
