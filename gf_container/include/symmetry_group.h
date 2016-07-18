
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmetries.h
 * 	contents:  	Defines a class that allows to initialize elements of the gf container
 * 			using symmetry relations in an efficient way
 * 
 ****************************************************************************************************/

#pragma once

#include <gf.h>
#include <complex>
#include <iostream>

/**
 *	Set of possible operations after symmetry operation
 */
class operation : public std::pair<bool, bool>
{
   public:
      using base_t = std::pair<bool, bool>; 
      using base_t::first; 
      using base_t::second; 

      operation( bool first_, bool second_ ) ///< Constructor taking two bools as argument
      {
	 first = first_ ;	// Apply minus sign?
	 second = second_;	// Apply complex conjugation?
      }

      operation  operator*( const operation& b )	///< Overload multiplication operator for successive application of operations
      {
	 return operation(  (*this ).first xor b.first, (*this ).second xor b.second ) ;
      }

      std::complex<double> operator()( const std::complex<double> val )		///< Apply operation
      {
	 if( first )
	 {
	    if( second )
	       return -conj( val ); 
	    return -val; 
	 }

	 if( second )
	    return conj( val ); 

	 return val; 
      }

   private:
};


/**
 *	Elements of vertex tensor. States index of independent coupling array and possible operations on it
 */
struct symm_idx_t
{
   public:
      unsigned long int idx;	///< Index in the raw data array of container
      operation oper; 	///< Possible operations that relate two tensor elements. ( sign change, complex conjugation )

      symm_idx_t():
	 idx( 0 ), oper( false, false )
   {}

      symm_idx_t( int idx_, bool first_ = false, bool second_ = false ):
	 idx( idx_ ), oper( first_, second_ )
   {}

      symm_idx_t( int idx_, operation oper_ ):
	 idx( idx_ ), oper( oper_ )
   {}
};

template< unsigned int ndims_ >
class symmetry_grp_t
{
   public:
      using dcomplex = std::complex<double>; 

      using idx_t = idx_obj_t<ndims_>; 
      using type = symmetry_grp_t<ndims_>; 
      using gf_t = gf<dcomplex, ndims_>; 
      using symm_func_t = boost::function<operation ( idx_t& idx )>; 
      using init_func_t = boost::function<dcomplex ( const idx_t& idx )>; 

      static constexpr unsigned int ndims = ndims_;                     ///< The number of dimensions        

      void init( gf<dcomplex,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp parallel for schedule( dynamic )
	 for( int i = 0; i < symm_classes.size(); ++i )
	 {
	    dcomplex val = init_func( gf_obj.get_idx( symm_classes[i][0].idx ) ); 
	    for( auto symm_idx : symm_classes[i] )
	       gf_obj( symm_idx.idx ) = symm_idx.oper( val ); 
	 }

      }

      symmetry_grp_t( const type& symm_grp ):
	 symm_lst( symm_grp.symm_lst ), symm_classes( symm_grp.symm_classes )
   {}

      symmetry_grp_t( type&& symm_grp ):
	 symm_lst( symm_grp.symm_lst ), symm_classes( std::move( symm_grp.symm_classes ) )
   {}

      type& operator=( const type& symm_grp )
      {
	 symm_lst = symm_grp.symm_lst; 
	 symm_classes = symm_grp.symm_classes; 
	 return *this; 
      }

      type& operator=( type&& symm_grp )
      {
	 symm_lst = symm_grp.symm_lst; 
	 symm_classes.operator=( std::move( symm_grp.symm_classes ) ); 
	 return *this; 
      }

      symmetry_grp_t( const gf_t& gf_obj, const std::vector< symm_func_t >& symm_lst_ ):
	 symm_lst(symm_lst_),
	 shape_arr( reinterpret_cast<const boost::array< boost::multi_array_types::size_type, ndims >& >( *(gf_obj.shape_arr) ) ), 
	 idx_bases( reinterpret_cast<const boost::array< boost::multi_array_types::index, ndims >& >( *(gf_obj.idx_bases) ) )
   {
      std::cout << " Initializing symmetry group for container of Length : " << ndims << std::endl; 
      std::cout << "   Shape: "; 
      for( int i = 0; i < ndims; ++i )
	 std::cout << shape_arr[i] << " "; 
      std::cout << "   Idx_bases: "; 
      for( int i = 0; i < ndims; ++i )
	 std::cout << idx_bases[i] << " ";

      // For frequencies, consider going out of range with symmetry functions
      for( int i = 0; i < ndims; ++i )
	 if( idx_bases[i] != 0 ){ shape_arr[i] *= 3; } 

      // Create bool array to track checked elements
      gf< bool, ndims > checked( shape_arr ); 
      checked.reindex( idx_bases ); 
      checked.init( []( const idx_t& idx ){ return false; } ); 

      for( long unsigned int iter = 0; iter < gf_obj.num_elements(); ++iter )
      {
	 idx_t idx( gf_obj.get_idx( iter ) ); 

	 if ( !(checked( idx )) )  								// if tensor object not yet related to any other
	 {
	    checked( idx ) = true;
	    std::vector< symm_idx_t > current_class { symm_idx_t(iter) }; 	// initialize new symmetry class
	    operation track_op( false, false );
	    iterate( idx, gf_obj, track_op, checked, current_class ); 		// start iterating on index 
	    symm_classes.push_back( current_class );				// Add current symmetry class to list
	 }
      }

      std::cout << std::endl << " Symmetries reduction " << 1.0 * symm_classes.size() / gf_obj.num_elements() << std::endl; 
   }

   private:
      std::vector< symm_func_t > symm_lst; 
      std::vector< std::vector< symm_idx_t > > symm_classes; 
      boost::array<size_t, ndims> shape_arr; 
      boost::array<long, ndims> idx_bases; 

      void iterate( const idx_t& idx_it, const gf<dcomplex,ndims>& gf_obj,  const operation& track_op, gf<bool,ndims>& checked, std::vector< symm_idx_t >& current_class  )
      {
	 for( auto symm: symm_lst ) 				// iterate over list of all symmetries specified
	 {
	    idx_t idx = idx_it;					// copy idx
	    operation curr_op = symm( idx ) * track_op;		// apply symmetry operation and track operations applied

	    if( !checked( idx ) )				// if index not yet checked
	    { 
	       checked( idx ) = true;

	       if( gf_obj.is_valid( idx ) ) 
		  current_class.push_back( symm_idx_t( gf_obj.get_pos_1d( idx ), curr_op ) ); // if valid index, add to current symmetry class

	       iterate( idx, gf_obj, curr_op, checked, current_class );	// iterate further	
	    }
	 }
      }
}; 
