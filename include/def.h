
/******************************************************************************************//** @file
 *  		
 * 	file: 		def.h
 * 	contents:  	Definition of used correlation function containers ( wrapper around gf container )
 * 
 ****************************************************************************************************/

#pragma once

#include <complex>
#include <boost/multi_array.hpp>
#include <gf.h>

const int PATCH_COUNT = 1; 

using dcomplex = std::complex<double>;						///< Complex double type

#define INSERT_COPY_AND_ASSIGN(X) 					\
X( const X & gf_obj ):    						\
   base_t( gf_obj )							\
{}       								\
X( X && gf_obj ):							\
   base_t( std::move(gf_obj) )						\
{}      								\
X & operator=( const X & gf_obj )					\
{									\
   base_t::operator=( gf_obj ); 					\
   return *this; 							\
} 									\
X & operator=( X && gf_obj )						\
{									\
   base_t::operator=( std::move( gf_obj) ); 				\
   return *this; 							\
} 

enum class I1P{ w, k, s_in, s_out }; 
class gf_1p_t : public  gf< dcomplex, 4 > 		///< Container type for one-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_1p_t( int pos_freq_count_, int qn_count_ ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][PATCH_COUNT][qn_count_][qn_count_] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_t)
}; 
using idx_1p_t = gf_1p_t::idx_t; 

enum class I2P{ w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_2p_t : public gf< dcomplex, 10 > 		///< Container type for two-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 10 >; 

      gf_2p_t( int pos_freq_count_, int qn_count_ ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][ffreq(pos_freq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][PATCH_COUNT][PATCH_COUNT]
	       [qn_count_][qn_count_][qn_count_][qn_count_] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_2p_t)
}; 
using idx_2p_t = gf_2p_t::idx_t; 
