#include <iostream>
#include <complex>
#include <boost/timer/timer.hpp>
#include <gf.h>
#include <symmetry_group.h>

using namespace std; 

const int N = 1; 

typedef complex<double> dcomplex; 

class gf_1p_t: public gf< dcomplex, 1 >
{
   public:
      using base_t = gf< dcomplex, 1 >; 

      gf_1p_t():                                                                                                                                                                                                                         
	 gf< dcomplex, 1 >( boost::extents[ffreq(N)] )
   {}
}; 

int main()
{
   cout << " Starting main " << endl; 

   // Enum that contains names of gf indeces
   enum class MIXED{ w, W };

   // Double-valued 2-dimensional gf
   using mygf_t = gf< dcomplex, 2 >;
   using idx_t = mygf_t::idx_t; 

   // Create object with one fermionic and one bosonic frequency
   mygf_t my_gf( boost::extents[ffreq(N)][bfreq(N)] ); 

   // Test initialization
   cout << endl << " --- Init Test --- " << endl << endl; 
   auto init_func = []( const idx_t& idx )->dcomplex{ cout << " call "; return idx(1); }; 

   cout << " Simple init: " << endl; 
   my_gf.init( init_func ); 
   cout << endl << " GF Values: " << my_gf << endl << endl; 

   my_gf.init( []( const idx_t& idx )->dcomplex{ return 0.0; } ); 

   cout << " Symmetry init: " << endl; 
   auto sym_func = []( idx_t& idx ){ idx(0) = -idx(0) - 1; return operation( false, false ); }; 
   symmetry_grp_t<2> sym_grp( my_gf, {sym_func} ); 
   sym_grp.init( my_gf, []( const idx_t& idx )->dcomplex{ cout << " call "; return idx(1); } ); 
   cout << endl << " GF Values: " << my_gf << endl << endl; 

   // Fill a vector with all possible indeces
   std::vector< idx_t > idx_lst; 
   my_gf.fill_idx_lst( idx_lst ); 

   // Access test
   cout << endl << " --- Acess Test --- " << endl << endl; 
   // Consider the element x,y
   int x = 0;
   int y = 0; 
   idx_t idx( { x, y } ); 
   cout << " idx(w) " << idx( MIXED::w ) << endl; 
   cout << " Direct access of gf[x][y] " << my_gf[x][y] << endl; 
   cout << " Acess with idx_t object " << my_gf( idx ) << endl; 
   cout << " Acess with corresponding pos1d " << my_gf( my_gf.get_pos_1d( idx ) ) << endl << endl; 

   // Output idx_t type object and check get_pos_1d and get_idx functions
   cout << " Output Test " << endl; 
   cout << " idx " << idx << endl; 
   cout << " get_idx( gf.pos_1d( idx ) ) " << my_gf.get_idx( my_gf.get_pos_1d( idx ) ) << endl << endl; 

   cout << " Test sum over all container elements " << endl; 
   dcomplex val ( 0.0, 0.0 ); 
   {  // .. with direct access
      boost::timer::auto_cpu_timer t;
      for( int w = -N; w < N; w++ )
	 for( int W = -N; W < N + 1; W++ )
	    val += my_gf[w][W]; 
   }
   cout << " .. using direct acess gf[][] : " << val << endl;  
   val = dcomplex( 0.0, 0.0 ); 
   {  // .. with idx_t access
      boost::timer::auto_cpu_timer t;
      for( auto idx : idx_lst )
	 val += my_gf( idx ); 
   }
   cout << " .. using idx acess gf() " << val << endl << endl; 

   // -- Access with idx_t of similar gf_t
   using mygf_other_t = gf< int, 2 >;
   mygf_other_t my_other_gf( boost::extents[ffreq(N+1)][bfreq(N+1)] ); 
   mygf_other_t::idx_t other_idx( { 0, 0} ); 
   idx_t cloned_idx( other_idx ); 
   my_gf ( other_idx ); 

   // --  Abs and Norm
   cout << " my_gf[-1][0] " << my_gf[-1][0] << endl; 
   cout << " abs(my_gf)[-1][0] " << abs(my_gf)[-1][0] << endl; 
   cout << " norm(my_gf) " << norm( my_gf )  << endl; 

   //// --  Two gf Operators   // generalize such that g<dcomplex> + g<int> possible?
   my_gf += my_gf;     	my_gf + my_gf; 
   my_gf -= my_gf; 	my_gf - my_gf; 
   my_gf *= my_gf;     	my_gf * my_gf; 
   my_gf /= my_gf; 	my_gf / my_gf; 

   //// -- Scalar operators
   my_gf += 1.0;     	my_gf + 1.0;	1.0 + my_gf; 
   my_gf -= 1.0; 	my_gf - 1.0;    1.0 - my_gf;
   my_gf *= 1.0;     	my_gf * 1.0;    1.0 * my_gf;
   my_gf /= 1.0; 	my_gf / 1.0;    

   // -- Test is_instance_of
   cout << " is_instance_of_gf< mygf_t > " << is_instance_of_gf< mygf_t >::value << endl; 
   cout << " is_instance_of_gf< gf_1p_t > " << is_instance_of_gf< gf_1p_t >::value << endl; 

   // -- norm( derived type ) 
   gf_1p_t gf_1p; 
   norm( gf_1p ); 

   // Example for two-particle vertex
   const int POS_FREQ_COUNT_VERT = 10; 
   const int PATCH_COUNT = 4; 
   const int QN_COUNT = 2; 
   enum class I2P{ w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out }; 
   class gf_2p_t : public gf< dcomplex, 10 >              ///< Container type for two-particle correlation functions
   {
      public:
	 gf_2p_t():
	    gf< dcomplex, 10 >( boost::extents[ffreq(POS_FREQ_COUNT_VERT)][ffreq(POS_FREQ_COUNT_VERT)][ffreq(POS_FREQ_COUNT_VERT)]
		  [PATCH_COUNT][PATCH_COUNT][PATCH_COUNT]
		  [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
      {}
   }; 
   using idx_2p_t = gf_2p_t::idx_t;  

}
