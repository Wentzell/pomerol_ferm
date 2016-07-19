//
// This file is a part of pomerol - a scientific ED code for obtaining 
// properties of a Hubbard model on a finite-size lattice 
//
// Copyright (C) 2010-2014 Andrey Antipov <Andrey.E.Antipov@gmail.com>
// Copyright (C) 2010-2014 Igor Krivenko <Igor.S.Krivenko@gmail.com>
//
// pomerol is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pomerol is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with pomerol.  If not, see <http://www.gnu.org/licenses/>.


/** \file prog/anderson.cpp
 ** \brief Diagonalization of the Hubbard 2d cluster with periodic boundary conditions
 **
 ** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
 */

#pragma clang diagnostic ignored "-Wc++11-extensions"
#pragma clang diagnostic ignored "-Wgnu"

#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>

#include <cstdlib>
#include <fstream>

#include <pomerol.h>
#include "mpi_dispatcher/mpi_dispatcher.hpp"

#include <def.h>
#include <H5Tools.h>
#include <grid.h>

using namespace Pomerol;
using namespace H5;

extern boost::mpi::environment env;
boost::mpi::communicator comm;

/* Auxiliary routines - implemented in the bottom. */
bool compare(ComplexType a, ComplexType b);
void print_section (const std::string& str);
template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance = 1e-7);
template <typename T1> void savetxt(std::string fname, T1 in);
struct my_logic_error;
double FMatsubara(int n, double beta){return M_PI/beta*(2.*n+1);}
double BMatsubara(int n, double beta){return M_PI/beta*(2.*n);}

int job_to_ffreq_index(int job, int wfmax) { return -wfmax + job; }

int main(int argc, char* argv[])
{
   boost::mpi::environment env(argc,argv);
   boost::mpi::communicator comm;

   print_section("Dimer Model");

   int size_x, size_y;
   RealType t, mu, U, beta, reduce_tol, coeff_tol;
   bool calc_gf, calc_2pgf, cluster; 
   int wf_max;
   double eta, hbw, step; // for evaluation of GF on real axis 

   try { // command line parser
      TCLAP::CmdLine cmd("Hubbard nxn diag", ' ', "");
      TCLAP::ValueArg<RealType> U_arg("U","U","Value of U",true,10.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> mu_arg("","mu","Global chemical potential",false,0.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> t_arg("t","t","Value of t",false,1.0,"RealType",cmd);

      TCLAP::ValueArg<RealType> beta_arg("b","beta","Inverse temperature",true,100.,"RealType");
      TCLAP::ValueArg<RealType> T_arg("T","T","Temperature",true,0.01,"RealType");
      cmd.xorAdd(beta_arg,T_arg);

      TCLAP::ValueArg<size_t> wf_arg("","wf","Number of positive fermionic Matsubara Freqs",false,64,"int",cmd);

      TCLAP::SwitchArg gf_arg("","calcgf","Calculate Green's functions",cmd, false);
      TCLAP::SwitchArg twopgf_arg("","calc2pgf","Calculate 2-particle Green's functions",cmd, false);
      TCLAP::ValueArg<RealType> reduce_tol_arg("","reducetol","Energy resonance resolution in 2pgf",false,1e-12,"RealType",cmd);
      TCLAP::ValueArg<RealType> coeff_tol_arg("","coefftol","Total weight tolerance",false,1e-12,"RealType",cmd);

      cluster = true; 

      TCLAP::ValueArg<RealType> eta_arg("","eta","Offset from the real axis for Green's function calculation",false,0.05,"RealType",cmd);
      TCLAP::ValueArg<RealType> hbw_arg("D","hbw","Half-bandwidth. Default = U",false,0.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> step_arg("","step","Step on a real axis. Default : 0.01",false,0.01,"RealType",cmd);

      cmd.parse( argc, argv ); // parse arguments

      U = U_arg.getValue();
      mu = (mu_arg.isSet()?mu_arg.getValue():U/2);
      boost::tie(t, beta, calc_gf, calc_2pgf, reduce_tol, coeff_tol) = boost::make_tuple( t_arg.getValue(), beta_arg.getValue(), 
	    gf_arg.getValue(), twopgf_arg.getValue(), reduce_tol_arg.getValue(), coeff_tol_arg.getValue());
      size_x = 2; size_y = 1; 
      boost::tie(wf_max) = boost::make_tuple(wf_arg.getValue());
      boost::tie(eta, hbw, step) = boost::make_tuple(eta_arg.getValue(), (hbw_arg.isSet()?hbw_arg.getValue():2.*U), step_arg.getValue());
      calc_gf = calc_gf || calc_2pgf;
   }
   catch (TCLAP::ArgException &e)  // catch parsing exceptions
   { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; exit(1);}
   catch (std::exception &e)  // catch standard exceptions
   { std::cerr << "error: " << e.what() << std::endl; exit(1);}

   int L = size_x*size_y;
   INFO("Diagonalization of " << L << "=" << size_x << "*" << size_y << " sites");
   Lattice Lat;

   /* Add sites */
   std::vector<std::string> names(L);
   auto SiteIndexF = [size_x](size_t x, size_t y){return y*size_x + x;};
   for (size_t y=0; y<size_y; y++)
      for (size_t x=0; x<size_x; x++)
      {
	 int i = SiteIndexF(x,y);
	 names[i]="S"+std::to_string(i);
	 Lat.addSite(new Lattice::Site(names[i],1,2));
      };

   INFO("Sites");
   Lat.printSites();

   /* Add hopping */
   for (size_t x=0; x<size_x; x++) 
      for (size_t y=0; y<size_y; y++) 
      {
	 auto pos = SiteIndexF(x,y);
	 auto pos_right = SiteIndexF( (x+1) % size_x, y ); 
	 auto pos_up = SiteIndexF( x, (y+1) % size_y ); 

	 if (x < size_x-1) LatticePresets::addHopping(&Lat, names[pos], names[pos_right], -t);
	 else if(!cluster) LatticePresets::addHopping(&Lat, names[pos_right], names[pos], -t); // periodicity in x
	    
	 if (y < size_y-1) LatticePresets::addHopping(&Lat, names[pos], names[pos_up], -t);
	 else if(!cluster) LatticePresets::addHopping(&Lat, names[pos_up], names[pos], -t); // periodicity in y
      }

   /* Add interaction on each site*/
   for (size_t i=0; i<L; i++) LatticePresets::addCoulombS(&Lat, names[i], U, -mu);

   auto rank = comm.rank();
   if (!rank) {
      INFO("Terms with 2 operators");
      Lat.printTerms(2);

      INFO("Terms with 4 operators");
      Lat.printTerms(4);
   };

   IndexClassification IndexInfo(Lat.getSiteMap());
   IndexInfo.prepare(false); // Create index space
   if (!rank) { print_section("Indices"); IndexInfo.printIndices(); };
   auto index_size = IndexInfo.getIndexSize();

   print_section("Matrix element storage");
   IndexHamiltonian Storage(&Lat,IndexInfo); 
   Storage.prepare(); // Write down the Hamiltonian as a symbolic formula
   print_section("Terms");
   if (!rank) INFO(Storage);

   Symmetrizer Symm(IndexInfo, Storage);
   Symm.compute(); // Find symmetries of the problem

   StatesClassification S(IndexInfo,Symm); // Introduce Fock space and classify states to blocks
   S.compute();

   Hamiltonian H(IndexInfo, Storage, S); // Hamiltonian in the basis of Fock Space
   H.prepare(); // enter the Hamiltonian matrices
   H.compute(); // compute eigenvalues and eigenvectors

   RealVectorType evals (H.getEigenValues());
   std::sort(evals.data(), evals.data() + H.getEigenValues().size());
   savetxt("spectrum.dat", evals); // dump eigenvalues

   DensityMatrix rho(S,H,beta); // create Density Matrix
   rho.prepare();
   rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z 

   // Register indeces for preparation (f)
   std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
   for( int i = 0; i < IndexInfo.getIndexSize(); i++ )
      f.insert(i);

   // Green's function calculation starts here
   FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

   Operators.prepareAll(f);
   Operators.computeAll(); // evaluate c, c^+ for chosen indices 

   ParticleIndex d0 = IndexInfo.getIndex("S0",0,down); 
   ParticleIndex u0 = IndexInfo.getIndex("S0",0,up);

   INFO("<N> = " << rho.getAverageOccupancy()); // get average total particle number
   INFO("<H> = " << rho.getAverageEnergy()); // get average energy
   INFO("<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0)); // get double occupancy
   for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++)
      INFO("<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i)); // get average total particle number
   savetxt("N_T.dat",rho.getAverageOccupancy());

   if (calc_gf) {
      INFO("1-particle Green's functions calc");
      std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

      // Register indeces for 1p Green function
      for ( int i = 0; i < L; i++ )
	 for ( int j = 0; j < L; j++ )
	 {
	    ParticleIndex ind_i = IndexInfo.getIndex(names[i],0,up);
	    ParticleIndex ind_j = IndexInfo.getIndex(names[j],0,up);
	    indices2.insert(IndexCombination2(ind_i,ind_j));
	 }

      GFContainer G(IndexInfo,S,H,rho,Operators);

      G.prepareAll(indices2); // identify all non-vanishing block connections in the Green's function
      G.computeAll(); // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

      gf_1p_t Giw( 4*wf_max, L ); 

      if (!comm.rank()) // dump gf into a file
      {
	 std::cout << "Saving imfreq G on "<< 4*wf_max << " Matsubara freqs. " << std::endl;
	 for (auto ind2 : indices2) // loops over all components (pairs of indices) of the Green's function 
	 { 
	    const GreensFunction & GF = G(ind2);

	    // Save Matsubara GF on grid
	    for (int wn = -4*wf_max; wn < 4*wf_max; wn++)	// Care for spin index!
	       Giw[wn][0][ind2.Index1 / 2][ind2.Index2 / 2] = GF(I*FMatsubara(wn,beta)); 
	 }
      }

      H5File file( "output.h5", H5F_ACC_TRUNC );
      Group grp_par( file.createGroup("/Params"));
      write( U, grp_par, "UINT" ); 
      write( beta, grp_par, "BETA" ); 

      Group grp_giw( file.createGroup("/Giw") );
      write( Giw, grp_giw, "" ); 
      write( F_Grid( 4*wf_max, 2.0*PI / beta ), grp_giw );

      // Start Two-particle GF calculation

      if (calc_2pgf) {   
	 print_section("2-Particle Green's function calc");
	 std::set<IndexCombination4> indices4; // a set of four indices to evaluate the 2pgf
	 // 2PGF = <T c c c^+ c^+>      i j k l
	 for ( int i = 0; i < L; i++ )
	    for ( int j = 0; j < L; j++ )
	       for ( int k = 0; k < L; k++ )
		  for ( int l = 0; l < L; l++ )
		  {
		     ParticleIndex ind_i = IndexInfo.getIndex(names[i],0,up);
		     ParticleIndex ind_j = IndexInfo.getIndex(names[j],0,down);
		     ParticleIndex ind_k = IndexInfo.getIndex(names[k],0,down);
		     ParticleIndex ind_l = IndexInfo.getIndex(names[l],0,up);
		     indices4.insert(IndexCombination4(ind_i, ind_j, ind_k, ind_l));
		  }

	 TwoParticleGFContainer Chi4(IndexInfo,S,H,rho,Operators);
	 /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
	 /**  Kronecker delta comparing two energies ... default machine precision */
	 Chi4.KroneckerSymbolTolerance = 1e-14;
	 /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
	 Chi4.ReduceResonanceTolerance = reduce_tol;
	 /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
	 Chi4.CoefficientTolerance = coeff_tol;
	 /** Knob that controls the caching frequency. */
	 Chi4.ReduceInvocationThreshold = 1e5;
	 /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
	 Chi4.MultiTermCoefficientTolerance = 1e-12;

	 Chi4.prepareAll(indices4); // find all non-vanishing block connections inside 2pgf

	 bool clearTerms = false; // Set if terms are kept in memory NOT YET WORKING

	 // ! The most important routine - actually calculate the 2PGF
	 Chi4.computeAll(clearTerms, comm, true); 

	 gf_2p_t G2arr( wf_max, L ); 

	 for (auto ind : indices4) { // dump 2PGF into files - loop through 2pgf components
	    std::string ind_str = std::to_string(ind.Index1 / 2) + std::to_string(ind.Index2 / 2) +std::to_string(ind.Index4 / 2) +std::to_string(ind.Index3 / 2);
	    if (!comm.rank()) std::cout << "Saving 2PGF " << ind_str << std::endl;
	    const TwoParticleGF &chi = Chi4(ind);

	    for (int w1_index = -wf_max; w1_index<wf_max; w1_index++) // loop over second fermionic frequency 
	       for (int w2_index = -wf_max; w2_index<wf_max; w2_index++) // loop over second fermionic frequency 
		  for (int w1p_index = -wf_max; w1p_index<wf_max; w1p_index++) 
		  { 
		     int w2p_index = w1_index + w2_index - w1p_index;
		     G2arr[w1_index][w2_index][w1p_index][0][0][0][ind.Index1 / 2][ind.Index2 / 2][ind.Index4 / 2][ind.Index3 / 2] = chi(I*FMatsubara( w1_index, beta ), I*FMatsubara( w2_index, beta ), I*FMatsubara( w2p_index, beta ));  
		  }
	 }

	 Group grp_G2( file.createGroup("/G2") );
	 write( G2arr, grp_G2, "" ); 
	 write( F_Grid( wf_max, 2.0*PI / beta ), grp_G2 );

      }
   }
}

bool compare(ComplexType a, ComplexType b)
{
   return abs(a-b) < 1e-5;
}

void print_section (const std::string& str)
{
   if (!comm.rank()) { 
      std::cout << std::string(str.size(),'=') << std::endl;
      std::cout << str << std::endl;
      std::cout << std::string(str.size(),'=') << std::endl;
   };
}

   template <typename F1, typename F2>
bool is_equal ( F1 x, F2 y, RealType tolerance)
{
   return (std::abs(x-y)<tolerance);
}

template <typename T1> void savetxt(std::string fname, T1 in){std::ofstream out(fname); out << in << std::endl; out.close();};

struct my_logic_error : public std::logic_error { my_logic_error (const std::string& what_arg):logic_error(what_arg){}; };
