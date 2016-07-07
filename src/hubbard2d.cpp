
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

#include<cstdlib>
#include <fstream>

#include <pomerol.h>
#include "mpi_dispatcher/mpi_dispatcher.hpp"

using namespace Pomerol;


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

   print_section("Hubbard nxn");

   int size_x, size_y, wn;
   RealType t, mu, U, beta, reduce_tol, coeff_tol;
   bool calc_gf, calc_2pgf, cluster; 
   int wf_max, wb_max;
   double eta, hbw, step; // for evaluation of GF on real axis 

   try { // command line parser
      TCLAP::CmdLine cmd("Hubbard nxn diag", ' ', "");
      TCLAP::ValueArg<RealType> U_arg("U","U","Value of U",true,10.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> mu_arg("","mu","Global chemical potential",false,0.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> t_arg("t","t","Value of t",false,1.0,"RealType",cmd);

      TCLAP::ValueArg<RealType> beta_arg("b","beta","Inverse temperature",true,100.,"RealType");
      TCLAP::ValueArg<RealType> T_arg("T","T","Temperature",true,0.01,"RealType");
      cmd.xorAdd(beta_arg,T_arg);

      TCLAP::ValueArg<size_t> x_arg("x","x","Size over x",false,2,"int",cmd);
      TCLAP::ValueArg<size_t> y_arg("y","y","Size over y",false,2,"int",cmd);

      TCLAP::ValueArg<size_t> wn_arg("","wf","Number of positive fermionic Matsubara Freqs",false,64,"int",cmd);

      TCLAP::SwitchArg gf_arg("","calcgf","Calculate Green's functions",cmd, false);
      TCLAP::SwitchArg twopgf_arg("","calc2pgf","Calculate 2-particle Green's functions",cmd, false);
      TCLAP::ValueArg<RealType> reduce_tol_arg("","reducetol","Energy resonance resolution in 2pgf",false,1e-12,"RealType",cmd);
      TCLAP::ValueArg<RealType> coeff_tol_arg("","coefftol","Total weight tolerance",false,1e-12,"RealType",cmd);

      TCLAP::SwitchArg cluster_arg("","cluster","Disable periodic boundrary conditions",cmd, false);

      TCLAP::ValueArg<RealType> eta_arg("","eta","Offset from the real axis for Green's function calculation",false,0.05,"RealType",cmd);
      TCLAP::ValueArg<RealType> hbw_arg("D","hbw","Half-bandwidth. Default = U",false,0.0,"RealType",cmd);
      TCLAP::ValueArg<RealType> step_arg("","step","Step on a real axis. Default : 0.01",false,0.01,"RealType",cmd);

      cmd.parse( argc, argv ); // parse arguments

      U = U_arg.getValue();
      mu = (mu_arg.isSet()?mu_arg.getValue():U/2);
      boost::tie(t, beta, calc_gf, calc_2pgf, cluster, reduce_tol, coeff_tol) = boost::make_tuple( t_arg.getValue(), beta_arg.getValue(), 
	    gf_arg.getValue(), twopgf_arg.getValue(), cluster_arg.getValue(), reduce_tol_arg.getValue(), coeff_tol_arg.getValue());
      boost::tie(size_x, size_y) = boost::make_tuple(x_arg.getValue(), y_arg.getValue());
      boost::tie(wf_max) = boost::make_tuple(wn_arg.getValue());
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
	 auto i = SiteIndexF(x,y);
	 std::stringstream s; s << i; 
	 names[i]="S"+s.str();
	 Lat.addSite(new Lattice::Site(names[i],1,2));
      };

   INFO("Sites");
   Lat.printSites();

   /* Add interaction on each site*/
   for (size_t i=0; i<L; i++) LatticePresets::addCoulombS(&Lat, names[i], U, -mu);

   /* Add hopping */
   for (size_t y=0; y<size_y-1; y++) {
      for (size_t x=0; x<size_x-1; x++) {
	 auto pos = SiteIndexF(x,y);
	 auto pos_right = SiteIndexF(x+1,y); 
	 auto pos_up = SiteIndexF(x,y+1); 

	 if (size_x > 1) LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_right]), std::max(names[pos], names[pos_right]), -t);
	 if (size_y > 1) LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_up]), std::max(names[pos], names[pos_up]), -t);
      };

      if (size_x > 1 && !cluster ) // periodic boundrary conditions in x direction
      {
	 auto pos = SiteIndexF(size_x-1,y);
	 auto pos_right = SiteIndexF(0,y); 

	 LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_right]), std::max(names[pos], names[pos_right]), -t); 
      }
   }

   if (size_y > 1 && !cluster ) // periodic boundrary conditions in y direction
      for (size_t x=0; x<size_x; x++) 
      {
	 auto pos = SiteIndexF(x,size_y-1);
	 auto pos_up = SiteIndexF(x,0); 

	 LatticePresets::addHopping(&Lat, std::min(names[pos], names[pos_up]), std::max(names[pos], names[pos_up]), -t); 
      }

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

   INFO("<N> = " << rho.getAverageOccupancy()); // get average total particle number
   INFO("<H> = " << rho.getAverageEnergy()); // get average energy
   ParticleIndex d0 = IndexInfo.getIndex("S0",0,down); // find the indices of the impurity, i.e. spin up index
   ParticleIndex u0 = IndexInfo.getIndex("S0",0,up);
   INFO("<N_{" << IndexInfo.getInfo(u0) << "}N_{"<< IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0,d0)); // get double occupancy

   for (ParticleIndex i=0; i<IndexInfo.getIndexSize(); i++) {  
      INFO("<N_{" << IndexInfo.getInfo(i) << "[" << i <<"]}> = " << rho.getAverageOccupancy(i)); // get average total particle number
   }

   savetxt("N_T.dat",rho.getAverageOccupancy());

   // Green's function calculation starts here

   FieldOperatorContainer Operators(IndexInfo, S, H); // Create a container for c and c^+ in the eigenstate basis

   if (calc_gf) {
      INFO("1-particle Green's functions calc");
      std::set<ParticleIndex> f; // a set of indices to evaluate c and c^+
      std::set<IndexCombination2> indices2; // a set of pairs of indices to evaluate Green's function

      // Register indeces for preparation (f) and Green function calculation (indeces2)
      f.insert(u0);
      f.insert(d0);
      for (size_t x=0; x<size_x; x++)
      {
	 ParticleIndex ind = IndexInfo.getIndex(names[SiteIndexF(x,0)],0,down);
	 f.insert(ind); 
	 indices2.insert(IndexCombination2(d0,ind));
      }

      Operators.prepareAll(f);
      Operators.computeAll(); // evaluate c, c^+ for chosen indices 

      GFContainer G(IndexInfo,S,H,rho,Operators);

      G.prepareAll(indices2); // identify all non-vanishing block connections in the Green's function
      G.computeAll(); // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

      if (!comm.rank()) // dump gf into a file
	 for (auto ind2 : indices2) { // loops over all components (pairs of indices) of the Green's function 
	    // Save Matsubara GF from pi/beta to pi/beta*(4*wf_max + 1)
	    std::cout << "Saving imfreq G" << ind2 << " on "<< 4*wf_max << " Matsubara freqs. " << std::endl;
	    std::ofstream gw_im("gw_imag"+std::to_string(ind2.Index1)+std::to_string(ind2.Index2)+".dat");
	    const GreensFunction & GF = G(ind2);
	    for (int wn = 0; wn < wf_max*4; wn++) {
	       ComplexType val = GF(I*FMatsubara(wn,beta)); // this comes from Pomerol - see GreensFunction::operator()
	       gw_im << std::scientific << std::setprecision(12) << FMatsubara(wn,beta) << "   " << real(val) << " " << imag(val) << std::endl;
	    };
	    gw_im.close();
	    // Save Retarded GF on the real axis
	    double e0 = U - 2.*mu;
	    std::ofstream gw_re("gw_real"+std::to_string(ind2.Index1)+std::to_string(ind2.Index2)+".dat"); 
	    std::cout << "Saving real-freq GF " << ind2 << " in energy space [" << e0-hbw << ":" << e0+hbw << ":" << step << "] + I*" << eta << "." << std::endl;
	    for (double w = e0-hbw; w < e0+hbw; w+=step) {
	       ComplexType val = GF(ComplexType(w) + I*eta);
	       gw_re << std::scientific << std::setprecision(12) << w << "   " << real(val) << " " << imag(val) << std::endl;
	    };
	    gw_re.close();
	 }

      // Start Two-particle GF calculation

      if (calc_2pgf) {   
	 print_section("2-Particle Green's function calc");
	 std::set<IndexCombination4> indices4; // a set of four indices to evaluate the 2pgf
	 // 2PGF = <T c c c^+ c^+>
	 indices4.insert(IndexCombination4(u0,d0,d0,u0)); // register 1010 in my notation

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
	 comm.barrier(); // MPI::BARRIER
	 bool clearTerms = false; // Set if terms are kept in memory NOT YET WORKING

	 // ! The most important routine - actually calculate the 2PGF
	 Chi4.computeAll(clearTerms, comm, true); 

	 for (auto ind : indices4) { // dump 2PGF into files - loop through 2pgf components
	    if (!comm.rank()) std::cout << "Saving 2PGF " << ind << std::endl;
	    std::string ind_str = std::to_string(ind.Index1) + std::to_string(ind.Index2) +std::to_string(ind.Index4) +std::to_string(ind.Index3);
	    const TwoParticleGF &chi = Chi4(ind);

	    // Save terms of two particle GF
	    //std::ofstream term_res_stream("terms_res"+ind_str+".pom");
	    //std::ofstream term_nonres_stream("terms_nonres"+ind_str+".pom");
	    //boost::archive::text_oarchive oa_res(term_res_stream);
	    //boost::archive::text_oarchive oa_nonres(term_nonres_stream);
	    //for(std::vector<TwoParticleGFPart*>::const_iterator iter = chi.parts.begin(); iter != chi.parts.end(); iter++) {
	       //oa_nonres << ((*iter)->getNonResonantTerms());
	       //oa_res << ((*iter)->getResonantTerms());
	    //};

	    // start output of vertex
	    // wf_max -  number of positive fermionic freqs
	    // give 2pgf in bosonic-fermionic-fermionic freq notation

	    { // dispatch and save two-particle GF data - MPI parallelization in first fermionic freq

	       // Master-slave scheme to distribute the first fermionic frequency on different processes
	       int ROOT = 0;
	       int ntasks = 2*wf_max;

	       std::unique_ptr<pMPI::MPIMaster> disp;

	       if (comm.rank() == ROOT) { 
		  DEBUG("Master at " << comm.rank());
		  disp.reset(new pMPI::MPIMaster(comm,ntasks,true));
	       };
	       comm.barrier();


	       for (pMPI::MPIWorker worker(comm,ROOT);!worker.is_finished();) {
		  if (rank == ROOT) disp->order(); 
		  worker.receive_order(); 
		  //DEBUG((worker.Status == WorkerTag::Pending));
		  if (worker.is_working()) { 
		     // this is what every process executes
		     auto p = worker.current_job();

		     int w1_index = job_to_ffreq_index(p, wf_max); 
		     std::cout << "["<<p+1<<"/" << ntasks << "] p" << comm.rank() << " w1_idx = " << w1_index << std::endl;

		     std::ofstream chi_stream ("chi_"+ind_str+"_w1_"+std::to_string(p)+".dat");

		     // Most important part 2 - loop over fermionic frequencies. Consider parallelizing one of the loop
		     for (int w2_index = -wf_max; w2_index<wf_max; w2_index++) { // loop over second fermionic frequency 
			for (int w1p_index = -wf_max; w1p_index<wf_max; w1p_index++) { // loop over third fermionic
			   int w2p_index = w1_index + w2_index - w1p_index;
			   ComplexType val = chi(I*FMatsubara( w1_index, beta ), I*FMatsubara( w2_index, beta ), I*FMatsubara( w2p_index, beta ));  
			   chi_stream << std::scientific << std::setprecision(12) 
			      << w1_index << " " << w2_index << " " << w1p_index << "   " << std::real(val) << " " << std::imag(val) << std::endl;
			}
			chi_stream << std::endl;
		     };
		     chi_stream.close();
		     worker.report_job_done(); 
		  };
		  if (rank == ROOT) disp->check_workers(); // check if there are free workers 
	       };
	       comm.barrier();
	       if (comm.rank() == ROOT) { disp.release(); DEBUG("Released master"); };
	    }
	 }
      };
   };
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



