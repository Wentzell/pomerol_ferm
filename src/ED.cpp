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
** \brief Diagonalization of the Anderson impurity model (1 impurity coupled to a set of non-interacting bath sites)
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
*/

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/mpi.hpp>

#include <string>
#include <iostream>
#include <algorithm>

#include <pomerol.h>
#include <triqs/gfs.hpp>
#include <translate.h>

namespace po = boost::program_options;

using namespace Pomerol;
using namespace triqs::gfs;
using namespace triqs::arrays;
using namespace triqs::utility;
using namespace triqs::h5;
using dcomplex = std::complex<double>;

void print_section(const std::string &str, boost::mpi::communicator comm = boost::mpi::communicator());

#define mpi_cout                                                                                                                                     \
  if (!comm.rank()) std::cout

template <typename T> po::options_description define(po::options_description &opts, std::string name, T def_val, std::string desc) {
  opts.add_options()(name.c_str(), po::value<T>()->default_value(T(def_val)), desc.c_str());
  return opts;
}

template <typename T> po::options_description define_vec(po::options_description &opts, std::string name, T &&def_val, std::string desc) {
  opts.add_options()(name.c_str(), po::value<T>()->multitoken()->default_value(T(def_val), ""), desc.c_str());
  return opts;
}

template <typename... A, std::size_t... Is> auto _make_boost_tuple_impl(std::tuple<A...> const &tpl, std::index_sequence<Is...>) {
  return boost::tuple<A...>{std::get<Is>(tpl)...};
}
template <typename... A> auto make_boost_tuple(std::tuple<A...> const &tpl) { return _make_boost_tuple_impl(tpl, std::index_sequence_for<A...>{}); }

// Command Line Parser
po::variables_map cmdline_params(int argc, char *argv[]) {
  po::options_description p("ED Solver based on Pomerol and Triqs");

  define_vec<std::vector<double>>(p, "levels", {}, "Energy levels of the bath sites");
  define_vec<std::vector<double>>(p, "hoppings", {}, "Hopping to the bath sites");

  define<double>(p, "U", 1.0, "Value of U");
  define<double>(p, "beta", 1.0, "Value of inverse temperature");
  define<double>(p, "e0", 0.0, "Value of impurity energy level");

  define<int>(p, "calc_gf", false, "Calculate Green's functions");
  define<int>(p, "calc_2pgf", false, "Calculate 2-particle Green's functions");
  define<int>(p, "n_iw", 32, "Number of positive fermionic Matsubara freq ( 4*n_iw for one-particle GF )");

  define<double>(p, "2pgf.reduce_tol", 1e-12, "Energy resonance resolution in 2pgf");
  define<double>(p, "2pgf.coeff_tol", 1e-12, "Tolerance on nominators in 2pgf");
  define<double>(p, "2pgf.multiterm_tol", 1e-12,
                 "Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms");
  define<size_t>(p, "2pgf.reduce_freq", 1e5, "How often to reduce terms in 2pgf");
  define_vec<std::vector<size_t>>(p, "2pgf.indices", std::vector<size_t>({0, 0, 0, 0}), "2pgf index combination");

  p.add_options()("help", "help");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);

  po::notify(vm);

  if (vm.count("help")) {
    std::cerr << p << "\n";
    MPI_Finalize();
    exit(0);
  }

  return vm;
}

int main(int argc, char *argv[]) {
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator comm;

  print_section("ED Solver based on Pomerol and Triqs");

  // Define Model Parameters
  RealType e0, U, beta;
  bool calc_gf, calc_2pgf;
  size_t L;
  std::vector<double> levels;
  std::vector<double> hoppings;

  // Read command line arguments
  po::variables_map p = cmdline_params(argc, argv);

  std::tie(U, beta, e0, calc_gf, calc_2pgf) =
     std::make_tuple(p["U"].as<double>(), p["beta"].as<double>(), p["e0"].as<double>(), p["calc_gf"].as<int>(), p["calc_2pgf"].as<int>());
  calc_gf = calc_gf || calc_2pgf;

  if (p.count("levels")) {
    levels   = p["levels"].as<std::vector<double>>();
    hoppings = p["hoppings"].as<std::vector<double>>();
  }

  if (levels.size() != hoppings.size()) {
    MPI_Finalize();
    throw(std::logic_error("number of levels != number of hoppings"));
  }

  L = levels.size();
  mpi_cout << "Diagonalization of 1+" << L << " sites" << std::endl;

  Lattice Lat;

  // Add the Impurity Site
  Lat.addSite(new Lattice::Site("A", 1, 2)); // ( Name, Orbital_degrees, Spin_degrees )
  LatticePresets::addCoulombS(&Lat, "A", U, e0);

  // Add the Bath Sites
  std::vector<std::string> names(L);
  for (size_t i = 0; i < L; i++) {
    names[i] = "b" + std::to_string(i);
    Lat.addSite(new Lattice::Site(names[i], 1, 2));
    LatticePresets::addHopping(&Lat, "A", names[i], hoppings[i]);
    LatticePresets::addLevel(&Lat, names[i], levels[i]);
  };

  // Print all Lattice Sites
  int rank = comm.rank();
  mpi_cout << "Sites" << std::endl;
  if (!rank) Lat.printSites();

  // Print all terms with 2 and 4 operators
  if (!rank) {
    mpi_cout << "Terms with 2 operators" << std::endl;
    Lat.printTerms(2);

    mpi_cout << "Terms with 4 operators" << std::endl;
    Lat.printTerms(4);
  };

  // Create index space
  IndexClassification IndexInfo(Lat.getSiteMap());
  IndexInfo.prepare(false);
  if (!rank) {
    print_section("Indices");
    IndexInfo.printIndices();
  };
  int index_size = IndexInfo.getIndexSize();

  print_section("Matrix element storage");
  IndexHamiltonian Storage(&Lat, IndexInfo);
  Storage.prepare();

  // Write down the Hamiltonian as a symbolic formula
  print_section("Terms");
  mpi_cout << Storage << std::endl;

  // Find symmetries of the problem
  Symmetrizer Symm(IndexInfo, Storage);
  Symm.compute();

  // Introduce Fock space and classify states to blocks
  StatesClassification S(IndexInfo, Symm);
  S.compute();

  // Hamiltonian in the basis of Fock Space
  Hamiltonian H(IndexInfo, Storage, S);
  // enter the Hamiltonian matrices
  H.prepare();
  // compute eigenvalues and eigenvectors
  H.compute();

  DensityMatrix rho(S, H, beta); // create Density Matrix
  rho.prepare();
  rho.compute(); // evaluate thermal weights with respect to ground energy, i.e exp(-beta(e-e_0))/Z

  ParticleIndex d0 = IndexInfo.getIndex("A", 0, down); // find the indices of the impurity, i.e. spin up index
  ParticleIndex u0 = IndexInfo.getIndex("A", 0, up);

  if (!rank) {
    // get average total particle number
    mpi_cout << "<N> = " << rho.getAverageOccupancy() << std::endl;
    // get average energy
    mpi_cout << "<H> = " << rho.getAverageEnergy() << std::endl;
    // get double occupancy
    mpi_cout << "<N_{" << IndexInfo.getInfo(u0) << "}N_{" << IndexInfo.getInfo(u0) << "}> = " << rho.getAverageDoubleOccupancy(u0, d0) << std::endl;
    // get average total particle number per index
    for (ParticleIndex i = 0; i < IndexInfo.getIndexSize(); i++) {
      std::cout << "<N_{" << IndexInfo.getInfo(i) << "[" << i << "]}> = " << rho.getAverageOccupancy(i) << std::endl;
    }
    double n_av = rho.getAverageOccupancy();
  }

  // ======== Start Green's function calculation ========

  // Container for c and c^+ in the eigenstate basis
  FieldOperatorContainer Operators(IndexInfo, S, H);

  if (calc_gf) {
    int n_iw = p["n_iw"].as<int>();

    // Open File for output
    triqs::h5::file h5file("pomerol_out.h5", 'w');

    print_section("1-Particle Green's function calculation");
    std::set<ParticleIndex> f;            // a set of indices to evaluate the c and c^+ operators on
    std::set<IndexCombination2> indices2; // a set of index pairs for the Green's function evalulation

    // Take only impurity spin up and spin down indices
    f.insert(u0);
    f.insert(d0);
    indices2.insert(IndexCombination2(d0, d0)); // evaluate only G_{\down \down}

    Operators.prepareAll(f);
    Operators.computeAll(); // evaluate c, c^+ for chosen indices

    GFContainer G(IndexInfo, S, H, rho, Operators);

    G.prepareAll(indices2); // Identify all non-vanishing block connections in the Green's function
    G.computeAll();         // Evaluate all GF terms, i.e. resonances and weights of expressions in Lehmans representation of the Green's function

    // Fill G_iw container and write to file
    if (!comm.rank())
      // Loop over all index pairs of the Green's function
      for (std::set<IndexCombination2>::const_iterator it = indices2.begin(); it != indices2.end(); ++it) {
        IndexCombination2 ind2   = *it;
        const GreensFunction &GF = G(ind2);

        mpi_cout << "Saving imfreq G" << ind2 << " on " << 4 * n_iw << " Matsubara freqs.\n";

        auto G_iw                                   = gf<imfreq, scalar_valued>{{beta, Fermion, 4 * n_iw}};
        for (auto const &iw : G_iw.mesh()) G_iw[iw] = GF(iw);

        std::string ind_str = boost::lexical_cast<std::string>(ind2.Index1) + boost::lexical_cast<std::string>(ind2.Index2);
        h5_write(h5file, "G_iw_" + ind_str, G_iw);
      }

    // Start Calculation of Two-particle Quantities ( chi4, chi3, chi2 )
    if (calc_2pgf) {
      print_section("2-Particle Green's function calculation (chi4)");

      std::vector<size_t> indices_2pgf = p["2pgf.indices"].as<std::vector<size_t>>();
      if (indices_2pgf.size() != 4) throw std::logic_error("Need 4 indices for 2pgf");

      // a set of four indices to evaluate the 2pgf
      IndexCombination4 index_comb(indices_2pgf[0], indices_2pgf[1], indices_2pgf[2], indices_2pgf[3]);

      std::set<IndexCombination4> indices4;
      // 2PGF = <T c c c^+ c^+>
      indices4.insert(index_comb);
      std::string ind_str = boost::lexical_cast<std::string>(index_comb.Index1) + boost::lexical_cast<std::string>(index_comb.Index2)
         + boost::lexical_cast<std::string>(index_comb.Index3) + boost::lexical_cast<std::string>(index_comb.Index4);

      AnnihilationOperator const &C1 = Operators.getAnnihilationOperator(index_comb.Index1);
      AnnihilationOperator const &C2 = Operators.getAnnihilationOperator(index_comb.Index2);
      CreationOperator const &CX3    = Operators.getCreationOperator(index_comb.Index3);
      CreationOperator const &CX4    = Operators.getCreationOperator(index_comb.Index4);
      TwoParticleGF G4(S, H, C1, C2, CX3, CX4, rho);

      /* Some knobs to make calc faster - the larger the values of tolerances, the faster is calc, but rounding errors may show up. */
      /** A difference in energies with magnitude less than this value is treated as zero - resolution of energy resonances. */
      G4.ReduceResonanceTolerance = p["2pgf.reduce_tol"].as<double>();
      /** Minimal magnitude of the coefficient of a term to take it into account - resolution of thermal weight. */
      G4.CoefficientTolerance = p["2pgf.coeff_tol"].as<double>();
      /** Knob that controls the caching frequency. */
      G4.ReduceInvocationThreshold = p["2pgf.reduce_freq"].as<size_t>();
      /** Minimal magnitude of the coefficient of a term to take it into account with respect to amount of terms. */
      G4.MultiTermCoefficientTolerance = p["2pgf.multiterm_tol"].as<double>();

      G4.prepare();
      comm.barrier(); // MPI::BARRIER

      // =============== CHI4 ===============

      // Prepare frequency triples to be calculated
      std::vector<boost::tuple<dcomplex, dcomplex, dcomplex>> freqs_2pgf;
      freqs_2pgf.reserve(8 * n_iw * n_iw * n_iw);

      auto iw_mesh = gf_mesh<imfreq>{beta, Fermion, n_iw};

      for (auto const &iw1 : iw_mesh)
        for (auto const &iw2 : iw_mesh)
          for (auto const &iw3 : iw_mesh) freqs_2pgf.push_back(make_boost_tuple(translate<POMEROL, FERM>({iw1, iw3, iw2})));

      mpi_cout << "2PGF : " << freqs_2pgf.size() << " freqs to evaluate" << std::endl;

      // Calculate chi4 based on the frequency triples
      auto chi4_freq_data = G4.compute(true, freqs_2pgf, comm);

      // Initialize chi4 container and write to file
      if (!comm.rank()) {
        mpi_cout << "Saving 2PGF " << index_comb << std::endl;
        auto chi4  = gf<cartesian_product<imfreq, imfreq, imfreq>, scalar_valued>{{iw_mesh, iw_mesh, iw_mesh}};
        size_t idx = 0;
        for (auto const &iw1 : iw_mesh)
          for (auto const &iw2 : iw_mesh)
            for (auto const &iw3 : iw_mesh) chi4[{iw1, iw2, iw3}] = chi4_freq_data[idx++];

        h5_write(h5file, "chi4_" + ind_str, chi4);
      }

      //// =============== CHI3 ===============

      //print_section("High-frequency asymptotics (chi3)");

      //// Prepare frequency triples to be calculated
      //std::vector<boost::tuple<dcomplex, dcomplex, dcomplex>> freqs_chi3_pp, freqs_chi3_ph, freqs_chi3_xph;
      //freqs_chi3_pp.reserve(2 * n_iw * (2 * n_iw + 1));
      //freqs_chi3_ph.reserve(2 * n_iw * (2 * n_iw + 1));
      //freqs_chi3_xph.reserve(2 * n_iw * (2 * n_iw + 1));

      //auto iOm_mesh = gf_mesh<imfreq>{beta, Boson, n_iw};

      //for (auto const &iOm : iOm_mesh)
        //for (auto const &iw : iw_mesh) {
          //freqs_chi3_pp.push_back(make_boost_tuple(translate<POMEROL, PP>(chi3_to_chi4({iOm, iw}))));
          //freqs_chi3_ph.push_back(make_boost_tuple(translate<POMEROL, PH>(chi3_to_chi4({iOm, iw}))));
          //freqs_chi3_xph.push_back(make_boost_tuple(translate<POMEROL, XPH>(chi3_to_chi4({iOm, iw}))));
        //}

      //// Calculate chi3 in all channels based on the frequency triples
      //auto chi3_pp_freq_data  = G4.compute(true, freqs_chi3_pp, comm);
      //auto chi3_ph_freq_data  = G4.compute(true, freqs_chi3_ph, comm);
      //auto chi3_xph_freq_data = G4.compute(true, freqs_chi3_xph, comm);

      //// Initialize chi3 containers and write to file
      //if (!comm.rank()) {
        //mpi_cout << "Saving chi3 " << index_comb << std::endl;

        //auto chi3_pp_iw  = gf<cartesian_product<imfreq, imfreq>, scalar_valued>{{iOm_mesh, iw_mesh}};
        //auto chi3_ph_iw  = gf<cartesian_product<imfreq, imfreq>, scalar_valued>{{iOm_mesh, iw_mesh}};
        //auto chi3_xph_iw = gf<cartesian_product<imfreq, imfreq>, scalar_valued>{{iOm_mesh, iw_mesh}};

        //size_t idx = 0;

        //for (auto const &iOm : iOm_mesh)
          //for (auto const &iw : iw_mesh) {
            //chi3_pp_iw[{iOm, iw}]  = chi3_pp_freq_data[idx];
            //chi3_ph_iw[{iOm, iw}]  = chi3_ph_freq_data[idx];
            //chi3_xph_iw[{iOm, iw}] = chi3_xph_freq_data[idx++];
          //}

        //h5_write(h5file, "chi3_pp_iw_" + ind_str, chi3_pp_iw);
        //h5_write(h5file, "chi3_ph_iw_" + ind_str, chi3_ph_iw);
        //h5_write(h5file, "chi3_xph_iw_" + ind_str, chi3_xph_iw);
      //}

      //// =============== CHI2 ===============

      //print_section("High-frequency asymptotics (chi2)");

      //// Prepare frequency triples to be calculated
      //std::vector<boost::tuple<dcomplex, dcomplex, dcomplex>> freqs_chi2_pp, freqs_chi2_ph, freqs_chi2_xph;
      //freqs_chi3_pp.reserve(2 * n_iw + 1);
      //freqs_chi3_ph.reserve(2 * n_iw + 1);
      //freqs_chi3_xph.reserve(2 * n_iw + 1);

      //for (auto const &iOm : iOm_mesh) {
        //freqs_chi2_pp.push_back(make_boost_tuple(translate<POMEROL, PP>(chi2_to_chi4(iOm))));
        //freqs_chi2_ph.push_back(make_boost_tuple(translate<POMEROL, PH>(chi2_to_chi4(iOm))));
        //freqs_chi2_xph.push_back(make_boost_tuple(translate<POMEROL, XPH>(chi2_to_chi4(iOm))));
      //}

      //// The most important routine - actually calculate the 2PGF
      //auto chi2_pp_freq_data  = G4.compute(true, freqs_chi2_pp, comm);
      //auto chi2_ph_freq_data  = G4.compute(true, freqs_chi2_ph, comm);
      //auto chi2_xph_freq_data = G4.compute(true, freqs_chi2_xph, comm);

      //// Initialize chi3 containers and write to file
      //if (!comm.rank()) {
        //mpi_cout << "Saving chi2 " << index_comb << std::endl;

        //auto chi2_pp_iw  = gf<imfreq, scalar_valued>{iOm_mesh};
        //auto chi2_ph_iw  = gf<imfreq, scalar_valued>{iOm_mesh};
        //auto chi2_xph_iw = gf<imfreq, scalar_valued>{iOm_mesh};

        //size_t idx = 0;

        //for (auto const &iOm : iOm_mesh) {
          //chi2_pp_iw[iOm]  = chi2_pp_freq_data[idx];
          //chi2_ph_iw[iOm]  = chi2_ph_freq_data[idx];
          //chi2_xph_iw[iOm] = chi2_xph_freq_data[idx];
        //}

        //h5_write(h5file, "chi2_pp_iw_" + ind_str, chi2_pp_iw);
        //h5_write(h5file, "chi2_ph_iw_" + ind_str, chi2_ph_iw);
        //h5_write(h5file, "chi2_xph_iw_" + ind_str, chi2_xph_iw);
      //}
    }
  }
}

void print_section(const std::string &str, boost::mpi::communicator comm) {
  mpi_cout << std::string(str.size(), '=') << std::endl;
  mpi_cout << str << std::endl;
  mpi_cout << std::string(str.size(), '=') << std::endl;
}
