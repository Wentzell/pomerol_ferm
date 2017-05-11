
/*******************************************************************************************/ /** @file
 *  		
 * 	file: 		translate.h
 * 	contents:  	Functions that allow for translations between different notations
 * 
 ****************************************************************************************************/

#pragma once

#include <tuple>
#include <complex>
#include <triqs/gfs.hpp>
#include <iostream>

using dcomplex = std::complex<double>;

using triqs::gfs::gf_mesh;
using triqs::gfs::imfreq;
using triqs::gfs::cartesian_product;
using triqs::gfs::matsubara_freq;

using chi4_mesh_t = gf_mesh<triqs::gfs::cartesian_product<triqs::gfs::imfreq, triqs::gfs::imfreq, triqs::gfs::imfreq>>;
using chi3_mesh_t = gf_mesh<triqs::gfs::cartesian_product<imfreq, imfreq>>;
using chi2_mesh_t = gf_mesh<imfreq>;

using chi4_iw_t = std::tuple<triqs::gfs::matsubara_freq, triqs::gfs::matsubara_freq, triqs::gfs::matsubara_freq>;
using chi3_iw_t = std::tuple<triqs::gfs::matsubara_freq, triqs::gfs::matsubara_freq>;
using chi2_iw_t = triqs::gfs::matsubara_freq;

std::ostream &operator<<(std::ostream &out, chi4_iw_t const &idx);
std::ostream &operator<<(std::ostream &out, chi3_iw_t const &idx);

const int n_large = 1e+5;

/********************* Index translations chi3->chi4, chi2->chi4, chi4->chi3, chi3->chi2 ********************/

chi4_iw_t chi3_to_chi4(const chi3_iw_t &idx);
chi4_iw_t chi2_to_chi4(const chi2_iw_t &idx);
chi3_iw_t chi4_to_chi3(const chi4_iw_t &idx);
chi2_iw_t chi3_to_chi2(const chi3_iw_t &idx);

/********************* Translations between different notations  ********************/

struct FERM {};
template <typename notation    = FERM> chi4_iw_t to_ferm(const chi4_iw_t &idx) { return idx; }
template <typename notation    = FERM> chi4_iw_t from_ferm(const chi4_iw_t &idx) { return idx; }
template <typename to_notation = FERM, typename from_notation = FERM> chi4_iw_t translate(const chi4_iw_t &idx) {
  return from_ferm<to_notation>(to_ferm<from_notation>(idx));
}

#define DECLARE_NOTATION(X)                                                                                                                          \
  struct X {};                                                                                                                                       \
                                                                                                                                                     \
  template <> chi4_iw_t to_ferm<X>(const chi4_iw_t &idx);                                                                                            \
                                                                                                                                                     \
  template <> chi4_iw_t from_ferm<X>(const chi4_iw_t &idx);

DECLARE_NOTATION(TRIQS)
DECLARE_NOTATION(POMEROL)
DECLARE_NOTATION(PP)
DECLARE_NOTATION(PH)
DECLARE_NOTATION(XPH)
DECLARE_NOTATION(PP_S)
DECLARE_NOTATION(PH_S)
DECLARE_NOTATION(XPH_S)
