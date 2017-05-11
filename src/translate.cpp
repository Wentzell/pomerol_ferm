
/*******************************************************************************************/ /** @file
 *  		
 * 	file: 		translate.cpp
 * 	contents:  	Functions that allow for translations between different notations
 * 
 ****************************************************************************************************/

#include <translate.h>
#include <cmath>

// Reference Fermionic notation
// ( w1_in, w2_in, w1_out )

// Idx structure for mixed notations
// ( W, w_in, w_out )

using namespace triqs::gfs;

chi4_iw_t chi3_to_chi4(const chi3_iw_t &iw_tpl) {
  return {std::get<0>(iw_tpl), std::get<1>(iw_tpl), matsubara_freq(n_large, std::get<0>(iw_tpl).beta, Fermion)};
}

chi4_iw_t chi2_to_chi4(const chi2_iw_t &iw) {
  auto iw_large = matsubara_freq(n_large, iw.beta, Fermion);
  return {iw, iw_large, iw_large};
}

chi3_iw_t chi4_to_chi3(const chi4_iw_t &iw_tpl) { return {std::get<0>(iw_tpl), std::get<1>(iw_tpl)}; }

chi2_iw_t chi3_to_chi2(const chi3_iw_t &iw_tpl) { return std::get<0>(iw_tpl); }

std::ostream &operator<<(std::ostream &out, chi4_iw_t const &iw_tpl) {
  return out << std::get<0>(iw_tpl) << " " << std::get<1>(iw_tpl) << " " << std::get<2>(iw_tpl) << "\n";
}

std::ostream &operator<<(std::ostream &out, chi3_iw_t const &iw_tpl) { return out << std::get<0>(iw_tpl) << " " << std::get<1>(iw_tpl) << "\n"; }

// =============== Notation TRIQS =================

template <> chi4_iw_t to_ferm<TRIQS>(const chi4_iw_t &iw_tpl_triqs) {
  auto w1_in  = std::get<0>(iw_tpl_triqs);
  auto w1_out  = std::get<1>(iw_tpl_triqs);
  auto w2_in = std::get<2>(iw_tpl_triqs);

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<TRIQS>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  return chi4_iw_t{w1_in, w1_out, w2_in};
}

// =============== Pomerol Notation =================

template <> chi4_iw_t to_ferm<POMEROL>(const chi4_iw_t &iw_tpl_pom) {
  auto w1_in  = std::get<0>(iw_tpl_pom);
  auto w2_in  = std::get<1>(iw_tpl_pom);
  auto w2_out = std::get<2>(iw_tpl_pom);

  auto w1_out = w1_in + w2_in - w2_out;

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<POMEROL>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto w2_out = w1_in + w2_in - w1_out;

  return chi4_iw_t{w1_in, w2_in, w2_out};
}

// =============== UNSHIFTED NOTATIONS =================

// -- PP ( Vertex paper, only outgoing legs flipped )
template <> chi4_iw_t to_ferm<PP>(const chi4_iw_t &iw_tpl_pp) {
  auto W     = std::get<0>(iw_tpl_pp);
  auto w_in  = std::get<1>(iw_tpl_pp);
  auto w_out = std::get<2>(iw_tpl_pp);

  auto w1_in  = w_in;
  auto w2_in  = W - w_in;
  auto w1_out = w_out;

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<PP>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w1_in + w2_in;
  auto w_in  = w1_in;
  auto w_out = w1_out;

  return chi4_iw_t{W, w_in, w_out};
}

// -- PH ( Vertex paper )
template <> chi4_iw_t to_ferm<PH>(const chi4_iw_t &iw_tpl_ph) {
  auto W     = std::get<0>(iw_tpl_ph);
  auto w_in  = std::get<1>(iw_tpl_ph);
  auto w_out = std::get<2>(iw_tpl_ph);

  auto w1_in  = w_in;
  auto w2_in  = w_out + W;
  auto w1_out = w_in + W;

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<PH>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w1_out - w1_in;
  auto w_in  = w1_in;
  auto w_out = w2_in - (w1_out - w1_in);

  return chi4_iw_t{W, w_in, w_out};
}

// -- XPH ( Vertex paper )
template <> chi4_iw_t to_ferm<XPH>(const chi4_iw_t &iw_tpl_xph) {
  auto W     = std::get<0>(iw_tpl_xph);
  auto w_in  = std::get<1>(iw_tpl_xph);
  auto w_out = std::get<2>(iw_tpl_xph);

  auto w1_in  = w_in;
  auto w2_in  = w_out + W;
  auto w1_out = w_out;

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<XPH>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w2_in - w1_out;
  auto w_in  = w1_in;
  auto w_out = w2_in - (w2_in - w1_out);

  return chi4_iw_t{W, w_in, w_out};
}

// =============== SHIFTED NOTATIONS =================

inline matsubara_freq div2_ceil(matsubara_freq const &W) { return {(W.n - 1000000) / 2 + 500000, W.beta, W.statistic}; }

inline matsubara_freq div2_floor(matsubara_freq const &W) { return {(W.n + 1000000) / 2 - 500000, W.beta, W.statistic}; }

// -- PP Shifted
template <> chi4_iw_t to_ferm<PP_S>(const chi4_iw_t &iw_tpl_pp) {
  auto W     = std::get<0>(iw_tpl_pp);
  auto w_in  = std::get<1>(iw_tpl_pp);
  auto w_out = std::get<2>(iw_tpl_pp);

  auto w1_in  = w_in + div2_ceil(W);
  auto w2_in  = div2_floor(W) - w_in; // care for shift
  auto w1_out = w_out + div2_ceil(W);

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<PP_S>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w1_in + w2_in;
  auto w_in  = w1_in - div2_ceil(w1_in + w2_in);
  auto w_out = w1_out - div2_ceil(w1_in + w2_in);

  return chi4_iw_t{W, w_in, w_out};
}

// -- PH Shifted
template <> chi4_iw_t to_ferm<PH_S>(const chi4_iw_t &iw_tpl_ph) {
  auto W     = std::get<0>(iw_tpl_ph);
  auto w_in  = std::get<1>(iw_tpl_ph);
  auto w_out = std::get<2>(iw_tpl_ph);

  auto w1_in  = w_in - div2_floor(W);
  auto w2_in  = w_out + div2_ceil(W);
  auto w1_out = w_in + div2_ceil(W);

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<PH_S>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w1_out - w1_in;
  auto w_in  = w1_in + div2_floor(w1_out - w1_in);
  auto w_out = w2_in - div2_ceil(w1_out - w1_in);

  return chi4_iw_t{W, w_in, w_out};
}

// -- XPH Shifted
template <> chi4_iw_t to_ferm<XPH_S>(const chi4_iw_t &iw_tpl_xph) {
  auto W     = std::get<0>(iw_tpl_xph);
  auto w_in  = std::get<1>(iw_tpl_xph);
  auto w_out = std::get<2>(iw_tpl_xph);

  auto w1_in  = w_in - div2_floor(W);
  auto w2_in  = w_out + div2_ceil(W);
  auto w1_out = w_out - div2_ceil(W);

  return chi4_iw_t{w1_in, w2_in, w1_out};
}

template <> chi4_iw_t from_ferm<XPH_S>(const chi4_iw_t &iw_tpl_ferm) {
  auto w1_in  = std::get<0>(iw_tpl_ferm);
  auto w2_in  = std::get<1>(iw_tpl_ferm);
  auto w1_out = std::get<2>(iw_tpl_ferm);

  auto W     = w2_in - w1_out;
  auto w_in  = w1_in + div2_floor(w2_in - w1_out);
  auto w_out = w2_in - div2_ceil(w2_in - w1_out);

  return chi4_iw_t{W, w_in, w_out};
}
