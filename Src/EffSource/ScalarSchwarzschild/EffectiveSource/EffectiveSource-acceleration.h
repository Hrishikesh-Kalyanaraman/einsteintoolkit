/*
 * Copyright 2015 Barry Wardell
 *
 * This file is distributed under the University of Illinois/NCSA
 * Open Source License. See LICENSE file for details.
 *
 */

#ifndef EFFECTIVESOURCE_H
#define EFFECTIVESOURCE_H

#include <complex>
#include <array>

class EffectiveSource
{
public:
  std::complex<double> operator()(const double r);
  void operator()(const int n, const double r[],
		  const double Win[], const double dWin[],
		  const double d2Win[],
		  double sre[], double sim[]);
  std::complex<double> Phi(const double r);
  std::complex<double> dPhi_dr(const double r, const int sign);
  std::complex<double> dPhi_dt(const double r, const int sign);
  EffectiveSource(const int l, const int m, const double M);

  void set_particle(const double r, const double phi, const double ur, const double E0, const double L, const double ar, const double aphi, const double dardt, const double daphidt, const double d2ardt2, const double d2aphidt2);
  void set_window(
    const double r1, const double w1, const double q1, const double s1,
    const double r2, const double w2, const double q2, const double s2);
  void set_time_window(const double T, const double dT_dt, const double d2T_dt2);
  void calc_window(
    const int n, const double r[], double Win[], double dWin[], double d2Win[]);
    int get_l();
    int get_m();


private:
  const double M;
  const int l, m;

  double r0, phi0;
  double r1, w1, q1, s1, r2, w2, q2, s2;
  double T, dT_dt, d2T_dt2;
  std::complex<double> w0, w1p, w1m, w2p, w2m;
  std::array<std::array<double,5>,3> coeffs, dcoeffs_dt, d2coeffs_dt2;
  std::array<std::array<double,4>,3> abscoeffs, dabscoeffs_dt, d2abscoeffs_dt2;
  double rt, rtt, phit, phitt, c, ct, ctt;
};

#endif /* EFFECTIVESOURCE_H */
