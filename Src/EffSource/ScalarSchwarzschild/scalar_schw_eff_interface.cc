#include <complex>
#include <cassert>
#include <iostream>
#include <vector>
#include "EffectiveSource-acceleration.h"

namespace source_interface {
  static std::vector<EffectiveSource*> effsource;

  extern "C" {
    void init_source ( const int* nmodes, const int l[],
		       const int m[], const double* M ) {
      assert(effsource.empty());
      for (int i=0; i<*nmodes; i++) {
        effsource.push_back(new EffectiveSource(l[i], m[i], *M));
      }
    }

/*    void init_source_q ( const int* nmodes, const int l[],
		       const int m[], const long double* M ) {
      assert(effsource.empty());
      double Md = *M;
      for (int i=0; i<*nmodes; i++) {
        effsource.push_back(new EffectiveSource(l[i], m[i], Md));
      }
    } */

    void set_window ( const double* r1, const double* w1, const double* q1,
		      const double* s1, const double* r2, const double* w2,
		      const double* q2, const double* s2, const int* nmodes ) {
      /* As I only test for equality, converting the unsigned size to integer
         should be safe */
      assert(static_cast<int>(effsource.size())==*nmodes);
      for (auto& x : effsource) 
	(*x).set_window ( *r1, *w1, *q1, *s1, *r2, *w2, *q2, *s2 );
    }

/*    void set_window_q ( const long double* r1, const long double* w1,
		        const long double* q1, const long double* s1,
			const long double* r2, const long double* w2,
		        const long double* q2, const long double* s2,
			const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double r1d = *r1;
      double w1d = *w1;
      double q1d = *q1;
      double s1d = *s1;
      double r2d = *r2;
      double w2d = *w2;
      double q2d = *q2;
      double s2d = *s2;
      for (auto& x : effsource) 
	(*x).set_window ( r1d, w1d, q1d, s1d, r2d, w2d, q2d, s2d );
    } */

    void calc_window ( const int* n, const double r[],
		       double Win[], double dWin[], double d2Win[] ) {
      assert(!effsource.empty());
      (*effsource.at(0)).calc_window ( *n, r, Win, dWin, d2Win );
    }

/*    void calc_window_q ( const int* n, const long double r[],
		       double long Win[], double long dWin[],
		       double long d2Win[] ) {
      assert(!effsource.empty());
      double rd[*n];
      double Wind[*n];
      double dWind[*n];
      double d2Wind[*n];
      for (int i=0; i++; i<*n) {
	rd[i] = r[i];
      }
      (*effsource.at(0)).calc_window ( *n, rd, Wind, dWind, d2Wind );
      for (int i=0; i++; i<*n) {
        Win[i] = Wind[i];
        dWin[i] = dWind[i];
        d2Win[i] = d2Wind[i];
      }
    } */
     
    void set_time_window ( const double* T, const double* dT_dt,
		           const double* d2T_dt2, 
			   const int* do_smooth_after_lmax, 
			   const int* nmodes ) {
      assert(static_cast<int>(effsource.size())==*nmodes);
      for (auto& x : effsource)
        if ( (*x).get_l() <= *do_smooth_after_lmax ) {
          (*x).set_time_window ( 1.0, 0.0, 0.0 );
	} else {
          (*x).set_time_window ( *T, *dT_dt, *d2T_dt2 );
        }
    }

/*    void set_time_window_q ( const long double* T, const long double* dT_dt,
		           const long double* d2T_dt2, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double Td = *T;
      double dT_dtd = *dT_dt;
      double d2T_dt2d = *d2T_dt2;
      for (auto& x : effsource)
        (*x).set_time_window ( Td, dT_dtd, d2T_dt2d );
    } */

    void set_particle ( const double* r, const double* phi,
		        const double* ur, const double* En,
			const double* Lz, const double* ar,
			const double* aphi, const double* dardt,
			const double* daphidt, const double* d2ardt2,
			const double* d2aphidt2, const int* nmodes ) {
      assert(static_cast<int>(effsource.size())==*nmodes);
#pragma omp parallel for shared(r,phi,ur,En,Lz,ar,aphi,dardt,daphidt,d2ardt2,d2aphidt2,effsource)
/*      for (auto& x : effsource) 
	(*x).set_particle ( *r, *phi, *ur, *En, *Lz, *ar, *aphi, *dardt,
                            *daphidt, *d2ardt2, *d2aphidt2); */
        for (int i=0; i<*nmodes;i++) {
          (*effsource.at(i)).set_particle ( *r, *phi, *ur, *En, *Lz, *ar, *aphi,
			                 *dardt, *daphidt, *d2ardt2, *d2aphidt2);
        }
    }

/*    void set_particle_q ( const long double* r, const long double* phi,
		        const long double* ur, const long double* En,
		        const long double* Lz, const long double* ar,
			const long double* aphi, const long double* dardt,
                        const long double* daphidt, const long double* d2ardt2,
                        const long double* d2aphidt2, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double rd = *r;
      double phid = *phi;
      double urd = *ur;
      double End = *En;
      double Lzd = *Lz;
      double ard = *ar;
      double aphid = *aphi;
      double dardtd = *dardt;
      double daphidtd = *daphidt;
      double d2ardt2d = *d2ardt2;
      double d2aphidt2d = *d2aphidt2;
      for (auto& x : effsource) 
	(*x).set_particle ( rd, phid, urd, End, Lzd, ard, aphid, dardtd,
                            daphidtd, d2ardt2d, d2aphidt2d );
    } */

/*    void set_particle ( const double* p, const double* e,
		        const double* chi, const double* phi,
			const double* r, const double* ur,
			const int* use_osc, const double* ar,
			const double* aphi, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      for (auto& x : effsource) 
	(*x).set_particle ( *p, *e, *chi, *phi, *r, *ur, *use_osc, *ar, *aphi);
    }

    void set_particle_q ( const long double* p, const long double* e,
		        const long double* chi, const long double* phi,
		        const long double* r, const long double* ur,
                        const int* use_osc, const long double* ar,
			const long double* aphi, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double pd = *p;
      double ed = *e;
      double chid = *chi;
      double phid = *phi;
      double rd = *r;
      double urd = *ur;
      double ard = *ar;
      double aphid = *aphi;
      for (auto& x : effsource) 
	(*x).set_particle ( pd, ed, chid, phid, rd, urd, *use_osc, ard, aphid );
    } */

/*    void set_particle ( const double* p, const double* e,
		        const double* chi, const double* phi,
			const double* r, const double* ur,
			const int* use_osc, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      for (auto& x : effsource) 
	(*x).set_particle ( *p, *e, *chi, *phi, *r, *ur, *use_osc );
    }

    void set_particle_q ( const long double* p, const long double* e,
		        const long double* chi, const long double* phi,
		        const long double* r, const long double* ur,
                        const int* use_osc, const int* nmodes ) {
      assert(effsource.size()==*nmodes);
      double pd = *p;
      double ed = *e;
      double chid = *chi;
      double phid = *phi;
      double rd = *r;
      double urd = *ur;
      for (auto& x : effsource) 
	(*x).set_particle ( pd, ed, chid, phid, rd, urd, *use_osc );
    } */

    void eval_source (const int* mode,  const double* r,
		      double* sre, double* sim ) {
      /* I first check that the mode number is zero or positive */
      assert((*mode)>=0);
      /* I then convert it to size_t and compare it with the size of the
         vector of effective sources. Hopefully, in the event of a compiler
         where the size of a vector is not of type size_t, the compiler will
         complain. */
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> src = (*(effsource.at(*mode)))(*r);
      *sre = std::real(src);
      *sim = std::imag(src);
    }

/*    void eval_source_q (const int* mode,  const long double* r,
		      long double* sre, long double* sim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> src = (*(effsource.at(*mode)))(rd);
      *sre = std::real(src);
      *sim = std::imag(src);
    } */

    void eval_source_all (const int* mode,  const int* n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], double sre[], double sim[] ) {
/*      printf("mode = %d\n", *mode);
      printf("n = %d\n", *n); */
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));

      (*(effsource.at(*mode)))(*n, r, Win, dWin, d2Win, sre, sim);
    }

/*    void eval_source_all_q (const int* mode,  const int* n, 
		            const long double r[], const long double Win[],
			    const long double dWin[], const long double d2Win[],
			    long double sre[], long double sim[] ) {
      assert(effsource.size()>=*mode);
      double rd[*n];
      double Wind[*n];
      double dWind[*n];
      double d2Wind[*n];
      double sred[*n];
      double simd[*n];
      for (int i=0; i++; i<*n) {
	rd[i] = r[i];
	Wind[i] = Win[i];
	dWind[i] = dWin[i];
	d2Wind[i] = d2Win[i];
      }
      (*(effsource.at(*mode)))(*n, rd, Wind, dWind, d2Wind, sred, simd);
      for (int i=0; i++; i<*n) {
	sre[i] = sred[i];
	sim[i] = simd[i];
      }
    } */

    void clean_source () {
      /* Use integers of type size_t in order to avoid type conversions */
      size_t s{effsource.size()};
      for (size_t i=0; i<s; i++) {
        assert(effsource.at(s-1-i)!=nullptr);
        delete effsource.at(s-1-i);
        effsource.pop_back(); 
      }
      assert(effsource.size()==0);
    }

    void Phi ( const int* mode, const double* r,
               double* phire, double* phiim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> phi{(*(effsource.at(*mode))).Phi(*r)};
      *phire = std::real(phi);
      *phiim = std::imag(phi);
    }

/*    void Phi_q ( const int* mode, const long double* r,
               long double* phire, long double* phiim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> phi{(*(effsource.at(*mode))).Phi(rd)};
      *phire = std::real(phi);
      *phiim = std::imag(phi);
    } */

    void dPhi_dr ( const int* mode, const double* r, const int* sign,
                   double* dphidrre, double* dphidrim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> dphidr{(*(effsource.at(*mode))).dPhi_dr(*r, *sign)};
      *dphidrre = std::real(dphidr);
      *dphidrim = std::imag(dphidr);
    }

/*    void dPhi_dr_q ( const int* mode, const long double* r,
                   long double* dphidrre, long double* dphidrim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> dphidr{(*(effsource.at(*mode))).dPhi_dr(rd)};
      *dphidrre = std::real(dphidr);
      *dphidrim = std::imag(dphidr);
    } */

    void dPhi_dt ( const int* mode, const double* r, const int* sign,
                   double* dphidtre, double* dphidtim ) {
      assert((*mode)>=0);
      assert(effsource.size()>=static_cast<size_t>(*mode));
      std::complex<double> dphidt{(*(effsource.at(*mode))).dPhi_dt(*r, *sign)};
      *dphidtre = std::real(dphidt);
      *dphidtim = std::imag(dphidt);
    }

/*    void dPhi_dt_q ( const int* mode, const long double* r,
                   long double* dphidtre, long double* dphidtim ) {
      assert(effsource.size()>=*mode);
      double rd = *r;
      std::complex<double> dphidt{(*(effsource.at(*mode))).dPhi_dt(rd)};
      *dphidtre = std::real(dphidt);
      *dphidtim = std::imag(dphidt);
    } */
  }
}
