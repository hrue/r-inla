#include <cstddef>
#include <cstring>
#include <set>
#include <map>
#include <sstream>
#include <cmath>
#include "gsl/gsl_sf_legendre.h"

#include "vector.hh"
#include "ioutils.hh"

#define WHEREAMI __FILE__ << "(" << __LINE__ << ")\t"

#ifdef DEBUG
#define _LOG(msg) std::cout << WHEREAMI << msg;
#else
#define _LOG(msg)
#endif

#define M_2_SQRT_PI 3.5449077018110320546


using std::ios;
using std::cout;
using std::cin;
using std::endl;

namespace fmesh {

  int sph_basis_n(int kmax, bool rot_sym)
  {
    if (kmax >= 0) {
      if (rot_sym)
	return (kmax + 1);
      else
	return (kmax + 1) * (kmax + 1);
    } else {
      return 0;
    }
  }

  Matrix<double> spherical_harmonics(Matrix3<double>& S,
				     int max_order,
				     bool rotationally_symmetric)
  {
    Matrix<double> sph(sph_basis_n(max_order,rotationally_symmetric));
    int i, k, m;
    double GSL_res_array[max_order+1];

    if (rotationally_symmetric) {    
      for (i=0; i<S.rows(); i++) {
	gsl_sf_legendre_sphPlm_array(max_order, 0, S[i][2], GSL_res_array);
	for (k = 0; k <= max_order; k++) {
	  sph(i,k) = M_2_SQRT_PI * GSL_res_array[k];
	}
      }
    } else {
        double phi, scaling_sin, scaling_cos;

	int Idxs2[max_order + 1];
	// 0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, 156, 182, 210, 240
	for (k = 0; k <= max_order; k++) {
	  Idxs2[k] = k*(k+1);
	}

	for (i=0; i<S.rows(); i++) {
	  phi = atan2(S[i][1], S[i][0]);
	  
	  gsl_sf_legendre_sphPlm_array(max_order, 0, S[i][2], GSL_res_array);
	  for (k = 0; k <= max_order; k++) {
	    sph(i,Idxs2[k]) = M_2_SQRT_PI * GSL_res_array[k];
	  }
	  for (m = 1; m <= max_order; m++) {
	    scaling_sin = M_2_SQRT_PI * M_SQRT2 * sin(-m * phi);
	    scaling_cos = M_2_SQRT_PI * M_SQRT2 * cos(m * phi);
	    gsl_sf_legendre_sphPlm_array(max_order, m, S[i][2], GSL_res_array);
	    for (k = m; k <= max_order; k++) {
	      sph(i,Idxs2[k] - m) = scaling_sin * GSL_res_array[k - m];
	      sph(i,Idxs2[k] + m) = scaling_cos * GSL_res_array[k - m];
	    }
	  }
	}

   }

    return sph;
  }

} /* namespace fmesh */
