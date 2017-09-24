
/* auxvar.c
 * 
 * Copyright (C) 2006 - 2007 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */

/*!
  
  \file auxvar.c

  \brief Functions to define auxiliary variables for the Poisson and Binomial observational model to
  maintain the Gaussian full conditionals.

  The functions in this file can be used to define auxiliary variables for the Poisson and Binomial observational model to
  maintain the Gaussian full conditionals, so that the MCMC-sampling is simplified at the cost of (slightly) more difficuly
  posterior 

  For Poisson observations, the observational model is assumed to be 
  \f[ y_i | x_i \;\sim\; \mbox{Poisson}\left(E_i \exp(x_i)\right) \f]

  For Binomial observations, the observational model is assumed to be 
  \f[ y_i | x_i \;\sim\; \mbox{Binomial}\left(n_i, p(x_i)\right) \f]
  where \f[ p(x_i) = \mbox{logit}^{-1}(x_i) = \exp(x_i)/(1+\exp(x_i)) \f]

  The typical usage of these routines, are as follows.
  - Call either \c GMRFLib_aux_setup_poisson_all() or \c GMRFLib_aux_setup_binomial_all() to setup the auxiliary variables
  - Call \c GMRFLib_aux_init_all() to inititalise the auxiliary variables

  Then for each iteration, 
  - Call \c GMRFLib_aux_update_all() to update all the auxiliary variables
  - Call \c GMRFLib_aux_gauss_approx_all() to obtain the equivalent Gaussian `observations' when conditioning on the auxiliary
  variables

  Finally,
  - Call \c GMRFLib_aux_free_all() to free the auxiliary variables workspace

  Here is a simple example how to use auxiliary variables with Poisson observations:
  \verbinclude example-doxygen-auxvar.txt

*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: auxvar.c,v 1.8 2009/08/26 06:12:47 hrue Exp $ */

#include <math.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

int GMRFLib_mixture_lgamma_testing__intern(void)
{
	/*
	 * This routine was initially written by Rudolf Fruehwirth (Institut fuer Hochenergiephysik Oesterreichische Akademie der
	 * Wissenschaften, Nikolsdorfer Gasse 18, 1050 Wien), and later adapted to GMRFLib by H.Rue 
	 */

	GMRFLib_lgamma_mixture_tp *mycoeffs = NULL;
	double n;

	n = 9999;
	int nc, ic;

	printf("This is a test\n");
	while (n > 0) {
		printf("Enter n:\n");
		scanf("%lf", &n);
		printf("n=%lf\n", n);
		GMRFLib_mixture_lgamma(&mycoeffs, n);
		nc = mycoeffs->ncomp;
		printf("nc=%d\n", nc);
		printf("Weights: ");
		for (ic = 0; ic < nc; ic++)
			printf("%E%s", mycoeffs->weights[ic], "  ");
		printf("\n");
		printf("Means: ");
		for (ic = 0; ic < nc; ic++)
			printf("%E%s", mycoeffs->means[ic], "  ");
		printf("\n");
		printf("Varis: ");
		for (ic = 0; ic < nc; ic++)
			printf("%E%s", mycoeffs->variances[ic], "  ");
		printf("\n");
		Free(mycoeffs);
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_mixture_lgamma(GMRFLib_lgamma_mixture_tp ** mixture, double n)
{
	/*
	 * This routine was initially written by Rudolf Fruehwirth (Institut fuer Hochenergiephysik Oesterreichische Akademie der
	 * Wissenschaften, Nikolsdorfer Gasse 18, 1050 Wien), and later adapted to GMRFLib by H.Rue
	 * 
	 * This routine approximate the log-gamma distribution, with density
	 * 
	 * exp(-n*x - exp(-x)) / Gamma(n)
	 * 
	 * by a mixture of normals
	 * 
	 * \sum_{i=1}^N weights[i] * N(x; means[i], variances[i])
	 * 
	 * with N components.
	 * 
	 */

	const int nrange = 7;
	const int range[7][2] = {
		{1, 4},					       /* 10 components */
		{5, 19},				       /* 9 components */
		{20, 49},				       /* 4 components */
		{50, 439},				       /* 3 components */
		{440, 1599},				       /* 2 components */
		{1600, 10000},				       /* 2 components */
		{10000, 30000}				       /* 2 components */
	};
	const int ncomp[] = { 10, 9, 4, 3, 2, 2, 2 };
	const double P10[4][10] = {
		{0.00396984425, 0.0396244597, 0.16776747, 0.147036501, 0.125306271,
		 0.1014852, 0.103758531, 0.115972617, 0.107065659, 0.0880134482},
		{0.00396860865, 0.039627049, 0.16777003, 0.147037759, 0.125304523,
		 0.101481714, 0.103759705, 0.115973128, 0.107065554, 0.0880119305},
		{0.00396904734, 0.0396280642, 0.167770514, 0.14703607, 0.12530043,
		 0.10148242, 0.103759287, 0.115974323, 0.107066971, 0.0880128738},
		{0.00396883344, 0.039627504, 0.167771274, 0.147036659, 0.125301189,
		 0.101481755, 0.103760036, 0.115974339, 0.107065718, 0.0880126919}
	};
	const double M10[4][10] = {
		{3.55887454, 2.11415904, 0.968124631, 0.51537638, 0.145465449,
		 -0.145346445, -0.416660312, -0.689002855, -0.974965634, -1.27310004},
		{2.78807754, 1.84979328, 0.94844169, 0.577673108, 0.223449219,
		 -0.0831666379, -0.387174155, -0.69613969, -1.01843553, -1.34112844},
		{2.43454312, 1.7315327, 0.942157786, 0.60208557, 0.251664821,
		 -0.0644746918, -0.379817508, -0.696781518, -1.0293035, -1.35705784},
		{2.30484474, 1.71231656, 0.952907078, 0.601128034, 0.252368847,
		 -0.059032783, -0.375605704, -0.699071542, -1.03734211, -1.3609072}
	};
	const double V10[4][10] = {
		{2.62603032, 1.21263644, 0.66586521, 0.256650604, 0.120071142,
		 0.0649909219, 0.0473513798, 0.046639443, 0.0576541602, 0.0888536903},
		{2.39619753, 1.16995764, 0.688870128, 0.307084756, 0.155644328,
		 0.0899360571, 0.0707828448, 0.0751755614, 0.0990773728, 0.15471843},
		{2.16215586, 1.11062998, 0.682294453, 0.324750601, 0.173204837,
		 0.108063698, 0.0917073596, 0.100257256, 0.131371692, 0.200024832},
		{1.92939547, 1.00671896, 0.638983371, 0.322852776, 0.18445103,
		 0.122217472, 0.106400052, 0.116936918, 0.154113316, 0.233525098}
	};
	const double P9[15][9] = {
		{0.0435820277, 0.167794347, 0.147040722, 0.125310654, 0.10147112,
		 0.10376347, 0.115973878, 0.107056197, 0.0880075845},
		{0.0435817033, 0.167795985, 0.1470426, 0.125311016, 0.101470666,
		 0.103763084, 0.115972864, 0.107055437, 0.0880066471},
		{0.0435798639, 0.167797087, 0.147042073, 0.125313291, 0.101470979,
		 0.103761847, 0.115973234, 0.107054351, 0.0880072753},
		{0.043578895, 0.167797426, 0.147041988, 0.125313875, 0.101470922,
		 0.103761581, 0.115973137, 0.107054001, 0.0880081751},
		{0.0435786725, 0.167797743, 0.1470428, 0.125313553, 0.101470946,
		 0.103761391, 0.115973188, 0.10705364, 0.0880080663},
		{0.0435779307, 0.167797788, 0.147042734, 0.125314068, 0.101471449,
		 0.10376142, 0.115973187, 0.107053473, 0.0880079505},
		{0.043576761, 0.167801375, 0.147042624, 0.125314075, 0.101470546,
		 0.103761069, 0.115973226, 0.107051966, 0.0880083593},
		{0.0435771819, 0.167801103, 0.147042441, 0.125313864, 0.101470305,
		 0.103761519, 0.11597319, 0.107052417, 0.0880079809},
		{0.0435778469, 0.167800486, 0.147041951, 0.125313914, 0.101470076,
		 0.103761707, 0.115973611, 0.107052756, 0.0880076518},
		{0.0435786417, 0.16779926, 0.147042119, 0.125313391, 0.101470554,
		 0.103762378, 0.115973792, 0.107052537, 0.0880073289},
		{0.043581505, 0.167797871, 0.147043608, 0.125312521, 0.101469081,
		 0.103762173, 0.115973414, 0.107054363, 0.0880054639},
		{0.0435811435, 0.167798952, 0.147043687, 0.125312616, 0.101468918,
		 0.103762052, 0.115973417, 0.107053968, 0.0880052462},
		{0.0435812603, 0.167798873, 0.147044518, 0.125312321, 0.101468879,
		 0.103761729, 0.115972692, 0.107054049, 0.0880056789},
		{0.0435808733, 0.167799002, 0.147044529, 0.125312675, 0.101468951,
		 0.103761472, 0.115972643, 0.107053883, 0.0880059719},
		{0.0435807283, 0.167799231, 0.14704464, 0.12531292, 0.101468814,
		 0.103761275, 0.115972628, 0.107053662, 0.088006103}
	};
	const double M9[15][9] = {
		{1.31113348, 0.963928895, 0.659198795, 0.240742429, -0.108844644,
		 -0.252087404, -0.6546691, -1.04146524, -1.37874376},
		{1.25919247, 0.957217299, 0.66710982, 0.251658342, -0.125234491,
		 -0.240137829, -0.64912733, -1.03921002, -1.37439461},
		{1.21602216, 0.94778507, 0.671484869, 0.265435387, -0.104709908,
		 -0.24708343, -0.653441223, -1.04076324, -1.36988994},
		{1.18027937, 0.939725546, 0.67760436, 0.293497817, -0.110879079,
		 -0.257696481, -0.655613756, -1.0406543, -1.36465528},
		{1.14996911, 0.934206664, 0.686267712, 0.311595579, -0.112948479,
		 -0.274222612, -0.653808807, -1.04092104, -1.35962481},
		{1.12841748, 0.932206841, 0.69102714, 0.319038554, -0.109581301,
		 -0.302963892, -0.641448217, -1.03858769, -1.35274157},
		{1.10451126, 0.925180162, 0.689947194, 0.309587296, -0.123979787,
		 -0.239246368, -0.658582798, -1.03932069, -1.347407},
		{1.08624068, 0.918081034, 0.697616213, 0.330655882, -0.106424319,
		 -0.290644969, -0.644517493, -1.04099153, -1.34370607},
		{1.0671125, 0.915784215, 0.70024231, 0.330800476, -0.125598534,
		 -0.244656951, -0.661886313, -1.04447342, -1.33948264},
		{1.05721516, 0.918099637, 0.698999193, 0.325014717, -0.153165358,
		 -0.225909041, -0.659788653, -1.03711782, -1.33064663},
		{1.02150943, 0.896206397, 0.702224917, 0.344137939, -0.119895501,
		 -0.256590721, -0.641185291, -1.03810889, -1.32943558},
		{1.02508782, 0.902555642, 0.699256309, 0.336391119, -0.121902141,
		 -0.242730179, -0.6538063, -1.0385784, -1.32415888},
		{0.997274184, 0.88197491, 0.696155279, 0.3460138, -0.128981232,
		 -0.227346713, -0.630322077, -1.03647508, -1.32316505},
		{0.995086849, 0.891409674, 0.70171109, 0.341992158, -0.127906113,
		 -0.245952673, -0.638792902, -1.03392281, -1.31486719},
		{0.997741814, 0.892112396, 0.698155553, 0.337731787, -0.122630195,
		 -0.240878604, -0.651951415, -1.02898878, -1.3062535}
	};
	const double V9[15][9] = {
		{1.5732832, 0.745075965, 0.340530976, 0.206325108, 0.206977107,
		 0.133034557, 0.123981078, 0.155417698, 0.247661591},
		{1.52550277, 0.745216293, 0.347702459, 0.213195645, 0.220928839,
		 0.147502243, 0.139478204, 0.17271313, 0.269719569},
		{1.48970429, 0.74910777, 0.35810967, 0.221696291, 0.216470192,
		 0.155837875, 0.148481868, 0.185394632, 0.28822907},
		{1.46105103, 0.752441091, 0.365198621, 0.220104509, 0.199190433,
		 0.167708126, 0.15761138, 0.197076001, 0.304425302},
		{1.43764551, 0.754190306, 0.367534375, 0.215776065, 0.185257157,
		 0.180142183, 0.165402413, 0.206954388, 0.318591695},
		{1.41468216, 0.75198881, 0.368357589, 0.215271168, 0.178178434,
		 0.198636491, 0.176790288, 0.218155881, 0.332156859},
		{1.39851898, 0.755429842, 0.377058085, 0.229287048, 0.214645547,
		 0.18489307, 0.178139004, 0.226237823, 0.343708183},
		{1.38111403, 0.759024378, 0.379809227, 0.222659694, 0.185443843,
		 0.206181273, 0.184773494, 0.231840962, 0.353714302},
		{1.36922993, 0.759197249, 0.381687395, 0.225704876, 0.199623554,
		 0.195711194, 0.18270427, 0.236837387, 0.363050264},
		{1.35398708, 0.753650144, 0.381381699, 0.231057971, 0.208319112,
		 0.210190241, 0.194957855, 0.249236388, 0.373774124},
		{1.35652837, 0.774748407, 0.400413698, 0.238592235, 0.199467639,
		 0.230239828, 0.19924794, 0.251600772, 0.380054821},
		{1.33546695, 0.763749521, 0.396745563, 0.241905327, 0.212176877,
		 0.218950701, 0.201882762, 0.257807637, 0.388524892},
		{1.33575722, 0.781739895, 0.417799104, 0.256728889, 0.211230256,
		 0.254750255, 0.208700024, 0.26087813, 0.393822372},
		{1.3227643, 0.771070524, 0.406631212, 0.249617029, 0.210958958,
		 0.249496089, 0.214362668, 0.270024593, 0.402615529},
		{1.30630549, 0.765952536, 0.407914566, 0.255018833, 0.226289944,
		 0.236887588, 0.221124118, 0.280039124, 0.411219814}
	};
	/*
	 * n from 20 to 49 
	 */
#define size_Coeff_p 4
#define size_Coeff_m 5
#define size_Coeff_v 5
	const double Coeff_p3[size_Coeff_p][4] = {
		{-5.644536495326e-009, 7.299190941772e-009, -1.788056445701e-008, 9.259794020097e-009},
		{-1.266992312621e-006, 1.387196986613e-006, -2.391966642312e-006, 1.224613603301e-006},
		{1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
		{4.730022618640e+000, 3.672627139064e+000, 4.871566292199e+000, 3.215154075256e+000}
	};
#define size_Coeff_m3 5
	const double Coeff_m3[size_Coeff_m][4] = {
		{4.552797222246e-005, 2.284729919322e-005, -3.900177124794e-005, -2.486737015928e-005},
		{4.638009105861e-002, -1.095058888700e-002, 4.731686443506e-002, 2.978371498898e-002},
		{1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
		{9.627160143020e-002, -1.690501643196e-002, -5.610095109269e-001, -4.643825308040e-002},
		{1.143772956136e+000, 1.944583776810e+000, -6.033854619021e+000, -1.105498133467e+000}
	};
#define size_Coeff_v3 5
	const double Coeff_v3[size_Coeff_v][4] = {
		{-2.191015160635e-005, 7.060864706965e-005, 1.823003483481e-004, 1.613752763707e-004},
		{9.939739739229e-002, 1.143203813438e-001, 1.675101633325e-001, 1.943336591437e-001},
		{1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000, 1.000000000000e+000},
		{9.208564449364e-002, 1.548617268518e-001, 2.735383307281e-001, 2.797653349940e-001},
		{7.148740127686e-001, 2.428636911969e+000, 4.861423133312e+000, 3.840341872065e+000}
	};
	/*
	 * n from 50 to 439 
	 */
	const double Coeff_p4[size_Coeff_p][3] = {
		{-5.639545796991280e-010, 2.651836392450035e-010, -2.384482520627535e-011},
		{4.698743002874532e-007, -1.280380156002802e-007, -1.227680572544847e-007},
		{1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
		{4.730482920811330e+000, 2.093982718501769e+000, 3.214956149674574e+000}
	};
	const double Coeff_m4[size_Coeff_m][3] = {
		{-1.653173201148335e-006, -8.298537364426537e-007, -1.431525987300163e-006},
		{1.036578627632170e-002, 5.017456263052972e-003, 8.386323466104712e-003},
		{1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
		{2.349390607953272e-002, 5.123168011502721e-002, -1.841057020139425e-002},
		{1.432904956685477e+000, 6.453910704667408e+000, -1.410602407670769e+000}
	};
	const double Coeff_v4[size_Coeff_v][3] = {
		{-2.726183914412441e-007, 1.118379212729684e-006, 2.197737873275589e-006},
		{2.788507874891710e-002, 2.433214514397419e-002, 3.186581505796005e-002},
		{1.000000000000000e+000, 1.000000000000000e+000, 1.000000000000000e+000},
		{2.777086294607445e-002, 2.778340896223197e-002, 3.808382220884354e-002},
		{8.369406298984288e-001, 1.489387981224663e+000, 1.958805931276004e+000}
	};
	/*
	 * n from 440 to 1599 
	 */
	const double Coeff_p5[size_Coeff_p][2] = {
		{1.034981541036597e-010, -2.291586556531707e-010},
		{-2.445177000398938e-007, 5.414543692806514e-007},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{1.451229377475864e+000, 3.216167113242079e+000}
	};
	const double Coeff_m5[size_Coeff_m][2] = {
		{-6.578325435644067e-008, -6.292364160498604e-008},
		{1.648723149067166e-003, 1.618047470065775e-003},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{1.594968525045459e-002, -7.091699113800587e-003},
		{5.566082591106806e+000, -2.516741952410371e+000}
	};
	const double Coeff_v5[size_Coeff_v][5] = {
		{-2.802162650788337e-009, 3.776558110733883e-008, 4.051503597935380e-003, 5.022619018941299e-003,
		 1.000000000000000e+000},
		{1.000000000000000e+000, 4.018981069179972e-003, 5.525253413878772e-003, 9.654061278849895e-001,
		 1.450507513327352e+000}
	};
	/*
	 * n from 1600 to 10000 
	 */
	const double Coeff_p6[size_Coeff_p][2] = {
		{-1.586037487404490e-013, 1.291237745205579e-013},
		{3.575996226727867e-009, -2.911316152726367e-009},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{2.228310599179340e+000, 1.814126328168031e+000}
	};
	const double Coeff_m6[size_Coeff_m][2] = {
		{-2.419956255676409e-009, -2.419411092563945e-009},
		{3.245753451748892e-004, 3.245669014250788e-004},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{1.895335618211674e-003, -2.327930564510444e-003},
		{3.388553853864067e+000, -4.162274939236667e+000}
	};
	const double Coeff_v6[size_Coeff_v][2] = {
		{-6.024563976875348e-011, 5.024317053887777e-010},
		{-6.540694956580495e-004, 8.898044793516080e-004},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{-6.582951415419203e-004, 9.246987493760628e-004},
		{1.006399508694657e+000, 1.149073788967684e+000}
	};
	/*
	 * n from 10000 to 30000 
	 */
	const double Coeff_p7[size_Coeff_p][2] = {
		{-1.663426552872397e-014, 1.354267905471566e-014},
		{1.141056828884990e-009, -9.289835742028532e-010},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{2.228285989630589e+000, 1.814142639751731e+000}
	};
	const double Coeff_m7[size_Coeff_m][2] = {
		{-8.929405559006038e-011, -8.931137480031157e-011},
		{6.319814700961324e-005, 6.320244393309693e-005},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{4.785131212048377e-004, -5.877524860249395e-004},
		{4.271922830906078e+000, -5.247218808668549e+000}
	};
	const double Coeff_v7[size_Coeff_v][2] = {
		{-1.418731402291282e-012, 1.576782807097003e-011},
		{-5.512224505288543e-006, 1.914006058179041e-004},
		{1.000000000000000e+000, 1.000000000000000e+000},
		{-5.638714069888806e-006, 1.959753272178233e-004},
		{1.006201804172733e+000, 1.087101027065273e+000}
	};

	*mixture = Calloc(1, GMRFLib_lgamma_mixture_tp);
	(*mixture)->weights[0] = 0;
	(*mixture)->means[0] = 0;
	(*mixture)->variances[0] = 0;
	(*mixture)->ncomp = 0;

	/*
	 * Check input 
	 */
	GMRFLib_ASSERT(n > 0, GMRFLib_EPARAMETER);
	int nr = (int) n;

	if (n < range[2][0] && nr != n) {
		GMRFLib_ASSERT(n < range[2][0] && (int) n != n, GMRFLib_EPARAMETER);
		return GMRFLib_EPARAMETER;
	}

	/*
	 * Mean and standard deviation of log-gamma 
	 */
	double mu = -gsl_sf_psi((double) n);
	double sigma = sqrt(gsl_sf_psi_1((double) n));

	/*
	 * Single component 
	 */
	if (n > range[6][1]) {
		(*mixture)->weights[0] = 1.0;
		(*mixture)->means[0] = mu;
		(*mixture)->variances[0] = SQR(sigma);
		(*mixture)->ncomp = 1;

		return GMRFLib_SUCCESS;
	}

	/*
	 * Find appropriate range 
	 */
	int rtake = 0, ir, ic, nc, jc;

	for (ir = 0; ir < nrange; ir++) {
		if (range[ir][0] <= n && n < range[ir][1] + 1) {
			rtake = ir;
		}
	}

	/*
	 * Set pointers to coefficients 
	 */
	const double *P = NULL, *M = NULL, *V = NULL, *ptr = NULL;

	nr = nr - range[rtake][0];
	switch (rtake) {
	case 0:
		P = &P10[nr][0];
		M = &M10[nr][0];
		V = &V10[nr][0];
		break;
	case 1:
		P = &P9[nr][0];
		M = &M9[nr][0];
		V = &V9[nr][0];
		break;
	case 2:
		P = &Coeff_p3[0][0];
		M = &Coeff_m3[0][0];
		V = &Coeff_v3[0][0];
		break;
	case 3:
		P = &Coeff_p4[0][0];
		M = &Coeff_m4[0][0];
		V = &Coeff_v4[0][0];
		break;
	case 4:
		P = &Coeff_p5[0][0];
		M = &Coeff_m5[0][0];
		V = &Coeff_v5[0][0];
		break;
	case 5:
		P = &Coeff_p6[0][0];
		M = &Coeff_m6[0][0];
		V = &Coeff_v6[0][0];
		break;
	case 6:
		P = &Coeff_p7[0][0];
		M = &Coeff_m7[0][0];
		V = &Coeff_v7[0][0];
		break;
	default:
		abort();
	}

	/*
	 * Store number of components 
	 */
	nc = ncomp[rtake];
	(*mixture)->ncomp = nc;

	/*
	 * Explicit mixture parameters 
	 */
	if (rtake <= 1) {
		for (ic = 0; ic < nc; ic++) {
			(*mixture)->weights[ic] = *P++;
			(*mixture)->means[ic] = (*M++) * sigma + mu;
			(*mixture)->variances[ic] = (*V++) * SQR(sigma);
		}

		return GMRFLib_SUCCESS;
	}

	/*
	 * Mixture parameters by rational approximation 
	 */
	double numer = 0, sum;

	for (ic = 0; ic < nc; ic++) {
		/*
		 * Weights 
		 */
		ptr = P++;
		for (sum = jc = 0; jc < size_Coeff_p; jc++) {
			sum = sum * n + (*ptr);
			if (*ptr == 1) {
				numer = sum;
				sum = 0;
			}
			ptr = ptr + nc;
		}
		(*mixture)->weights[ic] = numer / sum;

		/*
		 * Means 
		 */
		ptr = M++;
		for (sum = jc = 0; jc < size_Coeff_m; jc++) {
			sum = sum * n + (*ptr);
			if (*ptr == 1) {
				numer = sum;
				sum = 0;
			}
			ptr = ptr + nc;
		}
		(*mixture)->means[ic] = numer / sum * sigma + mu;

		/*
		 * Variances 
		 */
		ptr = V++;
		for (sum = jc = 0; jc < size_Coeff_v; jc++) {
			sum = sum * n + (*ptr);
			if (*ptr == 1) {
				numer = sum;
				sum = 0;
			}
			ptr = ptr + nc;
		}
		(*mixture)->variances[ic] = numer / sum * SQR(sigma);
	}

#undef size_Coeff_p
#undef size_Coeff_m
#undef size_Coeff_v
#undef size_Coeff_m3
#undef size_Coeff_v3

	return GMRFLib_SUCCESS;
}
int GMRFLib_aux_setup_poisson(GMRFLib_aux_tp ** aux, double *y, double *E, double *x)
{
	/*
	 * setup poisson: use ptrs...
	 */

	*aux = Calloc(1, GMRFLib_aux_tp);

	(*aux)->poisson = Calloc(1, GMRFLib_poisson_aux_tp);
	(*aux)->poisson->y = y;
	(*aux)->poisson->E = E;
	(*aux)->poisson->x = x;

	GMRFLib_aux_init(*aux);

	return GMRFLib_SUCCESS;
}
int GMRFLib_aux_setup_binomial(GMRFLib_aux_tp ** aux, double *y, double *nb, double *x)
{
	/*
	 * setup binomial: use ptrs...
	 */

	*aux = Calloc(1, GMRFLib_aux_tp);

	(*aux)->binomial = Calloc(1, GMRFLib_binomial_aux_tp);
	(*aux)->binomial->y = y;
	(*aux)->binomial->nb = nb;
	(*aux)->binomial->x = x;

	GMRFLib_aux_init(*aux);

	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Setup the auxiliary variables for Poisson observations \f$ y_i \sim \mbox{Poisson}(E_i\exp(x_i)) \f$.

  \param[out] auxs A pointer-pointer to the auxiliary variables objects.
  \param[in] n The dimension of the latent variable \c x
  \param[in] d A vector of length \c n of weights for the Poisson observations, similar to \c d in \c GMRFLib_blockupdate() etc.
  \param[in] y A vector of length \c n of observations
  \param[in] E A vector of length \c n of scale varibles
  \param[in] x The latent variable of length \c n

  Note that only the pointers to \c d, \c y and most importantly, \c x, are stored and later used, NOT the values of the vectors.
*/
int GMRFLib_aux_setup_poisson_all(GMRFLib_auxs_tp ** auxs, int n, double *d, double *y, double *E, double *x)
{
	/*
	 * OOPS: the pointers to y,E and x are used, so be aware!!!! 
	 */
	int i;

	*auxs = Calloc(n, GMRFLib_auxs_tp);
	(*auxs)->n = n;
	(*auxs)->d = d;
	(*auxs)->aux = Calloc(n, GMRFLib_aux_tp *);

	for (i = 0; i < n; i++) {
		if (d[i]) {
			GMRFLib_aux_setup_poisson(&((*auxs)->aux[i]), &y[i], &E[i], &x[i]);
		} else {
			(*auxs)->aux[i] = NULL;
		}
	}
	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Setup the auxiliary variables for Binomial observations \f$ y_i \sim \mbox{Binomial}(nb_i, p(x_i)) \f$.

  \param[out] auxs A pointer-pointer to the auxiliary variables objects.
  \param[in] n The dimension of the latent variable \c x
  \param[in] d A vector of length \c n of weights for the Poisson observations, similar to \c d in \c GMRFLib_blockupdate() etc.
  \param[in] y A vector of length \c n of observations
  \param[in] nb A vector of length \c n of number of trials in the Binomial
  \param[in] x The latent variable of length \c n

  Note that only the pointers to \c d, \c y and most importantly, \c x, are stored and later used, NOT the values of the vectors.
*/
int GMRFLib_aux_setup_binomial_all(GMRFLib_auxs_tp ** auxs, int n, double *d, double *y, double *nb, double *x)
{
	/*
	 * OOPS: the pointers to y, nb and x are used, so be aware!!!! 
	 */
	int i;

	*auxs = Calloc(n, GMRFLib_auxs_tp);
	(*auxs)->n = n;
	(*auxs)->d = d;
	(*auxs)->aux = Calloc(n, GMRFLib_aux_tp *);

	for (i = 0; i < n; i++) {
		if (d[i]) {
			GMRFLib_aux_setup_binomial(&((*auxs)->aux[i]), &y[i], &nb[i], &x[i]);
		} else {
			(*auxs)->aux[i] = NULL;
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_init(GMRFLib_aux_tp * aux)
{
	/*
	 * initialise the auxiliary variables in the aux-representation 
	 */

	if (!aux) {
		return GMRFLib_SUCCESS;
	}

	if (aux->poisson) {
		return GMRFLib_aux_init_poisson(aux->poisson);
	}
	if (aux->binomial) {
		return GMRFLib_aux_init_binomial(aux->binomial);
	}

	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Initialise the auxiliary variables defined in \c auxs

  \param[in] auxs A pointer-pointer to the auxiliary variables objects.
  
*/
int GMRFLib_aux_init_all(GMRFLib_auxs_tp * auxs)
{
	if (auxs) {
		int i;

		for (i = 0; i < auxs->n; i++) {
			if (auxs->d[i]) {
				GMRFLib_aux_init(auxs->aux[i]);
			}
		}
	}
	return GMRFLib_SUCCESS;
}
int GMRFLib_aux_gauss_approx(double *b, double *c, GMRFLib_aux_tp * aux)
{
	/*
	 * compute the 'b' and 'c' in the gaussian approximation, ie
	 * 
	 * -1/2 * c * x^2 + b * x 
	 */
	if (!aux) {
		if (b) {
			*b = 0.0;
		}
		if (c) {
			*c = 0.0;
		}
		return GMRFLib_SUCCESS;
	}

	if (aux->poisson) {
		return GMRFLib_aux_gauss_approx_poisson(b, c, aux->poisson);
	}
	if (aux->binomial) {
		return GMRFLib_aux_gauss_approx_binomial(b, c, aux->binomial);
	}

	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Compute the equivalent Gaussian `observations' when conditioning on the auxiliary variables

  When conditioning on the auxiliary variables we obtain from the auxiliary variables the Gaussian density \f[ \pi(x | w)
  \propto \exp\left( -\frac{1}{2} x^T \mbox{diag}(c) x + b^T x \right) \f]
  This function compute the vectors \c b and \c c.

  \param[in] auxs A pointer to the auxiliary variables objects.
  \param[out] b The vector \c b in the Gaussian density (with length \c GMRFLib_auxs_tp::n)
  \param[out] c The vector \c c in the Gaussian density (with length \c GMRFLib_auxs_tp::n)
  
*/
int GMRFLib_aux_gauss_approx_all(double *b, double *c, GMRFLib_auxs_tp * auxs)
{
	if (auxs) {
		int i;

		for (i = 0; i < auxs->n; i++) {
			if (auxs->d[i]) {
				GMRFLib_aux_gauss_approx(&b[i], &c[i], auxs->aux[i]);
				b[i] *= auxs->d[i];
				c[i] *= auxs->d[i];
			} else {
				b[i] = c[i] = 0.0;
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_update(GMRFLib_aux_tp * aux)
{
	/*
	 * update the auxiliary variables in the aux-representation 
	 */

	if (!aux) {
		return GMRFLib_SUCCESS;
	}

	if (aux->poisson) {
		return GMRFLib_aux_update_poisson(aux->poisson);
	}
	if (aux->binomial) {
		return GMRFLib_aux_update_binomial(aux->binomial);
	}

	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Update all the auxiliary variables conditioned on GMRFLib_auxs_tp::x

  \param[in] auxs A pointer to the auxiliary variables objects.
  
*/
int GMRFLib_aux_update_all(GMRFLib_auxs_tp * auxs)
{
	if (auxs) {
		int i;

		for (i = 0; i < auxs->n; i++) {
			if (auxs->d[i]) {
				GMRFLib_aux_update(auxs->aux[i]);
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_free(GMRFLib_aux_tp * aux)
{
	if (!aux) {
		return GMRFLib_SUCCESS;
	}
	if (aux->poisson) {
		Free(aux->poisson->mix1);
		Free(aux->poisson->mix2);
		Free(aux->poisson);
	}
	if (aux->binomial) {
		Free(aux->binomial->mix);
		Free(aux->binomial);
	}
	Free(aux);

	return GMRFLib_SUCCESS;
}

/*!
  
  \brief Free all auxiliary variables 

  \param[in] auxs A pointer to the auxiliary variables objects to be free'd.
  
*/
int GMRFLib_aux_free_all(GMRFLib_auxs_tp * auxs)
{
	if (auxs) {
		int i;

		for (i = 0; i < auxs->n; i++) {
			if (auxs->aux[i]) {
				GMRFLib_aux_free(auxs->aux[i]);
			}
		}
		Free(auxs);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_init_poisson(GMRFLib_poisson_aux_tp * poisson)
{
	/*
	 * init the aux 
	 */

	double xi, lambda;

	GMRFLib_mixture_lgamma(&(poisson->mix1), 1.0);
	poisson->r1 = (int) (GMRFLib_uniform() * poisson->mix1->ncomp);
	if (*(poisson->y) > 0.0) {
		GMRFLib_mixture_lgamma(&(poisson->mix2), *(poisson->y));
		poisson->r2 = (int) (GMRFLib_uniform() * poisson->mix2->ncomp);
	}

	lambda = DMAX(0.1, *(poisson->y));
	xi = gsl_ran_exponential(GMRFLib_rng, 1.0) / lambda;

	/*
	 * update first the inter-arrival times 
	 */
	if (ISZERO(*(poisson->y))) {
		poisson->tau1 = 1.0 + xi;
	} else {
		poisson->tau2 = gsl_ran_beta(GMRFLib_rng, *(poisson->y), 1.0);
		poisson->tau1 = 1.0 - poisson->tau2 + xi;
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_init_binomial(GMRFLib_binomial_aux_tp * binomial)
{
	/*
	 * init the aux 
	 */
	double lambda, u, v;

	GMRFLib_mixture_lgamma(&(binomial->mix), *(binomial->nb));
	binomial->r = (int) (GMRFLib_uniform() * binomial->mix->ncomp);
	lambda = 1.0;

	/*
	 * sample the aggregated utility 
	 */
	u = gsl_ran_gamma(GMRFLib_rng, *(binomial->nb), 1.0);
	if ((int) *(binomial->y) < (int) *(binomial->nb)) {
		v = gsl_ran_gamma(GMRFLib_rng, *(binomial->nb) - *(binomial->y), 1.0);
	} else {
		v = 0.0;
	}
	binomial->y_star = -log(u / (1.0 + lambda) + v / lambda);

	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_update_poisson(GMRFLib_poisson_aux_tp * poisson)
{
	double xi, lambda;

	lambda = *(poisson->E) * exp(*(poisson->x));
	xi = gsl_ran_exponential(GMRFLib_rng, 1.0) / lambda;

	/*
	 * update first the inter-arrival times 
	 */
	if (ISZERO(*(poisson->y))) {
		poisson->tau1 = 1.0 + xi;
	} else {
		poisson->tau2 = gsl_ran_beta(GMRFLib_rng, *(poisson->y), 1.0);
		poisson->tau1 = 1.0 - poisson->tau2 + xi;
	}

	/*
	 * ...then the mixture indicator 
	 */
	GMRFLib_mixture_update_indicator(&(poisson->r1), poisson->tau1, lambda, poisson->mix1);
	if (*(poisson->y) > 0.0) {
		GMRFLib_mixture_update_indicator(&(poisson->r2), poisson->tau2, lambda, poisson->mix2);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_aux_update_binomial(GMRFLib_binomial_aux_tp * binomial)
{
	double lambda = exp(*(binomial->x)), u, v;

	/*
	 * sample the aggregated utility 
	 */
	u = gsl_ran_gamma(GMRFLib_rng, *(binomial->nb), 1.0);
	if ((int) *(binomial->y) < (int) *(binomial->nb)) {
		v = gsl_ran_gamma(GMRFLib_rng, *(binomial->nb) - *(binomial->y), 1.0);
	} else {
		v = 0.0;
	}
	binomial->y_star = -log(u / (1.0 + lambda) + v / lambda);

	/*
	 * ...then the mixture indicator 
	 */
	GMRFLib_mixture_update_indicator(&(binomial->r), exp(-binomial->y_star), lambda, binomial->mix);

	return GMRFLib_SUCCESS;
}
int GMRFLib_mixture_update_indicator(int *r, double tau, double lambda, GMRFLib_lgamma_mixture_tp * mix)
{
	/*
	 * update the mixture indicator. this function use the notation for the Poisson mixture case, but this is easily adapted
	 * to the Binomial case, using tau = exp(-y_star)
	 */

	int i;
	double *p = Calloc(mix->ncomp, double);
	double log_lambda = log(lambda), log_tau = log(tau);

	/*
	 * do this safe: compute in log-scale, normalise, then convert 
	 */
	for (i = 0; i < mix->ncomp; i++) {
		p[i] = -0.5 * SQR(log_tau + log_lambda + mix->means[i]) / mix->variances[i] + log(mix->weights[i] / sqrt(mix->variances[i]));
	}
	GMRFLib_adjust_vector(p, mix->ncomp);
	for (i = 0; i < mix->ncomp; i++) {
		p[i] = exp(p[i]);
	}

	gsl_ran_discrete_t *table = gsl_ran_discrete_preproc((size_t) mix->ncomp, (const double *) p);

	*r = (int) gsl_ran_discrete(GMRFLib_rng, table);

	gsl_ran_discrete_free(table);
	Free(p);

	return GMRFLib_SUCCESS;
}
int GMRFLib_aux_gauss_approx_poisson(double *b, double *c, GMRFLib_poisson_aux_tp * poisson)
{
	double bb, cc, mean, var;

	var = poisson->mix1->variances[poisson->r1];
	mean = -(log(*(poisson->E)) + poisson->mix1->means[poisson->r1] + log(poisson->tau1));

	cc = 1.0 / var;
	bb = mean / var;

	if (*(poisson->y) > 0.0) {
		var = poisson->mix2->variances[poisson->r2];
		mean = -(log(*(poisson->E)) + poisson->mix2->means[poisson->r2] + log(poisson->tau2));
		cc += 1.0 / var;
		bb += mean / var;
	}

	if (b) {
		*b = bb;
	}
	if (c) {
		*c = cc;
	}

	return GMRFLib_SUCCESS;
}
int GMRFLib_aux_gauss_approx_binomial(double *b, double *c, GMRFLib_binomial_aux_tp * binomial)
{
	double bb, cc, mean, var;

	var = binomial->mix->variances[binomial->r];
	mean = binomial->y_star - binomial->mix->means[binomial->r];

	cc = 1.0 / var;
	bb = mean / var;

	if (b) {
		*b = bb;
	}
	if (c) {
		*c = cc;
	}

	return GMRFLib_SUCCESS;
}
