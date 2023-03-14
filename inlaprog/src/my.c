
/* my.c
 * 
 * Copyright (C) 2016-2023 Havard Rue
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "inla.h"
#include "my.h"
#include "my-fix.h"
#include "GMRFLib/GMRFLibP.h"
#include "GMRFLib/GMRFLib.h"

int my_file_exists(const char *filename)
{
	struct stat sb;
	return ((stat(filename, &sb) == 0 && S_ISREG(sb.st_mode)) ? INLA_OK : !INLA_OK);
}

int my_dir_exists(const char *dirname)
{
	struct stat sb;
	return ((stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode)) ? INLA_OK : !INLA_OK);
}

int my_setenv(char *str, int prefix)
{
	/*
	 * set a variable in the enviroment; if PREFIX prepend with inla_, so that a=b yields inla_a=b. 
	 */
	char *p = NULL, *var = NULL;
	const int debug = 0;

	if (debug)
		printf("enter my_setenv with [%s]\n", str);

	p = strchr(str, '=');
	if (!p) {
		fprintf(stderr, "*** Error: Argument is void [%s]\n", str);
		exit(EXIT_FAILURE);
	}
	*p = '\0';
#if defined(WINDOWS)
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s=%s", str, p + 1);
	} else {
		GMRFLib_sprintf(&var, "%s=%s", str, p + 1);
	}
	putenv(var);
	if (debug)
		printf("putenv \t%s\n", var);
#else
	if (prefix) {
		GMRFLib_sprintf(&var, "inla_%s", str);
	} else {
		GMRFLib_sprintf(&var, "%s", str);
	}
	setenv(var, p + 1, 1);
	if (debug)
		printf("\tsetenv %s=%s\n", var, p + 1);
#endif
	Free(var);
	return INLA_OK;
}

double my_gsl_sf_lnfact(int x)
{
	static int first = 1;
	static int nmax = 1048576;
	static double *lng = NULL;

	if (first) {
#pragma omp critical (Name_764ffe066cfbd16ba7b1096b9b762ad9b1f8e669)
		if (first) {
			lng = Calloc(nmax, double);
			lng[0] = 0.0;
			for (int i = 1; i < nmax; i++) {
				lng[i] = lng[i - 1] + log((double) i);
			}
			first = 0;
		}
	}

	if (x >= nmax) {
		return gsl_sf_lnfact((unsigned int) x);
	} else {
		return lng[x];
	}

	assert(0 == 1);
	return NAN;
}

double my_gsl_sf_lngamma(double x)
{
	if (round(x) != x) {
		return gsl_sf_lngamma(x);
	} else {
		// x is an int, then use the cached values

		static int first = 1;
		static int nmax = 1048576;
		static double *lng = NULL;

		if (first) {
#pragma omp critical (Name_72a7f789baa1bbf55989513ddf777ec4ee6c91df)
			if (first) {
				lng = Calloc(nmax, double);
				lng[0] = NAN;
				lng[1] = 0.0;
				for (int i = 2; i < nmax; i++) {
					lng[i] = lng[i - 1] + log((double) (i - 1));
				}
				first = 0;
			}
		}
		if (x >= nmax) {
			return gsl_sf_lngamma(x);
		} else {
			return lng[(int) round(x)];
		}
	}

	assert(0 == 1);
	return NAN;
}

int my_gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
	// copy of gsl_sf_lnfact_e

	result->val = my_gsl_sf_lnfact((int) n);
	result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	return GSL_SUCCESS;
}

int my_gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
	// copy of gsl_sf_lnchoose_e

	if (m > n) {
		assert(0 == 1);
	} else if (m == n || m == 0) {
		result->val = 0.0;
		result->err = 0.0;
	} else {
		gsl_sf_result nf;
		gsl_sf_result mf;
		gsl_sf_result nmmf;
		if (m * 2 > n)
			m = n - m;
		my_gsl_sf_lnfact_e(n, &nf);
		my_gsl_sf_lnfact_e(m, &mf);
		my_gsl_sf_lnfact_e(n - m, &nmmf);
		result->val = nf.val - mf.val - nmmf.val;
		result->err = nf.err + mf.err + nmmf.err;
		result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
	}
	return GMRFLib_SUCCESS;
}

double my_gsl_sf_lnbeta(double a, double b)
{
	double ab_min, ab_max;

	if (b < a) {
		ab_min = b;
		ab_max = a;
	} else {
		ab_min = a;
		ab_max = b;
	}

	ab_min = DMAX(DBL_EPSILON, ab_min);
	double n = ab_max / ab_min;
	const double threshold = 100.0 / DBL_EPSILON;

	if (n > threshold) {
		// log(Beta(a,a/n)) = as n->infinity + symmetry

/*			
 *			asympt(log(Beta(a,a/n)), n,2);
 *                                                         (Psi(a~) + gamma) a~      1
 *                                        ln(n) - ln(a~) - -------------------- + O(----)
 *                                                                  n                 2
 *                                                                                   n
 */

		return (log(n) - log(ab_max) - (gsl_sf_psi(ab_max) + 0.5772156649015328606065120) * ab_max / n);
	} else {
		return (gsl_sf_lnbeta(a, b));
	}
}

double my_betabinomial_helper(int n, double a)
{
	const int roll = 4L;				       /* dont change this */
	double s0 = 0.0;
	div_t d = div(n, roll);
	int m = d.quot * roll;

#pragma GCC ivdep
	for (int i = 0; i < m; i += roll) {
		double aa = i + a;
		s0 += log(aa * (aa + 1.0) * (aa + 2.0) * (aa + 3.0));
	}

	if (d.rem) {
		double aa = m + a;
		switch (d.rem) {
		case 1:
			s0 += log(aa);
			break;
		case 2:
			s0 += log(aa * (aa + 1.0));
			break;
		case 3:
			s0 += log(aa * (aa + 1.0) * (aa + 2.0));
			break;
		}
	}

	return (s0);
}
double my_betabinomial(int y, int n, double a, double b)
{
	double s1 = my_betabinomial_helper(y, a);
	double s2 = my_betabinomial_helper(n - y, b);
	double s3 = my_betabinomial_helper(n, a + b);
	return (s1 + s2 - s3);
}

double my_betabinomial2(int y, int n, double a, double b)
{
	// using Gamma(1+z)=z*Gamma(z), we can get this
	double ladd = 0.0;
	while (a > 1.0) {
		a--;
		ladd += log((y + a) * (a + b) / (n + a + b) / a);
	}
	while (b > 1.0) {
		b--;
		ladd += log((n - y + b) * (a + b) / (n + a + b) / b);
	}

	// here we have 0<a<1, 0<b<1, but NOT a+b<1.
	// this could be helpful creating approximations
	double s1 = my_betabinomial_helper(y, a);
	double s2 = my_betabinomial_helper(n - y, b);
	double s3 = my_betabinomial_helper(n, a + b);
	return (s1 + s2 - s3 + ladd);
}

double  my_lbell(int y) 
{
	// return log(Bell(y) / y!)

#define YMAX 256L
	static double lbell[YMAX+1L] = {
		0.0000000000000000, 
		0.0000000000000000, 
		0.0000000000000000, 
		-0.1823215567939547, 
		-0.4700036292457356, 
		-0.8362480242006187, 
		-1.2660452329683141, 
		-1.7486543686932310, 
		-2.2761518359258899, 
		-2.8425741496394838, 
		-3.4432826434555710, 
		-4.0745649241837487, 
		-4.7333872604910772, 
		-5.4172287809229847, 
		-6.1239644456467079, 
		-6.8517804693538071, 
		-7.5991120266883128, 
		-8.3645964067982472, 
		-9.1470371005827573, 
		-9.9453758204203968, 
		-10.7586703951706593, 
		-11.5860770873270500, 
		-12.4268362863283901, 
		-13.2802608141045102, 
		-14.1457262767463092, 
		-15.0226630367174998, 
		-15.9105494814919997, 
		-16.8089063389180708, 
		-17.7172918449129888, 
		-18.6352976106647183, 
		-19.5625450681058801, 
		-20.4986823966736509, 
		-21.4433818531615188, 
		-22.3963374411598686, 
		-23.3572628681625396, 
		-24.3258897476140916, 
		-25.3019660105298989, 
		-26.2852544972474291, 
		-27.2755317046698309, 
		-28.2725866682795015, 
		-29.2762199614103089, 
		-30.2862427969151113, 
		-31.3024762185581800, 
		-32.3247503712891273, 
		-33.3529038410819112, 
		-34.3867830563058803, 
		-35.4262417436782613, 
		-36.4711404327637112, 
		-37.5213460037659274, 
		-38.5767312740209078, 
		-39.6371746191702101, 
		-40.7025596254814701, 
		-41.7727747702043004, 
		-42.8477131272135523, 
		-43.9272720955071918, 
		-45.0113531484003389, 
		-46.0998616014954621, 
		-47.1927063977178705, 
		-48.2897999078880815, 
		-49.3910577454632715, 
		-50.4963985942212403, 
		-51.6057440477850164, 
		-52.7190184599960006, 
		-53.8361488052414927, 
		-54.9570645479289084, 
		-56.0816975203756130, 
		-57.2099818084524685, 
		-58.3418536443800733, 
		-59.4772513061312935, 
		-60.6161150229431982, 
		-61.7583868864850274, 
		-62.9040107672685167, 
		-64.0529322359226541, 
		-65.2050984889868630, 
		-66.3604582789054405, 
		-67.5189618479325731, 
		-68.6805608656807749, 
		-69.8452083700669135, 
		-71.0128587114300132, 
		-72.1834674996120214, 
		-73.3569915538095501, 
		-74.5333888550189414, 
		-75.7126185009103381, 
		-76.8946406629790999, 
		-78.0794165458335812, 
		-79.2669083484889541, 
		-80.4570792275460036, 
		-81.6498932621423990, 
		-82.8453154205718221, 
		-84.0433115284738221, 
		-85.2438482385037872, 
		-86.4468930013985357, 
		-87.6524140383588275, 
		-88.8603803146754245, 
		-90.0707615145298064, 
		-91.2835280169056489, 
		-92.4986508725507974, 
		-93.7161017819336166, 
		-94.9358530741411499, 
		-96.1578776866695222, 
		-97.3821491460604989, 
		-98.6086415493405752, 
		-99.8373295462219090, 
		-101.0681883220265007, 
		-102.3011935812976958, 
		-103.5363215320655002, 
		-104.7735488707320002, 
		-106.0128527675498020, 
		-107.2542108526614015, 
		-108.4976012026764067, 
		-109.7430023277579068, 
		-110.9903931591966995, 
		-112.2397530374497023, 
		-113.4910617006210032, 
		-114.7442992733660958, 
		-115.9994462562001019, 
		-117.2564835151908937, 
		-118.5153922720217992, 
		-119.7761540944058964, 
		-121.0387508868381019, 
		-122.3031648816692041, 
		-123.5693786304887993, 
		-124.8373749958035006, 
		-126.1071371429987948, 
		-127.3786485325713045, 
		-128.6518929126216904, 
		-129.9268543115969976, 
		-131.2035170312707066, 
		-132.4818656399544068, 
		-133.7618849659276066, 
		-135.0435600910799110, 
		-136.3268763447564993, 
		-137.6118192977978936, 
		-138.8983747567677938, 
		-140.1865287583603958, 
		-141.4762675639821055, 
		-142.7675776544986945, 
		-144.0604457251438930, 
		-145.3548586805812874, 
		-146.6508036301156892, 
		-147.9482678830468956, 
		-149.2472389441618930, 
		-150.5477045093591073, 
		-151.8496524614014049, 
		-153.1530708657917899, 
		-154.4579479667683870, 
		-155.7642721834139934, 
		-157.0720321058767865, 
		-158.3812164916965912, 
		-159.6918142622359937, 
		-161.0038144992093976, 
		-162.3172064413091107, 
		-163.6319794809238033, 
		-164.9481231609470910, 
		-166.2656271716720937, 
		-167.5844813477702928, 
		-168.9046756653517036, 
		-170.2262002391018996, 
		-171.5490453194967984, 
		-172.8732012900889004, 
		-174.1986586648656044, 
		-175.5254080856753092, 
		-176.8534403197207894, 
		-178.1827462571160936, 
		-179.5133169085061979, 
		-180.8451434027466007, 
		-182.1782169846418071, 
		-183.5125290127398046, 
		-184.8480709571821876, 
		-186.1848343976068918, 
		-187.5228110211025125, 
		-188.8619926202127033, 
		-190.2023710909898000, 
		-191.5439384310938067, 
		-192.8866867379387884, 
		-194.2306082068828061, 
		-195.5756951294600015, 
		-196.9219398916563932, 
		-198.2693349722242999, 
		-199.6178729410368078, 
		-200.9675464574808927, 
		-202.3183482688868935, 
		-203.6702712089945919, 
		-205.0233081964537973, 
		-206.3774522333598043, 
		-207.7326964038207961, 
		-209.0890338725587867, 
		-210.4464578835402904, 
		-211.8049617586391093, 
		-213.1645388963270022, 
		-214.5251827703947072, 
		-215.8868869286995107, 
		-217.2496449919414090, 
		-218.6134506524647065, 
		-219.9782976730859048, 
		-221.3441798859468008, 
		-222.7110911913918017, 
		-224.0790255568691123, 
		-225.4479770158550025, 
		-226.8179396668008962, 
		-228.1889076721026015, 
		-229.5608752570902027, 
		-230.9338367090396957, 
		-232.3077863762044899, 
		-233.6827186668668901, 
		-235.0586280484089059, 
		-236.4355090464017906, 
		-237.8133562437143098, 
		-239.1921642796382059, 
		-240.5719278490319937, 
		-241.9526417014812978, 
		-243.3343006404759024, 
		-244.7168995226028017, 
		-246.1004332567555934, 
		-247.4848968033587937, 
		-248.8702851736076127, 
		-250.2565934287219136, 
		-251.6438166792154902, 
		-253.0319500841779075, 
		-254.4209888505714900, 
		-255.8109282325409026, 
		-257.2017635307354908, 
		-258.5934900916452079, 
		-259.9861033069486211, 
		-261.3795986128724849, 
		-262.7739714895641896, 
		-264.1692174604751244, 
		-265.5653320917556925, 
		-266.9623109916611270, 
		-268.3601498099685045, 
		-269.7588442374042756, 
		-271.1583900050816283, 
		-272.5587828839485951, 
		-273.9600186842462790, 
		-275.3620932549753206, 
		-276.7650024833737916, 
		-278.1687422944029890, 
		-279.5733086502423816, 
		-280.9786975497942763, 
		-282.3849050281961013, 
		-283.7919271563417851, 
		-285.1997600404116042, 
		-286.6083998214090798, 
		-288.0178426747074241, 
		-289.4280848096021828, 
		-290.8391224688728016, 
		-292.2509519283505028, 
		-293.6635694964944037, 
		-295.0769715139742857, 
		-296.4911543532604128, 
		-297.9061144182199996, 
		-299.3218481437212972, 
		-300.7383519952427946, 
		-302.1556224684899803, 
		-303.5736560890185274, 
		-304.9924494118620828, 
		-306.4119990211684126};

	if (y >= 0 && y <= YMAX) {
		return (lbell[y]);
	}
	return ((y < 0) ? NAN : -INFINITY);
}
