
/* moments.c
 * 
 * Copyright (C) 2024-2024 Havard Rue
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


#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"

#define POW2(x_) gsl_pow_2(x_)
#define POW3(x_) gsl_pow_3(x_)
#define POW4(x_) gsl_pow_4(x_)
#define POW5(x_) gsl_pow_5(x_)
#define POW6(x_) gsl_pow_6(x_)
#define POW7(x_) gsl_pow_7(x_)
#define POW8(x_) gsl_pow_8(x_)
#define POW9(x_) gsl_pow_9(x_)
#define POW10(x_) gsl_pow_2(gsl_pow_5(x_))

// compute E(x^i * y^j), either central or non-central

double GMRFLib_noncentral_moment(int i, int j, double mi, double mj, double Sii, double Sij, double Sjj)
{
	double val = NAN;

	if (j > i) {
		return GMRFLib_noncentral_moment(j, i, mj, mi, Sjj, Sij, Sii);
	}

	if (j <= 0) {
		switch (i) {
		case 1:
			val = mi;
			break;

		case 2:
			val = SQR(mi) + Sii;
			break;

		case 3:
			val = POW3(mi) + 3.0 * mi * Sii;
			break;

		case 4:
			val = POW4(mi) + 6.0 * POW2(mi) * Sii + 3.0 * POW2(Sii);
			break;

		case 5:
			val = POW5(mi) + 10.0 * POW3(mi) * Sii + 5.0 * mi * 3.0 * POW2(Sii);
			break;

		case 6:
			val = POW6(mi) + 15 * POW4(mi) * Sii + 15.0 * POW2(mi) * 3.0 * POW2(Sii) + 15.0 * POW3(Sii);
			break;

		case 7:
			val = POW7(mi) + 21.0 * POW5(mi) * Sii + 35.0 * POW3(mi) * 3.0 * POW2(Sii) + 7.0 * mi * 15.0 * POW3(Sii);
			break;

		case 8:
			val =
			    POW8(mi) + 28.0 * POW6(mi) * Sii + 70.0 * POW4(mi) * 3.0 * POW2(Sii) + 28.0 * POW2(mi) * 15.0 * POW3(Sii) +
			    105.0 * POW4(Sii);
			break;

		case 9:
			val =
			    POW9(mi) + 36.0 * POW7(mi) * Sii + 126.0 * POW5(mi) * 3.0 * POW2(Sii) + 84.0 * POW3(mi) * 15.0 * POW3(Sii) +
			    9.0 * mi * 105.0 * POW4(Sii);
			break;

		case 10:
			val =
			    POW10(mi) + 45.0 * POW8(mi) * Sii + 210.0 * POW6(mi) * 3.0 * POW2(Sii) + 210.0 * POW4(mi) * 15.0 * POW3(Sii) +
			    45.0 * POW2(mi) * 105.0 * POW4(Sii) + 945.0 * POW5(Sii);
			break;

		default:
			assert(0 == 1);
		}

		return val;
	}

	switch (i) {
	case 1:
	{
		switch (j) {
		case 1:
			val = mi * mj + Sij;
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 2:
	{
		switch (j) {
		case 1:
			val = POW2(mi) * mj + mj * Sii + 2.0 * mi * Sij;
			break;

		case 2:
			val = POW2(mi) * POW2(mj) + POW2(mj) * Sii + 4.0 * mi * mj * Sij + POW2(mi) * Sjj + Sii * Sjj + 2.0 * POW2(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 3:
	{
		switch (j) {
		case 1:
			val = POW3(mi) * mj + 3.0 * mi * mj * Sii + 3.0 * POW2(mi) * Sij + 3.0 * Sii * Sij;
			break;

		case 2:
			val =
			    POW3(mi) * POW2(mj) + 3.0 * mi * POW2(mj) * Sii + 6.0 * POW2(mi) * mj * Sij + 2.0 * mj * 3.0 * Sii * Sij +
			    POW3(mi) * Sjj + 3.0 * mi * (Sii * Sjj + 2.0 * POW2(Sij));
			break;

		case 3:
			val =
			    POW3(mi) * POW3(mj) + 3.0 * mi * POW3(mj) * Sii + 9.0 * POW2(mi) * POW2(mj) * Sij + 3.0 * POW2(mj) * 3 * Sii * Sij +
			    3 * POW3(mi) * mj * Sjj + 9.0 * mi * mj * (Sii * Sjj + 2.0 * POW2(Sij)) + 3.0 * POW2(mi) * (3.0 * Sij * Sjj) +
			    9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 4:
	{
		switch (j) {
		case 1:
			val = POW4(mi) * mj + 6.0 * POW2(mi) * mj * Sii + mj * 3.0 * POW2(Sii) + 4.0 * POW3(mi) * Sij + 4.0 * mi * 3.0 * Sii * Sij;
			break;

		case 2:
			val =
			    POW4(mi) * POW2(mj) + 6.0 * POW2(mi) * POW2(mj) * Sii + POW2(mj) * 3.0 * POW2(Sii) + 8.0 * POW3(mi) * mj * Sij +
			    8.0 * mi * mj * 3.0 * Sii * Sij + POW4(mi) * Sjj + 6.0 * POW2(mi) * (Sii * Sjj + 2.0 * POW2(Sij)) +
			    3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij);
			break;

		case 3:
			val =
			    POW4(mi) * POW3(mj) + 6.0 * POW2(mi) * POW3(mj) * Sii + POW3(mj) * 3.0 * POW2(Sii) + 12.0 * POW3(mi) * POW2(mj) * Sij +
			    12.0 * mi * POW2(mj) * 3 * Sii * Sij + 3.0 * POW4(mi) * mj * Sjj + 18.0 * POW2(mi) * mj * (Sii * Sjj +
														       2.0 * POW2(Sij)) +
			    3 * mj * (3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij)) + 4.0 * POW3(mi) * 3.0 * Sij * Sjj +
			    4.0 * mi * (9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij));
			break;

		case 4:
			val =
			    POW4(mi) * POW4(mj) + 6.0 * POW2(mi) * POW4(mj) * Sii + POW4(mj) * 3.0 * POW2(Sii) + 16.0 * POW3(mi) * POW3(mj) * Sij +
			    16 * mi * POW3(mj) * 3.0 * Sii * Sij + 6.0 * POW4(mi) * POW2(mj) * Sjj + 36.0 * POW2(mi) * POW2(mj) * (Sii * Sjj +
																   2.0 *
																   POW2(Sij)) +
			    6.0 * POW2(mj) * (3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij)) + 16.0 * POW3(mi) * mj * 3.0 * Sij * Sjj +
			    16.0 * mi * mj * (9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij)) + POW4(mi) * 3.0 * POW2(Sjj) +
			    6.0 * POW2(mi) * (3.0 * Sii * POW2(Sjj) + 12.0 * POW2(Sij) * Sjj) + 9.0 * POW2(Sii) * POW2(Sjj) +
			    72.0 * Sii * POW2(Sij) * Sjj + 24.0 * POW4(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 5:
	{
		switch (j) {
		case 1:
			val =
			    POW5(mi) * mj + 10.0 * POW3(mi) * mj * Sii + 5.0 * mi * mj * 3.0 * POW2(Sii) + 5 * POW4(mi) * Sij +
			    10 * POW2(mi) * 3.0 * Sii * Sij + 15.0 * POW2(Sii) * Sij;
			break;

		case 2:
			val =
			    POW5(mi) * POW2(mj) + 10.0 * POW3(mi) * POW2(mj) * Sii + 5.0 * mi * POW2(mj) * 3.0 * POW2(Sii) +
			    10.0 * POW4(mi) * mj * Sij + 20.0 * POW2(mi) * mj * 3.0 * Sii * Sij + 2.0 * mj * 15.0 * POW2(Sii) * Sij +
			    POW5(mi) * Sjj + 10.0 * POW3(mi) * (Sii * Sjj + 2.0 * POW2(Sij)) + 5 * mi * (3.0 * POW2(Sii) * Sjj +
													 12.0 * Sii * POW2(Sij));
			break;

		case 3:
			val =
			    POW5(mi) * POW3(mj) + 10.0 * POW3(mi) * POW3(mj) * Sii + 5 * mi * POW3(mj) * 3.0 * POW2(Sii) +
			    15 * POW4(mi) * POW2(mj) * Sij + 30 * POW2(mi) * POW2(mj) * (3.0 * Sii * Sij) + 3 * POW2(mj) * 15.0 * POW2(Sii) * Sij +
			    3 * POW5(mi) * mj * Sjj + 30 * POW3(mi) * mj * (Sii * Sjj + 2.0 * POW2(Sij)) + 15 * mi * mj * (3.0 * POW2(Sii) * Sjj +
															   12.0 * Sii * POW2(Sij)) +
			    5 * POW4(mi) * 3.0 * Sij * Sjj + 10.0 * POW2(mi) * (9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij)) +
			    45.0 * POW2(Sii) * Sij * Sjj + 60.0 * Sii * POW3(Sij);
			break;

		case 4:
			val =
			    POW5(mi) * POW4(mj) + 10.0 * POW3(mi) * POW4(mj) * Sii + 5.0 * mi * POW4(mj) * 3.0 * POW2(Sii) +
			    20.0 * POW4(mi) * POW3(mj) * Sij + 40.0 * POW2(mi) * POW3(mj) * 3.0 * Sii * Sij +
			    4.0 * POW3(mj) * 15.0 * POW2(Sii) * Sij + 6.0 * POW5(mi) * POW2(mj) * Sjj + 60.0 * POW3(mi) * POW2(mj) * (Sii * Sjj +
																      2.0 *
																      POW2(Sij)) +
			    30.0 * mi * POW2(mj) * (3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij)) + 20.0 * POW4(mi) * mj * 3.0 * Sij * Sjj +
			    40 * POW2(mi) * mj * (9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij)) + 4 * mj * (45.0 * POW2(Sii) * Sij * Sjj +
												       60.0 * Sii * POW3(Sij)) +
			    POW5(mi) * 3.0 * POW2(Sjj) + 10.0 * POW3(mi) * (3.0 * Sii * POW2(Sjj) + 12.0 * POW2(Sij) * Sjj) +
			    5.0 * mi * (9.0 * POW2(Sii) * POW2(Sjj) + 72.0 * Sii * POW2(Sij) * Sjj + 24 * POW4(Sij));
			break;

		case 5:
			val =
			    POW5(mi) * POW5(mj) + 10.0 * POW3(mi) * POW5(mj) * Sii + 5 * mi * POW5(mj) * 3 * POW2(Sii) +
			    25 * POW4(mi) * POW4(mj) * Sij + 50.0 * POW2(mi) * POW4(mj) * 3.0 * Sii * Sij +
			    5.0 * POW4(mj) * 15.0 * POW2(Sii) * Sij + 10.0 * POW5(mi) * POW3(mj) * Sjj + 100.0 * POW3(mi) * POW3(mj) * (Sii * Sjj +
																	2.0 *
																	POW2(Sij))
			    + 50 * mi * POW3(mj) * (3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij)) + 50.0 * POW4(mi) * POW2(mj) * 3.0 * Sij * Sjj +
			    100.0 * POW2(mi) * POW2(mj) * (9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij)) +
			    10 * POW2(mj) * (45.0 * POW2(Sii) * Sij * Sjj + 60.0 * Sii * POW3(Sij)) + 5.0 * POW5(mi) * mj * 3 * POW2(Sjj) +
			    50.0 * POW3(mi) * mj * (3.0 * Sii * POW2(Sjj) + 12.0 * POW2(Sij) * Sjj) +
			    25.0 * mi * mj * (9.0 * POW2(Sii) * POW2(Sjj) + 72.0 * Sii * POW2(Sij) * Sjj + 24.0 * POW4(Sij)) +
			    5 * POW4(mi) * 15.0 * Sij * POW2(Sjj) + 10.0 * POW2(mi) * (45.0 * Sii * Sij * POW2(Sjj) + 60.0 * POW3(Sij) * Sjj) +
			    225.0 * POW2(Sii) * Sij * POW2(Sjj) + 600 * Sii * POW3(Sij) * Sjj + 120.0 * POW5(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	default:
		assert(0 == 1);
		break;
	}

	return val;
}

double GMRFLib_central_moment(int i, int j, double Sii, double Sij, double Sjj)
{
	double val = NAN;

	if (j > i) {
		return GMRFLib_central_moment(j, i, Sjj, Sij, Sii);
	}

	if (j <= 0) {
		switch (i) {
		case 1:
			val = 0.0;
			break;

		case 2:
			val = Sii;
			break;

		case 3:
			val = 0.0;
			break;

		case 4:
			val = 3.0 * POW2(Sii);
			break;

		case 5:
			val = 0.0;
			break;

		case 6:
			val = 15.0 * POW3(Sii);
			break;

		case 7:
			val = 0.0;
			break;

		case 8:
			val = 105.0 * POW4(Sii);
			break;

		case 9:
			val = 0.0;
			break;

		case 10:
			val = 945.0 * POW5(Sii);
			break;

		default:
			assert(0 == 1);
		}

		return val;
	}

	switch (i) {
	case 1:
	{
		switch (j) {
		case 1:
			val = Sij;
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 2:
	{
		switch (j) {
		case 1:
			val = 0.0;
			break;

		case 2:
			val = Sii * Sjj + 2.0 * POW2(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 3:
	{
		switch (j) {
		case 1:
			val = 3.0 * Sii * Sij;
			break;

		case 2:
			val = 0.0;
			break;

		case 3:
			val = 9.0 * Sii * Sij * Sjj + 6.0 * POW3(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 4:
	{
		switch (j) {
		case 1:
			val = 0.0;
			break;

		case 2:
			val = 3.0 * POW2(Sii) * Sjj + 12.0 * Sii * POW2(Sij);
			break;

		case 3:
			val = 0.0;
			break;

		case 4:
			val = 9.0 * POW2(Sii) * POW2(Sjj) + 72.0 * Sii * POW2(Sij) * Sjj + 24.0 * POW4(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	case 5:
	{
		switch (j) {
		case 1:
			val = 15.0 * POW2(Sii) * Sij;
			break;

		case 2:
			val = 0.0;
			break;

		case 3:
			val = 45.0 * POW2(Sii) * Sij * Sjj + 60.0 * Sii * POW3(Sij);
			break;

		case 4:
			val = 0.0;
			break;

		case 5:
			val = 225.0 * POW2(Sii) * Sij * POW2(Sjj) + 600 * Sii * POW3(Sij) * Sjj + 120.0 * POW5(Sij);
			break;

		default:
			assert(0 == 1);
			break;
		}
		return val;
	}
		break;

	default:
		assert(0 == 1);
		break;
	}

	return val;
}
