
/* example-sphere.c
 * 
 * Copyright (C) 2005-2006 Havard Rue
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * The author's contact information:
 *
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */

#if !defined(__FreeBSD__)
#include <math.h>
#include <malloc.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

static const char RCSId[] = "$Id: example-sphere.c,v 1.15 2008/10/26 03:46:18 hrue Exp $";

int main(int argc, char **argv, char **env)
{
	/*
	 * this example demonstrate the use the spherical igmrfs 
	 */

	int seed, i, n, resolution;
	double scale = 1.0, log_scale = 0.0, total, area_ref = 4.0 * M_PI;
	GMRFLib_sphere_tp *sphere = NULL;

	if (argc <= 2) {
		printf("Usage: %s SEED RESOLUTION\n", argv[0]);
	}

	seed = atoi(argv[1]);
	resolution = atoi(argv[2]);

	printf("make spherical igmrf with level=%1d using 'scale'\n", resolution);
	GMRFLib_make_spherical_igmrf(&sphere, resolution, &scale, NULL, NULL);

	/*
	 * now we can use
	 * 
	 * sphere->Qfunc sphere->Qfunc_arg sphere->graph
	 * 
	 * as usual.
	 * 
	 * further, information about each node is available as
	 * 
	 * sphere->area[i] sphere->coord[i].x sphere->coord[i].y sphere->coord[i].z
	 * 
	 * as now demonstrated 
	 */
	n = sphere->graph->n;
	for (i = 0; i < n; i++) {
		printf("%1d: area %f coord (%f,%f,%f)\n", i, sphere->area[i], sphere->coord[i].x, sphere->coord[i].y, sphere->coord[i].z);
	}

	/*
	 * the area should sum up to the surface area for a unit (radius = 1) ball in 3d, ie 4\pi 
	 */

	total = 0;
	for (i = 0; i < n; i++) {
		total += sphere->area[i];
	}
	printf("Relative error in the total area %g\n", (total - area_ref) / area_ref);

	/*
	 * the scale-parameters, constrols the scaling of the Q-matrix, ie, its a precision.
	 * 
	 * note, if &scale = NULL, in GMRFLib_make_spherical_igmrf(), then the scaling is always 1.
	 * 
	 * of'course it works... 
	 */
	scale = 0.01;
	printf("%.3f %.10f\n", scale, sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	scale = 0.1;
	printf("%.3f %.10f\n", scale, sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	scale = 1;
	printf("%.3f %.10f\n", scale, sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	scale = 10;
	printf("%.3f %.10f\n", scale, sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	scale = 100;
	printf("%.3f %.10f\n", scale, sphere->Qfunc(2, 2, sphere->Qfunc_arg));

	/*
	 * this will free the spherical igmrf 
	 */
	GMRFLib_free_spherical_igmrf(sphere);

	printf("make spherical igmrf with level=%1d using 'log_scale'\n", resolution);
	GMRFLib_make_spherical_igmrf(&sphere, resolution, NULL, &log_scale, NULL);

	log_scale = log(0.1);
	printf("%.3f %.10f\n", exp(log_scale), sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	log_scale = 0.0;
	printf("%.3f %.10f\n", exp(log_scale), sphere->Qfunc(2, 2, sphere->Qfunc_arg));
	log_scale = log(10.0);
	printf("%.3f %.10f\n", exp(log_scale), sphere->Qfunc(2, 2, sphere->Qfunc_arg));

	GMRFLib_free_spherical_igmrf(sphere);

	return 0;
}
