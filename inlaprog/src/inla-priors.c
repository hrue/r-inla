
/* inla-priors.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
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


int inla_read_prior(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "PRIOR", "PARAMETERS", "FROM.THETA", "TO.THETA", "HYPERID", default_prior, args);
}

int inla_read_prior_mix(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "MIX.PRIOR", "MIX.PARAMETERS", "MIX.FROM.THETA", "MIX.TO.THETA",
				       "MIX.HYPERID", default_prior, args);
}

int inla_read_prior_link(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "LINK.PRIOR", "LINK.PARAMETERS", "LINK.FROM.THETA", "LINK.TO.THETA",
				       "LINK.HYPERID", default_prior, args);
}

int inla_read_prior_link0(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "LINK.PRIOR0", "LINK.PARAMETERS0", "LINK.FROM.THETA0", "LINK.TO.THETA0",
				       "LINK.HYPERID0", default_prior, args);
}

int inla_read_prior_link1(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "LINK.PRIOR1", "LINK.PARAMETERS1", "LINK.FROM.THETA1", "LINK.TO.THETA1",
				       "LINK.HYPERID1", default_prior, args);
}

int inla_read_prior_link2(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "LINK.PRIOR2", "LINK.PARAMETERS2", "LINK.FROM.THETA2", "LINK.TO.THETA2",
				       "LINK.HYPERID2", default_prior, args);
}

int inla_read_prior_group(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR", "GROUP.PARAMETERS", "GROUP.FROM.THETA", "GROUP.TO.THETA",
				       "GROUP.HYPERID", default_prior, args);
}

int inla_read_prior_group0(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR0", "GROUP.PARAMETERS0", "GROUP.FROM.THETA0",
				       "GROUP.TO.THETA0", "GROUP.HYPERID0", default_prior, args);
}

int inla_read_prior_group1(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR1", "GROUP.PARAMETERS1", "GROUP.FROM.THETA1",
				       "GROUP.TO.THETA1", "GROUP.HYPERID1", default_prior, args);
}

int inla_read_prior_group2(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR2", "GROUP.PARAMETERS2", "GROUP.FROM.THETA2",
				       "GROUP.TO.THETA2", "GROUP.HYPERID2", default_prior, args);
}

int inla_read_prior_group3(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR3", "GROUP.PARAMETERS3", "GROUP.FROM.THETA3",
				       "GROUP.TO.THETA3", "GROUP.HYPERID3", default_prior, args);
}

int inla_read_prior_group4(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR4", "GROUP.PARAMETERS4", "GROUP.FROM.THETA4",
				       "GROUP.TO.THETA4", "GROUP.HYPERID4", default_prior, args);
}

int inla_read_prior_group5(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR5", "GROUP.PARAMETERS5", "GROUP.FROM.THETA5",
				       "GROUP.TO.THETA5", "GROUP.HYPERID5", default_prior, args);
}

int inla_read_prior_group6(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR6", "GROUP.PARAMETERS6", "GROUP.FROM.THETA6",
				       "GROUP.TO.THETA6", "GROUP.HYPERID6", default_prior, args);
}

int inla_read_prior_group7(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR7", "GROUP.PARAMETERS7", "GROUP.FROM.THETA7",
				       "GROUP.TO.THETA7", "GROUP.HYPERID7", default_prior, args);
}

int inla_read_prior_group8(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR8", "GROUP.PARAMETERS8", "GROUP.FROM.THETA8",
				       "GROUP.TO.THETA8", "GROUP.HYPERID8", default_prior, args);
}

int inla_read_prior_group9(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR9", "GROUP.PARAMETERS9", "GROUP.FROM.THETA9",
				       "GROUP.TO.THETA9", "GROUP.HYPERID9", default_prior, args);
}

int inla_read_prior_group10(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_prior_generic(mb, ini, sec, prior, "GROUP.PRIOR10", "GROUP.PARAMETERS10", "GROUP.FROM.THETA10",
				       "GROUP.TO.THETA10", "GROUP.HYPERID10", default_prior, args);
}

int inla_read_prior0(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 0, args);
}

int inla_read_prior1(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 1, args);
}

int inla_read_prior2(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 2, args);
}

int inla_read_prior3(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 3, args);
}

int inla_read_prior4(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 4, args);
}

int inla_read_prior5(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 5, args);
}

int inla_read_prior6(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 6, args);
}

int inla_read_prior7(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 7, args);
}

int inla_read_prior8(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 8, args);
}

int inla_read_prior9(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 9, args);
}

int inla_read_prior10(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, void *args)
{
	return inla_read_priorN(mb, ini, sec, prior, default_prior, 10, args);
}

int inla_read_priorN(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *default_prior, int N, void *args)
{
	int val;
	char *a = NULL, *b = NULL, *c = NULL, *d = NULL, *e = NULL;

	GMRFLib_sprintf(&a, "PRIOR%1d", N);
	GMRFLib_sprintf(&b, "PARAMETERS%1d", N);
	GMRFLib_sprintf(&c, "FROM.THETA%1d", N);
	GMRFLib_sprintf(&d, "TO.THETA%1d", N);
	GMRFLib_sprintf(&e, "HYPERID%1d", N);
	val =
	    inla_read_prior_generic(mb, ini, sec, prior, (const char *) a, (const char *) b, (const char *) c, (const char *) d,
				    (const char *) e, default_prior, args);
	Free(a);
	Free(b);
	Free(c);
	Free(d);
	Free(e);

	return val;
}

int inla_read_prior_generic(inla_tp *mb, dictionary *ini, int sec, Prior_tp *prior, const char *prior_tag,
			    const char *param_tag, const char *from_theta, const char *to_theta, const char *hyperid, const char *default_prior,
			    void *UNUSED(args))
{
	char *secname = NULL, *param = NULL;
	secname = Strdup(iniparser_getsecname(ini, sec));
	prior->name = Strdup(iniparser_getstring(ini, inla_string_join(secname, prior_tag), Strdup(default_prior)));

	if (!prior->name) {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, prior_tag, NULL);
	}

	prior->priorfunc = NULL;
	prior->expression = NULL;

	if (mb->verbose) {
		/*
		 * remove trailing -[a-zA-Z]*$ 
		 */
		char *p, *new_name;
		new_name = Strdup(prior->name);
		p = GMRFLib_rindex((const char *) new_name, '-');
		if (p) {
			*p = '\0';
		}
		p = GMRFLib_rindex((const char *) new_name, ':');
		if (p) {
			*p = '\0';
		}
		printf("\t\t%s->name=[%s]\n", prior_tag, new_name);
		Free(new_name);
	}

	prior->hyperid = inla_create_hyperid(iniparser_getint(ini, inla_string_join(secname, hyperid), 0), secname);
	prior->from_theta = Strdup(iniparser_getstring(ini, inla_string_join(secname, from_theta), NULL));
	prior->to_theta = Strdup(iniparser_getstring(ini, inla_string_join(secname, to_theta), NULL));
	if (mb->verbose) {
		printf("\t\thyperid=[%s]\n", prior->hyperid);
		printf("\t\t%s->from_theta=[%s]\n", prior_tag, prior->from_theta);
		printf("\t\t%s->to_theta = [%s]\n", prior_tag, prior->to_theta);
	}
	// this is special
	if (!strcasecmp(prior->name, "WISHARTKD")) {
		prior->name = Strdup(default_prior);
	}

	param = Strdup(iniparser_getstring(ini, inla_string_join(secname, param_tag), NULL));
	if (!strcasecmp(prior->name, "LAPLACE")) {
		prior->id = P_LAPLACE;
		prior->priorfunc = priorfunc_laplace;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0;
			prior->parameters[1] = sqrt(2.0 * DEFAULT_NORMAL_PRIOR_PRECISION);
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAMMA")) {
		prior->id = P_GAMMA;
		prior->priorfunc = priorfunc_gamma;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = DEFAULT_GAMMA_PRIOR_A;
			prior->parameters[1] = DEFAULT_GAMMA_PRIOR_B;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "LOGGAMMA")) {
		prior->id = P_LOGGAMMA;
		prior->priorfunc = priorfunc_loggamma;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = DEFAULT_GAMMA_PRIOR_A;
			prior->parameters[1] = DEFAULT_GAMMA_PRIOR_B;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "LOGGAMMA-ALPHA")) {
		prior->id = P_LOGGAMMA;
		prior->priorfunc = priorfunc_loggamma;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 25.0;
			prior->parameters[1] = 25.0;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN") || !strcasecmp(prior->name, "NORMAL")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;	       /* mean */
			prior->parameters[1] = DEFAULT_NORMAL_PRIOR_PRECISION;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN-1") || !strcasecmp(prior->name, "NORMAL-1")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 1.0;	       /* mean */
			prior->parameters[1] = 1.0;	       /* precision */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN-a") || !strcasecmp(prior->name, "NORMAL-a")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;	       /* mean */
			prior->parameters[1] = 6.25;	       /* precision */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN-std") || !strcasecmp(prior->name, "NORMAL-std")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;	       /* mean */
			prior->parameters[1] = 1.0;	       /* precision */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN-rho") || !strcasecmp(prior->name, "NORMAL-rho")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;	       /* mean */
			prior->parameters[1] = 0.2;	       /* precision */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "GAUSSIAN-group") || !strcasecmp(prior->name, "NORMAL-group")) {
		prior->id = P_GAUSSIAN;
		prior->priorfunc = priorfunc_gaussian;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;	       /* mean */
			prior->parameters[1] = 0.2;	       /* precision */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g, %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "MVNORM") || !strcasecmp(prior->name, "MVGAUSSIAN")) {
		int nparam = 0, i, dim;
		double *tmp;

		prior->id = P_MVNORM;
		prior->priorfunc = priorfunc_mvnorm;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
		/*
		 * add the dimension of the mvgaussian as the first argument of the parameter 
		 */
		dim = -1;
		for (i = 0;; i++) {			       /* yes, an infinite loop */
			if (nparam == i + ISQR(i)) {
				dim = i;
				break;
			}
			if (nparam < i + ISQR(i)) {
				inla_error_general("nparam does not match with any dimension of the mvnorm");
				exit(EXIT_FAILURE);
			}
		}
		tmp = Calloc(nparam + 1, double);
		tmp[0] = dim;
		Memcpy(&(tmp[1]), prior->parameters, nparam * sizeof(double));
		Free(prior->parameters);
		prior->parameters = tmp;
	} else if (!strcasecmp(prior->name, "PCSPDEGA")) {
		int nparam, i;

		prior->id = P_PC_SPDE_GA;
		prior->priorfunc = priorfunc_pc_spde_ga;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 4);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCMATERN")) {
		int nparam, i;

		prior->id = P_PC_MATERN;
		prior->priorfunc = priorfunc_pc_matern;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 3);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCRANGE")) {
		int nparam, i;

		prior->id = P_PC_RANGE;
		prior->priorfunc = priorfunc_pc_range;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 2);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCGEVTAIL")) {
		int nparam, i;

		prior->id = P_PC_GEVTAIL;
		prior->priorfunc = priorfunc_pc_gevtail;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 3);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCGAMMA")) {
		int nparam, i;

		prior->id = P_PC_GAMMA;
		prior->priorfunc = priorfunc_pc_gamma;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 1);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCALPHAW")) {
		int nparam, i;

		prior->id = P_PC_ALPHAW;
		prior->priorfunc = priorfunc_pc_alphaw;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 1);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCMGAMMA")) {
		int nparam, i;

		prior->id = P_PC_MGAMMA;
		prior->priorfunc = priorfunc_pc_mgamma;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 1);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "PCGAMMACOUNT")) {
		int nparam, i;

		prior->id = P_PC_GAMMACOUNT;
		prior->priorfunc = priorfunc_pc_gammacount;
		inla_sread_doubles_q(&(prior->parameters), &nparam, param);
		assert(nparam == 1);
		if (mb->verbose) {
			for (i = 0; i < nparam; i++) {
				printf("\t\t%s->%s[%1d]=[%g]\n", prior_tag, param_tag, i, prior->parameters[i]);
			}
		}
	} else if (!strcasecmp(prior->name, "MINUSLOGSQRTRUNCNORMAL") || !strcasecmp(prior->name, "MINUSLOGSQRTRUNCGAUSSIAN") ||
		   // easier names...
		   !strcasecmp(prior->name, "LOGTNORMAL") || !strcasecmp(prior->name, "LOGTGAUSSIAN")) {
		prior->id = P_MINUSLOGSQRTRUNCGAUSSIAN;
		prior->priorfunc = priorfunc_minuslogsqrtruncnormal;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;
			prior->parameters[1] = DEFAULT_NORMAL_PRIOR_PRECISION;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "FLAT") || !strcasecmp(prior->name, "UNIFORM")) {
		/*
		 * do not care about the range for the FLAT/UNI prior, as the parameters are already transformed to R. 
		 */
		prior->id = P_FLAT;
		prior->priorfunc = priorfunc_flat;
		prior->parameters = NULL;
	} else if (!strcasecmp(prior->name, "INVALID")) {
		prior->id = P_INVALID;
		prior->priorfunc = priorfunc_invalid;
		prior->parameters = NULL;
	} else if (!strcasecmp(prior->name, "WISHART1D") ||
		   !strcasecmp(prior->name, "WISHART2D") || !strcasecmp(prior->name, "WISHART3D") || !strcasecmp(prior->name, "WISHART4D")
		   || !strcasecmp(prior->name, "WISHART5D")) {

		if (!strcasecmp(prior->name, "WISHART1D")) {
			prior->id = P_WISHART1D;
			prior->priorfunc = priorfunc_wishart1d;
		} else if (!strcasecmp(prior->name, "WISHART2D")) {
			prior->id = P_WISHART2D;
			prior->priorfunc = priorfunc_wishart2d;
		} else if (!strcasecmp(prior->name, "WISHART3D")) {
			prior->id = P_WISHART3D;
			prior->priorfunc = priorfunc_wishart3d;
		} else if (!strcasecmp(prior->name, "WISHART4D")) {
			prior->id = P_WISHART4D;
			prior->priorfunc = priorfunc_wishart4d;
		} else if (!strcasecmp(prior->name, "WISHART5D")) {
			prior->id = P_WISHART5D;
			prior->priorfunc = priorfunc_wishart5d;
		} else {
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}

		double *xx = NULL;
		int nxx;
		int idim = (!strcasecmp(prior->name, "WISHART1D") ? 1 :
			    (!strcasecmp(prior->name, "WISHART2D") ? 2 :
			     (!strcasecmp(prior->name, "WISHART3D") ? 3 :
			      (!strcasecmp(prior->name, "WISHART4D") ? 4 : (!strcasecmp(prior->name, "WISHART5D") ? 5 : -1)))));
		assert(idim > 0);

		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx == inla_iid_wishart_nparam(idim) + 1);	/* this must be TRUE */

		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK2D")) {
		prior->id = P_WISHARTK_2D;
		prior->priorfunc = priorfunc_wishartk_2d;

		double *xx = NULL;
		int nxx;
		int idim = 2;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK3D")) {
		prior->id = P_WISHARTK_3D;
		prior->priorfunc = priorfunc_wishartk_3d;

		double *xx = NULL;
		int nxx;
		int idim = 3;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK4D")) {
		prior->id = P_WISHARTK_4D;
		prior->priorfunc = priorfunc_wishartk_4d;

		double *xx = NULL;
		int nxx;
		int idim = 4;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK5D")) {
		prior->id = P_WISHARTK_5D;
		prior->priorfunc = priorfunc_wishartk_5d;

		double *xx = NULL;
		int nxx;
		int idim = 5;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK6D")) {
		prior->id = P_WISHARTK_6D;
		prior->priorfunc = priorfunc_wishartk_6d;

		double *xx = NULL;
		int nxx;
		int idim = 6;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK7D")) {
		prior->id = P_WISHARTK_7D;
		prior->priorfunc = priorfunc_wishartk_7d;

		double *xx = NULL;
		int nxx;
		int idim = 7;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK8D")) {
		prior->id = P_WISHARTK_8D;
		prior->priorfunc = priorfunc_wishartk_8d;

		double *xx = NULL;
		int nxx;
		int idim = 8;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK9D")) {
		prior->id = P_WISHARTK_9D;
		prior->priorfunc = priorfunc_wishartk_9d;

		double *xx = NULL;
		int nxx;
		int idim = 9;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK10D")) {
		prior->id = P_WISHARTK_10D;
		prior->priorfunc = priorfunc_wishartk_10d;

		double *xx = NULL;
		int nxx;
		int idim = 10;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK11D")) {
		prior->id = P_WISHARTK_11D;
		prior->priorfunc = priorfunc_wishartk_11d;

		double *xx = NULL;
		int nxx;
		int idim = 11;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK12D")) {
		prior->id = P_WISHARTK_12D;
		prior->priorfunc = priorfunc_wishartk_12d;

		double *xx = NULL;
		int nxx;
		int idim = 12;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK13D")) {
		prior->id = P_WISHARTK_13D;
		prior->priorfunc = priorfunc_wishartk_13d;

		double *xx = NULL;
		int nxx;
		int idim = 13;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK14D")) {
		prior->id = P_WISHARTK_14D;
		prior->priorfunc = priorfunc_wishartk_14d;

		double *xx = NULL;
		int nxx;
		int idim = 14;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK15D")) {
		prior->id = P_WISHARTK_15D;
		prior->priorfunc = priorfunc_wishartk_15d;

		double *xx = NULL;
		int nxx;
		int idim = 15;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK16D")) {
		prior->id = P_WISHARTK_16D;
		prior->priorfunc = priorfunc_wishartk_16d;

		double *xx = NULL;
		int nxx;
		int idim = 16;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK17D")) {
		prior->id = P_WISHARTK_17D;
		prior->priorfunc = priorfunc_wishartk_17d;

		double *xx = NULL;
		int nxx;
		int idim = 17;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK18D")) {
		prior->id = P_WISHARTK_18D;
		prior->priorfunc = priorfunc_wishartk_18d;

		double *xx = NULL;
		int nxx;
		int idim = 18;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK19D")) {
		prior->id = P_WISHARTK_19D;
		prior->priorfunc = priorfunc_wishartk_19d;

		double *xx = NULL;
		int nxx;
		int idim = 19;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "WISHARTK20D")) {
		prior->id = P_WISHARTK_20D;
		prior->priorfunc = priorfunc_wishartk_20d;

		double *xx = NULL;
		int nxx;
		int idim = 20;
		inla_sread_doubles_q(&xx, &nxx, param);
		assert(xx);
		prior->parameters = xx;
		assert(nxx >= INLA_WISHARTK_NPARAM(idim));
		nxx = INLA_WISHARTK_NPARAM(idim);
		for (int j = 1; j < nxx; j++) {
			if (INLA_IS_SPECIAL(xx[j])) {
				if (j < idim + 1) {
					xx[j] = 1.0;
				} else {
					xx[j] = 0.0;
				}
			}
		}
		if (mb->verbose) {
			int ii;
			for (ii = 0; ii < nxx; ii++) {
				printf("\t\t%s->%s prior_parameter[%1d] = %g\n", prior_tag, param_tag, ii, prior->parameters[ii]);
			}
		}
	} else if (!strcasecmp(prior->name, "LOGFLAT")) {
		prior->id = P_LOGFLAT;
		prior->priorfunc = priorfunc_logflat;
		prior->parameters = NULL;
		if (mb->verbose) {
			printf("\t\t%s->%s=[]\n", prior_tag, param_tag);
		}
	} else if (!strcasecmp(prior->name, "LOGIFLAT")) {
		prior->id = P_LOGIFLAT;
		prior->priorfunc = priorfunc_logiflat;
		prior->parameters = NULL;
		if (mb->verbose) {
			printf("\t\t%s->%s=[]\n", prior_tag, param_tag);
		}
	} else if (!strcasecmp(prior->name, "NONE")) {
		prior->id = P_NONE;
		prior->priorfunc = NULL;
		prior->parameters = NULL;
		if (mb->verbose) {
			printf("\t\t%s->%s=[]\n", prior_tag, param_tag);
		}
	} else if (!strcasecmp(prior->name, "PCFGNH")) {
		prior->id = P_PC_FGN_H;
		prior->priorfunc = priorfunc_fgn_priorH;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.9;
			prior->parameters[1] = 0.1;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "BETACORRELATION")) {
		prior->id = P_BETACORRELATION;
		prior->priorfunc = priorfunc_betacorrelation;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 5.0;
			prior->parameters[1] = 5.0;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "LOGITBETA")) {
		prior->id = P_LOGITBETA;
		prior->priorfunc = priorfunc_betacorrelation;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 1.0;
			prior->parameters[1] = 1.0;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCPREC")) {
		prior->id = P_PC_PREC;
		prior->priorfunc = priorfunc_pc_prec;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.1;	       /* u */
			prior->parameters[1] = 0.001;	       /* alpha */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCDOF")) {
		prior->id = P_PC_DOF;
		prior->priorfunc = priorfunc_pc_dof;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double); /* Prob(dof < u) = alpha */
			prior->parameters[0] = 10.0;	       /* u */
			prior->parameters[1] = 0.5;	       /* alpha */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCSN")) {
		prior->id = P_PC_SN;
		prior->priorfunc = priorfunc_pc_sn;
		if (param && inla_is_NAs(1, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(1, double);
			if (inla_sread_doubles(prior->parameters, 1, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(1, double);
			prior->parameters[0] = 10.0;	       /* lambda */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g]\n", prior_tag, param_tag, prior->parameters[0]);
		}
	} else if (!strcasecmp(prior->name, "LINKSNINTERCEPT")) {
		prior->id = P_SN_INTERCEPT;
		prior->priorfunc = priorfunc_linksnintercept;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.0;
			prior->parameters[1] = 0.0;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "LOGITBETA")) {
		prior->id = P_LOGITBETA;
		prior->priorfunc = priorfunc_logitbeta;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 5.0;
			prior->parameters[1] = 5.0;
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCCOR0")) {
		prior->id = P_PC_COR0;
		prior->priorfunc = priorfunc_pc_cor0;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.5;	       /* u */
			prior->parameters[1] = 0.5;	       /* alpha */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCCOR1")) {
		prior->id = P_PC_COR1;
		prior->priorfunc = priorfunc_pc_cor1;
		if (param && inla_is_NAs(2, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			if (inla_sread_doubles(prior->parameters, 2, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 0.5;	       /* u */
			prior->parameters[1] = 0.5;	       /* alpha */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g %g]\n", prior_tag, param_tag, prior->parameters[0], prior->parameters[1]);
		}
	} else if (!strcasecmp(prior->name, "PCAR")) {
		prior->id = P_PC_AR;
		prior->priorfunc = priorfunc_pc_ar;
		if (param && inla_is_NAs(1, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(2, double);
			double tmp = 0.0;
			if (inla_sread_doubles(&tmp, 1, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
			prior->parameters[0] = tmp;	       /* lambda */
			prior->parameters[1] = -1;	       /* ORDER: to be decided */
		} else {
			prior->parameters = Calloc(2, double);
			prior->parameters[0] = 1.0;	       /* lambda */
			prior->parameters[1] = -1;	       /* ORDER: to be decided */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g]\n", prior_tag, param_tag, prior->parameters[0]);
		}
	} else if (!strcasecmp(prior->name, "DIRICHLET")) {
		prior->id = P_DIRICHLET;
		prior->priorfunc = priorfunc_dirichlet;
		if (param && inla_is_NAs(1, param) != GMRFLib_SUCCESS) {
			prior->parameters = Calloc(3, double); /* yes, 3 */
			if (inla_sread_doubles(prior->parameters, 1, param) == INLA_FAIL) {
				inla_error_field_is_void(__GMRFLib_FuncName, secname, param_tag, param);
			}
			prior->parameters[1] = NAN;	       /* number of classes, added later */
			prior->parameters[2] = NAN;	       /* cdf, added later */
		} else {
			prior->parameters = Calloc(3, double);
			prior->parameters[0] = 0.5;	       /* alpha */
			prior->parameters[1] = NAN;	       /* number of classes, added later */
			prior->parameters[2] = NAN;	       /* number of classes, added later */
		}
		if (mb->verbose) {
			printf("\t\t%s->%s=[%g]\n", prior_tag, param_tag, prior->parameters[0]);
		}
	} else if (!strcasecmp(prior->name, "REFAR")) {
		prior->id = P_REF_AR;
		prior->priorfunc = priorfunc_ref_ar;
		prior->parameters = Calloc(1, double);
		prior->parameters[0] = -1;		       /* ORDER: to be decided */
		if (mb->verbose) {
			printf("\t\t%s->%s=[NULL]\n", prior_tag, param_tag);
		}
	} else if (!strncasecmp(prior->name, "EXPRESSION:", strlen("EXPRESSION:"))) {
		prior->id = P_EXPRESSION;
		prior->expression = Strdup(prior->name);
		prior->name[strlen("EXPRESSION")] = '\0';
		prior->parameters = NULL;

		if (mb->verbose) {
			printf("\t\t%s->%s=[%s]\n", prior_tag, prior->name, prior->expression);
		}
	} else if (!strncasecmp(prior->name, "TABLE:", strlen("TABLE:"))) {
		prior->id = P_TABLE;
		prior->expression = Strdup(prior->name);       /* yes, use the same storage */
		prior->name[strlen("TABLE")] = '\0';
		prior->parameters = NULL;

		if (mb->verbose) {
			printf("\t\t%s->%s=[%s]\n", prior_tag, prior->name, prior->expression);
		}
	} else if (!strncasecmp(prior->name, "RPRIOR:", strlen("RPRIOR:"))) {
		prior->id = P_RPRIOR;
		char *tag[3] = { NULL, NULL, NULL };
		inla_sread_str_str(tag, 3, prior->name);       // rprior:RPRIOR_FUNCTION:FILENAME
		prior->rprior = Strdup(tag[1]);
		GMRFLib_sprintf(&(prior->name), "%s:%s", tag[0], tag[1]);
		GMRFLib_sprintf(&(prior->expression), "%s:%s", tag[0], tag[1]);
		prior->parameters = NULL;

#pragma omp critical (Name_efff84617ceb85320978c2deffcb5b8433bc888d)
		{
			inla_R_init_();
			inla_R_load(tag[2]);
		}

		if (mb->verbose) {
			printf("\t\t%s->%s=[%s]\n", prior_tag, prior->name, prior->rprior);
		}
	} else if (!strcasecmp(prior->name, "JEFFREYSTDF")) {
		prior->id = P_JEFFREYS_T_DF;
		prior->priorfunc = priorfunc_jeffreys_df_student_t;
		prior->parameters = NULL;
		if (mb->verbose) {
			printf("\t\t%s->%s=[NULL]\n", prior_tag, param_tag);
		}
	} else {
		inla_error_field_is_void(__GMRFLib_FuncName, secname, prior_tag, prior->name);
	}
	return INLA_OK;
}
double priorfunc_linksnintercept(double *x, double *parameters)
{
	// input is theta, need to find the corresponding mu.
	double theta = *x;
	double mu;
	double step = 1.0e-4;
	double theta_lim = 25.0;

	// we should not be outside the limits, really...
	if (ABS(theta) < theta_lim) {
		mu = GMRFLib_cdfnorm_inv(1.0 / (1.0 + exp(-theta)));
	} else if (theta >= theta_lim) {
		mu = GMRFLib_cdfnorm_inv(1.0 / (1.0 + exp(-theta_lim)));
	} else {
		mu = GMRFLib_cdfnorm_inv(1.0 / (1.0 + exp(-(-theta_lim))));
	}

	// d_mu/d_theta = 1 / (d_theta/d_mu) 
	double deriv = 1.0 / ((inla_logit_Phi(mu + step) - inla_logit_Phi(mu - step)) / (2.0 * step));
	return (log(ABS(deriv)) + priorfunc_normal(&mu, parameters));
}

double priorfunc_pc_gevtail(double *x, double *parameters)
{
#define DIST(_xi) ((_xi)*sqrt(2.0/(1.0-(_xi))))
	double lambda = parameters[0], low = parameters[1], high = parameters[2], xi, xi_deriv, p_low, p_high, ld, d, d_deriv;

	// the internal parameterisation is in the interval [low,high]
	xi = map_interval(*x, MAP_FORWARD, (void *) &(parameters[1]));
	xi_deriv = map_interval(*x, MAP_DFORWARD, (void *) &(parameters[1]));
	d = DIST(xi);
	d_deriv = (2.0 - xi) * sqrt(2.0) / 2.0 / sqrt(1.0 / (1.0 - xi)) / SQR(1.0 - xi);
	p_low = (low > 0.0 ? ONE_mexp(-lambda * DIST(low)) : 0.0);
	p_high = (high < 1.0 ? ONE_mexp(-lambda * DIST(high)) : 1.0);
	ld = -log(p_high - p_low) + log(lambda) - lambda * d + log(ABS(d_deriv)) + log(ABS(xi_deriv));

#undef DIST
	return ld;
}

double priorfunc_pc_spde_ga(double *x, double *parameters)
{
	double theta1 = x[0], theta2 = x[1], *par = parameters, ldens, lam1, lam2;
	const int debug = 0;

	lam1 = -par[0] * log(par[1]);
	lam2 = -log(par[3]) / par[2];
	ldens = (log(lam1) - theta1 - lam1 * exp(-theta1)) + (log(lam2) + theta2 - lam2 * exp(theta2));

	if (debug) {
		fprintf(stderr, "pc_spde_ga: x = %g %g\n", x[0], x[1]);
		fprintf(stderr, "            param = %g %g %g %g\n", par[0], par[1], par[2], par[3]);
		fprintf(stderr, "            lam1 = %g  lam2 = %g  ldens = %g\n", lam1, lam2, ldens);
	}

	return ldens;
}

double priorfunc_pc_matern(double *x, double *parameters)
{
	double theta1 = x[0], theta2 = x[1], *par = parameters, ldens, lam1, lam2, dHalf;
	const int debug = 0;

	lam1 = parameters[0];
	lam2 = parameters[1];
	dHalf = parameters[2] / 2.0;

	ldens = 0.0;

	// Check if range is fixed
	if (!ISNAN(theta1))
		ldens += log(lam1 * dHalf) - dHalf * theta1 - lam1 * exp(-dHalf * theta1);

	// Check if standard deviation is fixed
	if (!ISNAN(theta2))
		ldens += log(lam2) + theta2 - lam2 * exp(theta2);

	if (debug) {
		fprintf(stderr, "pc_matern: x = %g %g\n", x[0], x[1]);
		fprintf(stderr, "           param = %g %g %g\n", par[0], par[1], par[2]);
		fprintf(stderr, "           lam1 = %g  lam2 = %g  ldens = %g\n", lam1, lam2, ldens);
	}

	return ldens;
}

double priorfunc_pc_range(double *x, double *parameters)
{
	double theta1 = x[0], ldens, lam, dHalf;
	const int debug = 0;

	lam = parameters[0];
	dHalf = parameters[1] / 2.0;

	ldens = log(lam * dHalf) - dHalf * theta1 - lam * exp(-dHalf * theta1);

	if (debug) {
		fprintf(stderr, "pc_range: x = %g\n", x[0]);
		fprintf(stderr, "          param = %g %g\n", parameters[0], parameters[1]);
		fprintf(stderr, "          lam = %g  dHalf = %g  ldens = %g\n", lam, dHalf, ldens);
	}

	return ldens;
}

double priorfunc_pc_alphaw(double *x, double *parameters)
{
	// the inla.pc.alphaw prior. alpha = exp(sc * x), where sc is given
#define _GAM (0.577215664901532860606512090082)
#define KLD(_a) ((MATHLIB_FUN(gammafn)((1.0+(_a))/(_a)) * (_a) + (_a) * log((_a)) - (_a) * _GAM + _GAM - (_a))/(_a))
#define D(_a) sqrt(2.0*DMAX(0.0, KLD((_a))))

	double d0, ld, diff, h = 1.0E-06;
	double alpha = exp(INLA_WEIBULL_ALPHA_SCALE * x[0]), lambda = parameters[0];

	d0 = D(alpha);
	if (alpha >= 1.0) {
		diff = (D(alpha + h) - d0) / h;
	} else {
		diff = (d0 - D(alpha - h)) / h;
	}

	ld = log(0.5) + lambda - lambda * d0 + log(ABS(diff)) + log(alpha * INLA_WEIBULL_ALPHA_SCALE);
#undef _GAM
#undef KLD
#undef D
	return (ld);
}

double priorfunc_pc_gamma(double *x, double *parameters)
{
	// in Maple, this is
	// asympt(log(Psi(1, n)-1/n), n,4);
	// see also testing no 48 and 49
#define SPECIAL(x) ((x > 1.0E4 ?					\
		     -2.0 * log(x) - log(2.0) + 1.0/(3.0*(x)) - 1.0/(18.0*SQR(x)) - 22.0/(405.0 * POW3(x)) : \
		     log(gsl_sf_psi_1(x) - 1.0/(x))))

	// the inla.pc.dgamma prior, which is the prior for 'a' in Gamma(1/a, 1/a) where a=0 is the base model. Here we have the
	// argument log(a). Almost the same function as priorfunc_pc_mgamma
	double ldens, d, a, a_inv, lambda;

	lambda = parameters[0];
	a = exp(x[0]);
	a_inv = 1.0 / a;
	d = sqrt(2.0 * (log(a_inv) - gsl_sf_psi(a_inv)));
	ldens = log(lambda) - lambda * d - 2.0 * log(a) - log(d) + SPECIAL(a_inv) + x[0];
#undef SPECIAL

	return ldens;
}

double priorfunc_pc_mgamma(double *x, double *parameters)
{
	// the inla.pc.dgamma prior, which is the prior for 'a' in Gamma(1/a, 1/a) where a=0 is the base model. Here we have the
	// argument log(a). this function is the pc_gamma prior for 'x' when a=exp(-x). Almost the same function as
	// priorfunc_pc_gamma
	double xx = -x[0];

	return (priorfunc_pc_gamma(&xx, parameters));
}

double priorfunc_pc_gammacount(double *x, double *parameters)
{
	// the inla.pc.dgammacount prior, which is the prior for 'a' in Gamma(a, 1) where a=1 is the base model.
	// argument is theta=log(a), so the its the density for log(a) and not a

	double lambda = parameters[0], xx, ldens, t1, t3, t4, t5, t8, t12, t14, t15, t16;

	if (ISZERO(x[0])) {
		xx = exp(INLA_REAL_SMALL);
	} else {
		xx = exp(x[0]);
	}
	t1 = log(lambda);
	t3 = gsl_sf_lngamma(xx);
	t4 = xx - 1.0;
	t5 = gsl_sf_psi(xx);
	t8 = sqrt(2.0) * sqrt(t5 * t4 - t3);
	t12 = gsl_sf_psi_1(xx);
	t14 = ABS(t12 * t4 / t8);
	t15 = log(t14);
	t16 = -t8 * lambda + t1 + t15;
	ldens = t16 - log(2.0) + log(xx);

	return ldens;
}

double priorfunc_pc_dof(double *x, double *parameters)
{
#define _NP 5
	int k;
	const int debug = 0;
	double u = parameters[0], alpha = parameters[1], lambda, dof, val, deriv;
	double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0 };
	double step, dofs[_NP], f[_NP];

	dof = map_dof(*x, MAP_FORWARD, NULL);
	lambda = -log(alpha) / inla_pcp_dof_d(u);
	// be somewhat careful evaluating the derivative
	step = sqrt(dof) * 1e-3;			       /* step-size is found empirically */
	if (dof - 2.0 * step < 2.003) {			       /* if we're to close to the lower limit */
		step = (dof - 2.003) / 2.0;
	}

	dofs[0] = dof - 2.0 * step;
	dofs[1] = dof - step;
	dofs[2] = dof;
	dofs[3] = dof + step;
	dofs[4] = dof + 2.0 * step;
	for (k = 0; k < _NP; k++) {
		f[k] = inla_pcp_dof_d(dofs[k]);
	}
	deriv = (wf[0] * f[0] + wf[1] * f[1] + wf[2] * f[2] + wf[3] * f[3] + wf[4] * f[4]) / step;
	val = log(lambda) - lambda * f[2] + log(ABS(deriv)) + log(ABS(map_dof(*x, MAP_DFORWARD, NULL)));

	if (debug) {
		P(dof);
		P(f[2]);
		P(deriv);
		P(lambda);
		P(val);
	}

	return val;
#undef _NP
}

double priorfunc_pc_sn(double *x, double *parameters)
{
	double lambda = parameters[0], val, dist, dist_max, deriv, xx, xxd, skew_max = GMRFLib_SN_SKEWMAX;

	xx = map_phi(*x, MAP_FORWARD, (void *) &skew_max);
	xxd = map_phi(*x, MAP_DFORWARD, (void *) &skew_max);
	dist = inla_pc_sn_d(xx, &deriv);
	dist_max = inla_pc_sn_d(skew_max, NULL);
	val = log(0.5) + log(lambda) - lambda * dist - LOG_1mp(exp(-lambda * dist_max)) + log(ABS(deriv)) + log(ABS(xxd));

	return val;
}

double priorfunc_pc_prec(double *x, double *parameters)
{
	double u = parameters[0], alpha = parameters[1], theta, val, xx2;
	theta = -log(alpha) / u;
	xx2 = (*x) / 2.0;
	val = log(theta / 2.0) - theta * exp(-xx2) - xx2;

	return val;
}

double priorfunc_pc_cor0(double *x, double *parameters)
{
	// alpha = Prob(|rho| > u)
	const int debug = 0;
	double u = parameters[0], alpha = parameters[1];
	double rho = map_rho(*x, MAP_FORWARD, NULL);
	double ldens, lambda, ljac, val, mu;

	mu = sqrt(-LOG_1mp(SQR(rho)));
	if (alpha <= 0 || alpha >= 1.0) {
		FIXME1("*****************************  USING lambda = u **********************************");
		lambda = u;
	} else {
		lambda = -log(alpha) / sqrt(-LOG_1mp(SQR(u)));
	}
	// add the EPS to ensure its not INFINITY...
	ldens = log(lambda) - lambda * mu + log((ABS(rho) + INLA_REAL_SMALL) / (1.0 - SQR(rho))) - log(mu + INLA_REAL_SMALL);
	ljac = log(ABS(map_rho(*x, MAP_DFORWARD, NULL)));
	val = ldens + ljac;

	if (debug) {
		fprintf(stderr, "priorfunc_pc_cor0: mu %g lambda %g ldens %g ljac %g\n", mu, lambda, ldens, ljac);
		fprintf(stderr, "priorfunc_pc_cor0: theta %g val %g\n", *x, val);
	}

	return val;
}

double priorfunc_pc_cor1(double *x, double *parameters)
{
	// alpha = Prob(rho > u)
	const int debug = 0;
	double u = parameters[0], alpha = parameters[1];
	double lambda, rho, ljac, ldens, val, mu;

	// solve for lambda
#define _Fsolve(_lam) (((ONE_mexp(-(_lam)*sqrt(1.0-u)))/(ONE_mexp(-(_lam)*M_SQRT2))) - alpha)

	int count = 0, count_max = 1000;
	double lambda_initial = -1.0, lambda_step = 1.1, h = GSL_ROOT3_DBL_EPSILON, eps_lambda = GSL_SQRT_DBL_EPSILON, df;

	if (!(u > -1.0 && u < 1.0 && alpha > sqrt((1.0 - u) / 2.0) && alpha < 1.0)) {
		char *msg;
		GMRFLib_sprintf(&msg, "Wrong cor1 prior-parameters. We must have alpha > sqrt((1-u)/2); see the documentation.");
		inla_error_general(msg);
		exit(1);
	}

	lambda = 1.0;
	if (_Fsolve(lambda) > 0.0) {
		while (_Fsolve(lambda) > 0.0) {
			lambda /= lambda_step;
			assert(count++ < count_max);
		}
	} else {
		while (_Fsolve(lambda) < 0.0 && lambda >= eps_lambda) {
			lambda *= lambda_step;
			assert(count++ < count_max);
		}
	}

	if (debug) {
		printf("priorfunc_pc_cor1: u=%g alpha=%g  initial value for lambda=%g\n", u, alpha, lambda);
	}

	count = 0;
	while (ABS(lambda - lambda_initial) > eps_lambda) {
		lambda_initial = lambda;
		df = (_Fsolve(lambda_initial + h) - _Fsolve(lambda_initial - h)) / (2.0 * h);
		lambda = lambda_initial - _Fsolve(lambda) / df;
		if (debug) {
			printf("priorfunc_pc_cor1: iteration=%d lambda=%g\n", count, lambda);
		}
		assert(count++ < count_max);
	}
	if (debug)
		printf("priorfunc_pc_cor1: function value %g\n", _Fsolve(lambda));

#undef Fsolve
	rho = map_rho(*x, MAP_FORWARD, NULL);
	mu = sqrt(1.0 - rho);
	ldens = log(lambda) - lambda * mu - LOG_1mp(exp(-lambda * M_SQRT2)) - log(2.0 * mu);
	ljac = log(ABS(map_rho(*x, MAP_DFORWARD, NULL)));
	val = ldens + ljac;

#undef _Fsolve
	return (val);
}

double priorfunc_jeffreys_df_student_t(double *x, double *UNUSED(parameters))
{
	double df = map_dof(x[0], MAP_FORWARD, NULL);
	double value, log_jacobian;

	if (1) {
#define _DIGAMMA(xx) gsl_sf_psi(xx)
#define _TRIGAMMA(xx) gsl_sf_psi_1(xx)

		value =
		    0.5 * log(df / (df + 3.0)) + 0.5 * log(_TRIGAMMA(df / 2.0) - _TRIGAMMA((df + 1.0) / 2.0) -
							   2.0 * (df + 3.0) / (df * SQR(df + 1.0)))
		    - log(0.7715233664);		       /* normalising constant: computed in R from 2 to infinity */

		log_jacobian = x[0];
		return value + log_jacobian;
#undef _DIGAMMA
#undef _TRIGAMMA
	} else {
		value = -log(df);
		log_jacobian = x[0];
		return value + log_jacobian;
	}
}

double priorfunc_bymjoint(double *logprec_besag, double *p_besag, double *logprec_iid, double *p_iid)
{
	/*
	 * This is the Jon Wakefield prior. Notation: U Spatial, V unstruct.  Conceptually, there is a prior for total variance (U + V) and a prior for the
	 * proportion of the total that is spatial. The spatial variance is approximated as the mean of the marginal spatial variances in an ICAR model.  This mean
	 * is a function of the (known) neighborhood structure and the conditional variance parameter.  Thus, the joint distribution is expressed of the conditional
	 * spatial variance parameter and the non-spatial parameter. The value returned, is the log of the joint for the log-precisions, (-\log Conditional.Var(U),
	 * -\log Var(V))
	 */

	double a, b, c, d, var_u, var_v, val, mean_ievalue;

	a = p_besag[0];					       /* Parameters for the IGamma(a,b) */
	b = p_besag[1];
	mean_ievalue = p_besag[2];			       /* mean(1/eigenvalues) for the reference precision matrix */
	c = p_iid[0];					       /* Parameters for the Beta(c,d) */
	d = p_iid[1];

	var_u = 1.0 / exp(*logprec_besag);		       /* var_u = Conditional.Var(U) i.e. Cond.var_u_i = var_u/n_i */
	var_v = 1.0 / exp(*logprec_iid);

	/*
	 * the log prior for Conditonal.Var(U) and Var(V) 
	 */
	val = a * log(b) + gsl_sf_lngamma(c + d) - gsl_sf_lngamma(a) - gsl_sf_lngamma(c) - gsl_sf_lngamma(d)
	    - (a + c + d) * log((var_u * mean_ievalue) + var_v) + (c - 1.0) * log(var_u * mean_ievalue)
	    + (d - 1.0) * log(var_v) - b / ((var_u * mean_ievalue) + var_v) + log(mean_ievalue);

	/*
	 * the Jacobian for converting to (-\log Conditional.Var(U), -\log Var(V)) 
	 */
	val += log(var_u) + log(var_v);

	return val;
}

double priorfunc_invalid(double *UNUSED(x), double *UNUSED(parameters))
{
	inla_error_general("Prior 'invalid' is used, but it is not ment to be.");
	exit(EXIT_FAILURE);
	return 0.0;
}

double priorfunc_betacorrelation(double *x, double *parameters)
{
	/*
	 * The prior for the correlation coefficient \rho is Beta(a,b), scaled so that it is defined on (-1,1)
	 * The function returns the log prior for \rho.intern = log((1 +\rho)/(1-\rho))
	 */
	double p = exp(*x) / (1.0 + exp(*x)), a = parameters[0], b = parameters[1];
	return log(gsl_ran_beta_pdf(p, a, b)) + (*x) - 2.0 * log1p(exp(*x));
}

double priorfunc_logitbeta(double *x, double *parameters)
{
	/*
	 * The prior for the the logit of a Beta(a,b), logit(p) = log(p/(1-p))
	 */
	double p = exp(*x) / (1.0 + exp(*x)), a = parameters[0], b = parameters[1];
	return log(gsl_ran_beta_pdf(p, a, b)) + (*x) - 2.0 * log1p(exp(*x));
}

double priorfunc_logflat(double *x, double *UNUSED(parameters))
{
	return exp(*x);
}

double priorfunc_logiflat(double *x, double *UNUSED(parameters))
{
	return exp(-*x);
}

double priorfunc_flat(double *UNUSED(x), double *UNUSED(parameters))
{
	return 0.0;
}

double priorfunc_minuslogsqrtruncnormal(double *x, double *parameters)
{
	/*
	 * Requested feature from Sophie Ancelet.
	 * 
	 * This is the prior for \sigma ~ TrucNormal(mean, 1/prior_prec), where \sigma > 0, and the internal parameter is \log\tau = \log(1/\sigma^2) = -\log
	 * \sigma^2, which explans the name.
	 */
	double sd = exp(-0.5 * (*x)), val;

	val = priorfunc_normal(&sd, parameters) - log(gsl_cdf_gaussian_Q(-parameters[0], 1.0 / sqrt(parameters[1]))) +
	    // log(Jacobian)
	    log(ABS(-0.5 * sd));

	return val;
}

double priorfunc_pc_ar(double *x, double *parameters)
{
	int i, p;
	double lambda, *b, *gamma, *pacf, ldens, logjac;

	p = (int) parameters[1];
	lambda = parameters[0];
	b = Calloc(3 * p, double);
	gamma = &(b[p]);
	pacf = &(b[2 * p]);

	for (i = 0, logjac = 0.0; i < p; i++) {
		b[i] = 0.5;
		// x is internal and this gives us the pacf. 
		pacf[i] = ar_map_pacf(x[i], MAP_FORWARD, NULL);
		// but the pc-simplex prior is given in terms of 'gamma'
		gamma[i] = -LOG_1mp(SQR(pacf[i]));
		// hence we need two jacobians, one for x->pacf and one for pacf->gamma. recall that we have a singularity for
		// x[i]=0
		double xtmp = (ISZERO(pacf[i]) ? INLA_REAL_SMALL : pacf[i]);
		logjac += log(ABS(ar_map_pacf(x[i], MAP_DFORWARD, NULL))) + log(ABS(xtmp / (1.0 - SQR(pacf[i]))));
	}
	ldens = inla_pc_simplex_d(gamma, b, p, lambda) + logjac;

	Free(b);
	return (ldens);
}

double priorfunc_ref_ar(double *x, double *parameters)
{
	int i, p;
	double pacf[3], ldens;

	p = (int) parameters[0];
	assert(p >= 0 && p <= 3);
	for (i = 0; i < p; i++) {
		pacf[i] = map_phi(x[i], MAP_FORWARD, NULL);
	}

	switch (p) {
	case 0:
	{
		ldens = 0.0;
	}
		break;
	case 1:
	{
		ldens = -log(M_PI)
		    - 0.5 * LOG_1mp(SQR(pacf[0]))
		    + log(ABS(map_phi(x[0], MAP_DFORWARD, NULL)));
	}
		break;
	case 2:
	{
		ldens = -2.0 * log(M_PI)
		    - 0.5 * LOG_1mp(SQR(pacf[0]))
		    - 0.5 * LOG_1mp(SQR(pacf[1]))
		    + log(ABS(map_phi(x[0], MAP_DFORWARD, NULL)))
		    + log(ABS(map_phi(x[1], MAP_DFORWARD, NULL)));
	}
		break;
	case 3:
	{
		ldens = -2.0 * log(M_PI)
		    - log(1.12)
		    - 0.5 * LOG_1mp(SQR(pacf[0]))
		    - 0.5 * LOG_1mp(SQR(pacf[1]))
		    - 0.5 * LOG_1mp(SQR(pacf[2]))
		    - M_PI * SQR(pacf[2])
		    + log(ABS(map_phi(x[0], MAP_DFORWARD, NULL)))
		    + log(ABS(map_phi(x[1], MAP_DFORWARD, NULL)))
		    + log(ABS(map_phi(x[2], MAP_DFORWARD, NULL)));
	}
		break;
	default:
		ldens = NAN;
		GMRFLib_ASSERT(p <= 3 && p >= 0, GMRFLib_EPARAMETER);
	}

	return (ldens);
}

double priorfunc_laplace(double *x, double *parameters)
{
	double mean = parameters[0], lambda = sqrt(2.0 * parameters[1]);
	return (-M_LN2 + log(lambda) - lambda * ABS(*x - mean));
}

double priorfunc_beta(double *x, double *parameters)
{
	double xx = *x, a = parameters[0], b = parameters[1];

	return log(gsl_ran_beta_pdf(xx, a, b));
}

double priorfunc_loggamma(double *x, double *parameters)
{
	/*
	 * return log(loggamma(x,a,b)). NOTE: if y ~ gamma(a,b), then log(y) ~ loggamma(a,b). 
	 */
	double val = exp(*x);
	return priorfunc_gamma(&val, parameters) + (*x);
}

double priorfunc_dirichlet(double *x, double *parameters)
{
#define _F_logit(_x) (1.0/(1.0+exp(-(_x))))
#define _f_logit(_x) (exp(-(_x)) / SQR(1.0 + exp(-(_x))))

#define _F_probit(_x) inla_Phi(_x)
#define _f_probit(_x) (0.39894228040143270286 * exp(-0.5 * SQR(_x)))

#define _F(_x) (cdf_logit ? _F_logit(_x) : _F_probit(_x))
#define _f(_x) (cdf_logit ? _f_logit(_x) : _f_probit(_x))

	double alpha = parameters[0], nclasses = parameters[1], ld;
	int cdf_logit = ((inla_pom_cdf_tp) parameters[2] == POM_CDF_LOGIT);
	int K = (int) nclasses, k;
	const int debug = 0;
	double *work = Calloc(4 * K, double);
	double *xx = work, *alphas = work + K, *qs = work + 2 * K, *v = work + 3 * K;

	// from the internal representation, x, to cutpoints, xx 
	xx[0] = x[0];
	for (k = 1; k < K - 1; k++) {
		xx[k] = xx[k - 1] + exp(x[k]);
	}
	// from cutpoints, xx, to quantiles, qs

	for (k = 0; k < K - 1; k++) {
		qs[k] = _F(xx[k]);
	}

	// from quantiles, qs, to Dirichlet variables, v
	v[0] = qs[0];
	for (k = 1; k < K - 1; k++) {
		v[k] = qs[k] - qs[k - 1];
	}
	v[K - 1] = 1.0 - qs[K - 2];

	for (k = 0; k < K; k++) {			       /* also to provide this */
		alphas[k] = alpha;
	}
	ld = gsl_ran_dirichlet_lnpdf((size_t) K, alphas, v);

	// finally, add the log-jacobian of the transformation. This is approximately correct, as there are K variables and a
	// sum to 1 constr for the Dirichlet, but K-1 variables for the theta without the constr. Maybe look into this later

	ld += log(_f(xx[0]));
	for (k = 1; k < K - 1; k++) {
		ld += log(_f(xx[k])) + x[k];
	}

	if (debug) {
		printf("priorfunc_dirichlet: alpha = %g  K=%1d\n", alpha, K);
		for (int i = 0; i < K - 1; i++) {
			printf("input theta[%1d] = %g\n", i, x[i]);
		}
		for (int i = 0; i < K; i++) {
			printf("input v[%1d] = %g\n", i, v[i]);
		}
		P(ld);
	}

	Free(work);
#undef _F
#undef _f
#undef _F_logit
#undef _f_logit
#undef _F_probit
#undef _f_probit

	return (ld);
}

double priorfunc_gamma(double *x, double *parameters)
{
	/*
	 * return log(gamma(x, a, b))
	 *
	 * parametes (a,b) are such that E(x) = a/b and Var(x) = a/b^2
	 */
	double a = parameters[0], b = 1.0 / parameters[1];

	if (0) {
		/*
		 * this crash of'course for extreme values 
		 */
		return log(gsl_ran_gamma_pdf(*x, a, b));
	} else {
		/*
		 * while this is ok. (code copied from GSL...) 
		 */
		if (*x < 0.0) {
			assert(*x >= 0.0);
		} else if (*x == 0.0) {
			if (a == 1.0) {
				return 1.0 / b;
			} else {
				GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			}
		} else if (a == 1.0) {
			return -(*x) / b - log(b);
		} else {
			return (a - 1.0) * log((*x) / b) - (*x) / b - gsl_sf_lngamma(a) - log(b);
		}
	}
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	return 0.0;
}

double priorfunc_gaussian(double *x, double *parameters)
{
	return priorfunc_normal(x, parameters);
}

double priorfunc_normal(double *x, double *parameters)
{
	/*
	 * return log(normal(x, mean, precision)) 
	 */
	if (ISZERO(parameters[1])) {
		return 0.0;				       // = log(1)
	} else {
		double mean = parameters[0], prec = parameters[1];
		return (-0.9189385332 + 0.5 * log(prec) - 0.5 * prec * SQR(*x - mean));
	}
}

double priorfunc_mvnorm(double *x, double *parameters)
{
	/*
	 * this is the multivariate normal 
	 */
	int n = (int) parameters[0], i, j;

	if (n == 0) {
		return 0.0;
	}

	double *mean, *Q, *chol, *xx, q = 0.0, logdet = 0.0;

	mean = &(parameters[1]);
	Q = &(parameters[1 + n]);
	chol = NULL;
	xx = Calloc(n, double);

	if (0) {
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%g ", Q[i + j * n]);
			printf("\n");
		}
		for (i = 0; i < n; i++)
			printf("%g\n", mean[i]);
	}

	GMRFLib_comp_chol_general(&chol, Q, n, &logdet, 0);
	for (i = 0; i < n; i++) {
		xx[i] = x[i] - mean[i];
	}

	/*
	 * q = xx^T * Q * xx. I dont have any easy function for matrix vector except the messy BLAS-FORTRAN-INTERFACE.... so I just do this manually now. FIXME
	 * later. 
	 */
	q = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			q += xx[i] * Q[i + n * j] * xx[j];
		}
	}

	Free(xx);
	Free(chol);

	return (-n / 2.0 * log(2 * M_PI) + 0.5 * logdet - 0.5 * q);
}

int inla_iid_wishart_adjust(int dim, double *theta)
{
	/*
	 * adjust rho's in theta with factor f until the matrix is SPD. 
	 */
#define _IDX(_i, _j) ((_i) + (_j)*(dim))

	int i, j, k, ok = 0;
	const int debug = 0;
	int n_theta = inla_iid_wishart_nparam(dim);
	double f = 0.95, *S = NULL, *chol = NULL;

	S = Calloc(ISQR(dim), double);
	while (!ok) {
		k = 0;
		for (i = 0; i < dim; i++) {
			S[_IDX(i, i)] = 1.0 / theta[k];
			k++;
		}
		for (i = 0; i < dim; i++) {
			for (j = i + 1; j < dim; j++) {
				S[_IDX(i, j)] = S[_IDX(j, i)] = theta[k] / sqrt(theta[i] * theta[j]);
				k++;
			}
		}
		assert(k == n_theta);

		if (debug) {
			printf("\n");
			for (i = 0; i < dim; i++) {
				for (j = 0; j < dim; j++) {
					printf(" %.12f", S[_IDX(i, j)]);
				}
				printf("\n");
			}
			printf("\n");
		}

		if (GMRFLib_comp_chol_general(&chol, S, dim, NULL, !GMRFLib_SUCCESS) != GMRFLib_SUCCESS) {
			/*
			 * only adjust the rho's 
			 */
			if (debug) {
				printf("matrix is not spd, adjust with factor %f\n", f);
			}
			for (i = dim; i < n_theta; i++) {
				theta[i] *= f;
			}
		} else {
			ok = 1;
		}
	}
	Free(S);
	Free(chol);

#undef _IDX
	return (ok ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

double priorfunc_wishart1d(double *x, double *parameters)
{
	// its the same and quicker

	// return priorfunc_wishart_generic(1, x, parameters);

	double p[2];
	p[0] = parameters[0] / 2.0;
	p[1] = parameters[1] / 2.0;

	return priorfunc_gamma(x, p);
}

double priorfunc_wishart2d(double *x, double *parameters)
{
	return priorfunc_wishart_generic(2, x, parameters);
}

double priorfunc_wishart3d(double *x, double *parameters)
{
	return priorfunc_wishart_generic(3, x, parameters);
}

double priorfunc_wishart4d(double *x, double *parameters)
{
	return priorfunc_wishart_generic(4, x, parameters);
}

double priorfunc_wishart5d(double *x, double *parameters)
{
	return priorfunc_wishart_generic(5, x, parameters);
}

double priorfunc_wishart_generic(int idim, double *x, double *parameters)
{
	/*
	 * 
	 * Q ~ Wishart(r, R^{-1} )
	 * 
	 * input is x = [ tau0, tau1, tau2, rho01, rho02, rho12 ] which parameterise Q.
	 * 
	 * prior parameters are p = [ r, R00, R11, R22, R12, R13, R23 ]
	 * 
	 * output is the logdensity for x!!!! (BE AWARE)
	 */
	gsl_matrix *R = NULL, *Q = NULL;
	double r, val;
	const int debug = 0;
	int fail;
	size_t i, ii, j, k, dim = (size_t) idim;

	size_t n_x = (size_t) inla_iid_wishart_nparam(idim);
	size_t n_param = n_x + 1;

	if (debug) {
		for (i = 0; i < n_param; i++) {
			printf("parameters[%d] = %g\n", (int) i, parameters[i]);
		}
		for (i = 0; i < n_x; i++) {
			printf("x[%d] = %g\n", (int) i, x[i]);
		}
	}

	r = parameters[0];
	R = gsl_matrix_calloc(dim, dim);
	Q = gsl_matrix_calloc(dim, dim);

	fail = inla_iid_wishart_adjust(idim, x);

	/*
	 * offset of 1, since parameters[0] = r
	 */
	k = 1;
	for (i = 0; i < dim; i++) {
		gsl_matrix_set(R, i, i, parameters[k]);
		k++;
	}
	for (i = 0; i < dim; i++) {
		for (j = i + 1; j < dim; j++) {
			gsl_matrix_set(R, i, j, parameters[k]);
			gsl_matrix_set(R, j, i, parameters[k]);
			k++;
		}
	}
	assert(k == n_param);

#define COMPUTE_Q(Q_)							\
	if (1) {							\
		k = 0;							\
		for(i = 0; i<dim; i++) {				\
			gsl_matrix_set(Q_, i, i, 1/x[k]);		\
			k++;						\
		}							\
		for(i=0; i<dim; i++) {					\
			for(j=i+1; j < dim; j++) {			\
				double value = x[k] / sqrt(x[i] * x[j]); \
				gsl_matrix_set(Q_, i, j, value);	\
				gsl_matrix_set(Q_, j, i, value);	\
				k++;					\
			}						\
		}							\
		assert(k == n_x);					\
		if (debug) printf("Covmatrix:\n");			\
		if (debug) GMRFLib_printf_gsl_matrix(stdout, Q_, NULL); \
		if (debug) printf("Precision:\n");			\
		if (debug) GMRFLib_printf_gsl_matrix(stdout, Q_, NULL); \
		GMRFLib_gsl_spd_inverse(Q_);				\
	}

	COMPUTE_Q(Q);
	val = GMRFLib_Wishart_logdens(Q, r, R) + (fail ? PENALTY : 0.0);

	/*
	 * tau1 and tau1 are the *MARGINAL* precisions, which it should be.  The jacobian is computed like this:
	 * 
	 * with(LinearAlgebra); with(Student[VectorCalculus]); S := matrix(3,3, [ 1/tau0, rho01/sqrt(tau0*tau1), rho02/sqrt(tau0*tau2), rho01/sqrt(tau0*tau1),
	 * 1/tau1, rho12/sqrt(tau1*tau2), rho02/sqrt(tau0*tau2), rho12/sqrt(tau1*tau2), 1/tau2]); Q := inverse(S); simplify(Determinant(Jacobian([Q[1,1], Q[2,2],
	 * Q[3,3], Q[1,2], Q[1,3],Q[2,3]], [tau0,tau1,tau2,rho01,rho02,rho12])));
	 * 
	 * this gives a very long answer of'course; so we have to do this numerically 
	 */
	gsl_matrix *QQ = NULL, *J = NULL;
	double f, save;

	/*
	 * for the numerical derivatives: compute the `population' variance: Det(Sigma), and set f = (Det(Sigma))^1/dim. 
	 */
	f = 1.0e-6 * pow(exp(-GMRFLib_gsl_spd_logdet(Q)), 1.0 / (double) idim);	/* Yes, its a minus... */
	QQ = GMRFLib_gsl_duplicate_matrix(Q);
	J = gsl_matrix_calloc(n_x, n_x);

	/*
	 * the precision terms *can* get negative since we're using a central estimate, so we need to check the step-size 
	 */
	for (ii = 0; ii < dim; ii++) {
		f = DMIN(f, 0.5 * x[ii]);
	}

	for (ii = 0; ii < n_x; ii++) {
		save = x[ii];

		x[ii] += f;
		COMPUTE_Q(Q);
		x[ii] = save;

		x[ii] -= f;
		COMPUTE_Q(QQ);
		x[ii] = save;

		k = 0;
		for (i = 0; i < dim; i++) {
			gsl_matrix_set(J, ii, k, (gsl_matrix_get(Q, k, k) - gsl_matrix_get(QQ, k, k)) / (2.0 * f));
			assert(k == i);
			k++;
		}
		for (i = 0; i < dim; i++) {
			for (j = i + 1; j < dim; j++) {
				gsl_matrix_set(J, ii, k, (gsl_matrix_get(Q, i, j) - gsl_matrix_get(QQ, i, j)) / (2.0 * f));
				k++;
			}
		}
	}
	assert(k == n_x);

	gsl_permutation *p;
	int signum;
	double logdet;

	p = gsl_permutation_alloc(n_x);
	gsl_linalg_LU_decomp(J, p, &signum);
	logdet = gsl_linalg_LU_lndet(J);		       /* log(abs|J|) */

	if (debug) {
		P(logdet);
	}

	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_matrix_free(QQ);
	gsl_matrix_free(J);
	gsl_permutation_free(p);

	val += logdet;

#undef COMPUTE_Q

	return val;
}

double priorfunc_wishartk_2d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(2, x, parameters);
}

double priorfunc_wishartk_3d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(3, x, parameters);
}

double priorfunc_wishartk_4d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(4, x, parameters);
}

double priorfunc_wishartk_5d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(5, x, parameters);
}

double priorfunc_wishartk_6d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(6, x, parameters);
}

double priorfunc_wishartk_7d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(7, x, parameters);
}

double priorfunc_wishartk_8d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(8, x, parameters);
}

double priorfunc_wishartk_9d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(9, x, parameters);
}

double priorfunc_wishartk_10d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(10, x, parameters);
}
double priorfunc_wishartk_11d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(11, x, parameters);
}
double priorfunc_wishartk_12d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(12, x, parameters);
}
double priorfunc_wishartk_13d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(13, x, parameters);
}
double priorfunc_wishartk_14d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(14, x, parameters);
}
double priorfunc_wishartk_15d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(15, x, parameters);
}
double priorfunc_wishartk_16d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(16, x, parameters);
}
double priorfunc_wishartk_17d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(17, x, parameters);
}
double priorfunc_wishartk_18d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(18, x, parameters);
}
double priorfunc_wishartk_19d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(19, x, parameters);
}
double priorfunc_wishartk_20d(double *x, double *parameters)
{
	return priorfunc_wishartk_generic(20, x, parameters);
}

double priorfunc_wishartk_generic(int idim, double *x, double *parameters)
{
	/*
	 * 
	 * Q ~ Wishart(r, R^{-1} )
	 * 
	 * output is the logdensity for x!!!! 
	 */
	gsl_matrix *R = NULL, *Q = NULL, *QQ, *L = NULL;
	double r, val;
	const int debug = 0;
	size_t i, ii, j, k, dim = (size_t) idim;

	size_t n_x = (size_t) INLA_WISHARTK_NTHETA(idim);
	size_t n_param = n_x + 1;

	if (debug) {
		for (i = 0; i < n_param; i++) {
			printf("parameters[%d] = %g\n", (int) i, parameters[i]);
		}
		for (i = 0; i < n_x; i++) {
			printf("x[%d] = %g\n", (int) i, x[i]);
		}
	}

	r = parameters[0];
	R = gsl_matrix_calloc(dim, dim);
	L = gsl_matrix_calloc(dim, dim);
	Q = gsl_matrix_calloc(dim, dim);
	QQ = gsl_matrix_calloc(dim, dim);

	/*
	 * offset of 1, since parameters[0] = r
	 */
	k = 1;
	for (i = 0; i < dim; i++) {
		gsl_matrix_set(R, i, i, parameters[k]);
		k++;
	}
	for (i = 0; i < dim; i++) {
		for (j = i + 1; j < dim; j++) {
			gsl_matrix_set(R, i, j, parameters[k]);
			gsl_matrix_set(R, j, i, parameters[k]);
			k++;
		}
	}
	assert(k == n_param);

	inla_wishartk_build_Q(idim, x, Q, L);
	val = GMRFLib_Wishart_logdens(Q, r, R);

	gsl_matrix *J = NULL;
	double h, save, *xx = NULL;

	xx = Calloc(n_x, double);
	Memcpy(xx, x, n_x * sizeof(double));

	h = 0.005;
	J = gsl_matrix_calloc(n_x, n_x);
	// do not need for h=0 as the weight is zero
	double wf[] = { 1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0 };
	double hh[] = { -2.0 * h, -h, h, 2.0 * h };

	gsl_matrix_set_zero(QQ);
	for (ii = 0; ii < n_x; ii++) {

		for (size_t ih = 0; ih < sizeof(hh) / sizeof(double); ih++) {
			save = xx[ii];
			xx[ii] += hh[ih] * h;
			inla_wishartk_build_Q(idim, xx, Q, L);
			xx[ii] = save;
			for (i = 0; i < dim; i++) {
				for (j = 0; j < dim; j++) {
					gsl_matrix_set(QQ, i, j, gsl_matrix_get(QQ, i, j) + gsl_matrix_get(Q, i, j) * wf[ih]);
				}
			}
		}

		k = 0;
		for (i = 0; i < dim; i++) {
			gsl_matrix_set(J, ii, k, gsl_matrix_get(QQ, k, k) / SQR(h));
			k++;
		}
		for (j = 0; j < dim; j++) {
			for (i = j + 1; i < dim; i++) {
				gsl_matrix_set(J, ii, k, gsl_matrix_get(QQ, i, j) / SQR(h));
				k++;
			}
		}
		assert(k == n_x);
	}

	gsl_permutation *p = NULL;
	int signum;
	double logdet;

	p = gsl_permutation_alloc(n_x);
	gsl_linalg_LU_decomp(J, p, &signum);
	logdet = gsl_linalg_LU_lndet(J);		       /* log(abs|J|) */
	val += logdet;

	if (debug) {
		P(logdet);
	}

	Free(xx);
	gsl_matrix_free(R);
	gsl_matrix_free(Q);
	gsl_matrix_free(QQ);
	gsl_matrix_free(L);
	gsl_matrix_free(J);
	gsl_permutation_free(p);

	return val;
}
