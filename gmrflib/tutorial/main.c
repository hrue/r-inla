// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.
/*!
\file main.c
\brief Run all of the examples / tests in the GMRFLib tutorial.
*/

/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: main.c,v 1.20 2010/03/21 14:21:30 bradbell Exp $"; 

# include <stdio.h>

extern int block_update            (void);
extern int fixed_condition         (void);
extern int graph_from_int_vec_test (void);
extern int linear_condition        (void);
extern int ln_gamma_test           (void);
extern int ln_gamma_density_test   (void);
extern int ln_poisson_density_test (void);
extern int sim_poisson_test        (void);
extern int stochastic_condition    (void);
extern int unconditional           (void);

/*!
Run one of the examples / tests in the GMRFLib tutorial.

\param[in] name
identifies (for the user) which test is being run.

\param test
is the function that preforms the test.

\return
returns true (1) if the test passes
and false (0) if it fails.
*/
static int run(const char *name, int test (void) )
{	int ok = 1;
	if( test() )
		printf("OK:     %s\n", name);
	else
	{	printf("Error:  %s\n", name);
		ok = 0;
	}
	return ok;
}
/*!
Run all of the examples / tests in the GMRFLib tutorial.
*/
int main(void)
{	int ok = 1;
	ok    &= run("block_update",              block_update);
	ok    &= run("fixed_condition",           fixed_condition);
	ok    &= run("graph_from_int_vec_test",   graph_from_int_vec_test);
	ok    &= run("linear_condition",          linear_condition);
	ok    &= run("ln_gamma_test",             ln_gamma_test);
	ok    &= run("ln_gamma_density_test",     ln_gamma_density_test);
	ok    &= run("ln_poisson_density_test",   ln_poisson_density_test);
	ok    &= run("sim_poisson_test",          sim_poisson_test);
	ok    &= run("stochastic_condition",      stochastic_condition);
	ok    &= run("unconditional",             unconditional);

	return ! ok;
}
