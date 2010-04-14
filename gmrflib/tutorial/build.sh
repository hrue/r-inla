#! /bin/bash -eu
sources="
	block_update
	check_matrix
	fixed_condition
	graph_from_int_vec
	ln_gamma
	ln_gamma_density
	ln_gauss_density
	ln_poisson_density
	linear_condition
	main
	multiply_matrix
	sim_poisson
	stochastic_condition
	uniform_precision
	unconditional
"
objects=""
main_ok="yes"
for name in $sources
do
	if [ -e $name.o ]
	then
		echo "Using existing $name.o"
	else
		main_ok="no"
		echo "Creating $name.o"
		gcc -std=gnu99 \
			-Wall \
			-g \
			-march=native \
			-mfpmath=sse \
			-msse4 \
			-funroll-loops \
			-ftracer \
			-fopenmp \
			-pipe \
			-I $HOME/prefix/gmrf/include \
		 	-c -o $name.o \
			$name.c
	fi
	objects="$objects $name.o"
done
list="
	gen_svd
	pos_inv
"
for name in $list
do
	if [ -e $name.o ]
	then
		echo "Using existing $name.o"
	else
		main_ok="no"
		echo "Creating $name.o"
		gfortran \
			-Wall \
			-g \
			-march=native \
			-mfpmath=sse \
			-msse4 \
			-funroll-loops \
			-ftracer \
			-fopenmp \
			-pipe \
			-I $HOME/prefix/gmrf/include \
	 		-c -o $name.o \
			$name.f
	fi
	objects="$objects $name.o"
done
#
if [ "$main_ok" != "yes" ]
then
	if [ -e main ]
	then
		rm  main
	fi
	echo "Linking main"
	gfortran \
		-Wall \
		-g \
		-march=native \
		-mfpmath=sse \
		-msse4 \
		-funroll-loops \
		-ftracer \
		-fopenmp \
		-pipe \
		-o main \
		$objects \
		-L/home/bradbell/prefix/gmrf/lib \
		-lGMRFLib \
		-lgsl -ltaucs -lmetis -llapack -lblas -lgsl -lgslcblas -lz -lm
fi
#
echo "doxygen >& doxygen.log"
doxygen >& doxygen.log
#
if grep 'Warning:' doxygen.log > /dev/null
then
	echo "There are warnings in doxygen.log"
	exit 1
fi 
exit 0
