// GMRFLib-Tutorial: Copyright (C) 2010 Bradley M. Bell
// Distributed under the terms of the GNU Public License version 2 or higher.

/*!
\file graph_from_int_vec.c
\brief Create a graph from a single integer vector.
*/
/// identify who made the last change to this file and when
static const char RCSId[] = "$Id: graph_from_int_vec.c,v 1.5 2010/03/17 02:52:53 bradbell Exp $";

# include "GMRFLib/GMRFLib.h"

/*!
Create a graph structure from neighbour information in an integer vector. 

\param[in] n
is the number of nodes in the graph.

\param[in] len
is the number of integers in the vector \c vec.

\param[in] vec
is the vector containting the neighbour list for each node.
To be specific, it has the following structure:
\verbatim
int vec[] = {
0,     nbs[0][0],   nbs[0][1], ... ,     nbs[0][ nnbs[0]-1 ], -1
1,     nbs[1][0],   nbs[1][1], ... ,     nbs[1][ nnbs[1]-1 ], -1
 .              .                  .                        .
 .              .                  .                        .
 .              .                  .                        .
n-1, nbs[n-1][0], nbs[n-1][1], ... , nbs[n-1][ nnbs[n-1]-1 ], -1
}
\endverbatim
where 
<tt>nnbs[i]</tt> is the number of neighbours for node \c i and
<tt>nbs[i][j]</tt> is the \c j-th neighbour for node \c i.

\par Assumptions.
For <tt>i = 0 , ... , n-1</tt> and for <tt>j = 0 , ... , nnbs[i]-1</tt>,
we have <tt>nnbs[i][j] != i</tt> and <tt>0 <= nnbs[i][j] < n</tt>. 
In addition, \c len is the total number integers in the diagram above; 
i.e., <tt>vec[len-1]</tt> is the final <tt>-1</tt> in the diagram above.

\returns
a pointer to the graph data structure.
When it is no longer needed, this should be freed using
\code
	GMRFLib_free_graph(graph)
\endcode
where \c graph is the value retured by \c graph_from_int_vec.
*/

GMRFLib_graph_tp*  graph_from_int_vec(
	int         n        ,
	int         len      ,
	int*        vec      )
{	int  sum_nnbs;
	int  vec_index;
	int  nbs_index;
	int* nbs_all;
	int i, j;
	GMRFLib_graph_tp* graph;

	// initialize the graph data structure
	GMRFLib_make_empty_graph(&graph);
	graph->n    = n;
	graph->nnbs = Calloc(n+1, int);
	graph->nbs  = Calloc(n+1, int*);

	// Check integer vector data structure and
	// determine the total number of neighbours in the graph
	vec_index = 0;
	sum_nnbs  = 0;
	for(i = 0; i < n; i++)
	{	// ensure vec is long enough to contain the value i
		// and terminator for this list of neighbours
		assert( vec_index + 1 < len );
	
		// check that this is the correct neighbour list
		assert( vec[vec_index] == i );

		// initialize number of neighbours for node i
		graph->nnbs[i] = 0;
		while( vec[++vec_index] >= 0 )
		{	// check vec is long enough for this neighbour and
			// terminator for this list of neighbours
			assert( vec_index + 1 < len  );

			// check that this neighbour has a valid value
			assert( vec[vec_index] < n  );

			// check that it is not the same as node i
			assert( vec[vec_index] != i );

			// increment the count for this neighbour list
			++(graph->nnbs[i]);
		}
		// check the terminator for this list of neighbours
		assert( vec[vec_index] == -1 );

		// add number of neighbours for this node to total sum
		sum_nnbs      += graph->nnbs[i];

		// index where data for next list of neighbours starts;
		++vec_index;
	}
	assert( vec_index == len );

	nbs_all    = Calloc(sum_nnbs, int);
	nbs_index  = 0;
	for (i = 0; i < n; i++)
	{	if( graph->nnbs[i] == 0 )
			graph->nbs[i] = NULL;
		else
		{	// pointer in graph where this neighbour list begins
			graph->nbs[i] = nbs_all + nbs_index; 
			// index in vec where this neighbour list begins
			vec_index    = nbs_index + 2 * i;
			// copy this neighbour list 
			for(j = 0; j < graph->nnbs[i]; j++)
				graph->nbs[i][j] = vec[vec_index + 1 + j];
			nbs_index   += graph->nnbs[i];
		}
	}

	// prepare the graph for computations
	GMRFLib_prepare_graph(graph);

	return graph;
}

/*!
Reads an unsigned integer from a data file.

This is a \c static function and should only be used by
\ref graph_from_int_vec_test .

\param fp
is a pointer to the file (which has been opened with \c fopen)
the file name.

\returns
the value of the integer.
*/

static int read_number_from_file(FILE *fp)
{	char number[20];
	int i, ch;

	// skip to start of next number
	ch = fgetc(fp);
	while( ch < '0' || '9' < ch )
	{	assert( ch != EOF );
		ch = fgetc(fp);
	}
	// store ascii in a buffer
	i = 0;
	number[i++] = ch;
	ch          = fgetc(fp);
	while( '0' <= ch && ch <= '9' )
	{	assert( i < 19 );
		number[i++] = ch;
		ch          = fgetc(fp);
	}
	number[i] = '\0';

	return atoi(number);
}

/*!
Example and test of \ref graph_from_int_vec .

\returns
one if the test passes and zero otherwise.
*/
	
int graph_from_int_vec_test(void)
{	int               ok, nnbs, index, i, j;
	FILE             *fp;	
	GMRFLib_graph_tp *graph;
	int               n      = 4;
	// for node equal 0 to n-1:
	// node, followed by list of neighbours, followed by -1 
	int vec[] = {
	   0,    1,          -1,
	   1,    0,     2,   -1,
	   2,    1,     3,   -1,
	   3,    2,          -1
	};
	// length of the vector
	int len = sizeof(vec) / sizeof(vec[0]);

	// create a graph data structure from the vector
	graph = graph_from_int_vec(n, len, vec);

	// write the graph to a file
	fp    = fopen("graph_from_int_vec.out", "w");
	GMRFLib_print_graph(fp, graph);
	fclose(fp);

	// open the file for reading
	fp    = fopen("graph_from_int_vec.out", "r");

	// check number of nodes in the graph
	ok    = (n == read_number_from_file(fp) );

	// check the list of neighbours for each node
	index = 0;
	for(i = 0; i < n; i++)
	{	// check this node index (both values should be equal i)
		ok &= ( vec[index++] == read_number_from_file(fp) ); 	
		// number of neighbours for this node
		nnbs = read_number_from_file(fp);
		for(j = 0; j < nnbs; j++)
		{	// check this neighbour
			ok &= (vec[index++] == read_number_from_file(fp) );
		}
		// check for -1 to make sure nnbs is correct
		ok &= (vec[index++] == -1);
	}
	// check lenght
	ok &= (index == len);

	// free memory used by the graph
	GMRFLib_free_graph(graph);

	return ok;
}
