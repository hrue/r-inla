
/* rw.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
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
  \file rw.c
  \brief Handy functions when using RW1, RW2, CRW1, CRW2 and approximate CRW2 models.
*/

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: rw.c,v 1.62 2010/03/10 18:18:08 hrue Exp $ */

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!
  \brief This function returns element Q(i,j), with i=node and j=nnode, of the precision matrix
  for the RW defined in \c rwdef.

  This function allows also arguments ij that are not neighbours and can therefore also be used for
  pruning (using \ref GMRFLib_prune_graph()).

  \param[in] node   First node
  \param[in] nnode  Second node
  \param[in] def  The definition of the RW1 or RW2

  \sa \ref GMRFLib_rwdef_tp, \ref GMRFLib_make_rw_graph, \ref GMRFLib_prune_graph
*/
double GMRFLib_rw(int node, int nnode, void *def)
{
	int imax, imin, idiff, edge;
	double prec;
	GMRFLib_rwdef_tp *rwdef = (GMRFLib_rwdef_tp *) def;

	prec = GMRFLib_SET_PREC(rwdef);
	prec *= (rwdef->prec_scale ? rwdef->prec_scale[0] : 1.0);

	/*
	 * this is the easy case. Note that this case has an additional 'scale0' parameter
	 */
	if (rwdef->order == 0) {
		return (node == nnode ? ((rwdef->scale0 ? rwdef->scale0[node] : 1.0) * prec) : 0.0);
	}

	if (rwdef->cyclic) {
		/*
		 * cyclic, simpler.. 
		 */

		idiff = IMIN(IABS(node - nnode), IABS(IMIN(node, nnode) - IMAX(node, nnode) + rwdef->n));

		switch (rwdef->order) {
		case 1:
			switch (idiff) {
			case 0:
				return 2.0 * prec;
			case 1:
				return -prec;
			default:
				return 0.0;
			}
		case 2:
			switch (idiff) {
			case 0:
				return 6.0 * prec;
			case 1:
				return -4.0 * prec;
			case 2:
				return prec;
			default:
				return 0.0;
			}
		}
	} else {
		imax = IMAX(node, nnode);
		imin = IMIN(node, nnode);
		idiff = imax - imin;

		if (idiff > rwdef->order) {
			return 0.0;			       /* fast return */
		}

		edge = rwdef->order;

		if (imax > edge && imax < rwdef->n - rwdef->order - 1) {
			/*
			 * internal node 
			 */
			switch (rwdef->order) {
			case 1:
				return prec * (idiff == 0 ? 2.0 : -1.0);
			case 2:
				return prec * (idiff == 0 ? 6.0 : (idiff == 1 ? -4.0 : 1.0));
			default:
				return 0.0;
			}
		} else {
			/*
			 * the edge. 
			 */
			if (imax > edge) {
				/*
				 * map to the left egde 
				 */
				imax = rwdef->n - 1 - IMIN(node, nnode);
				imin = rwdef->n - 1 - IMAX(node, nnode);
			}
			switch (rwdef->order) {
			case 1:
				return prec * (idiff == 1 ? -1.0 : (imax == 0 ? 1.0 : 2.0));
			case 2:
				switch (idiff) {
				case 0:
					switch (imax) {
					case 0:
						return prec * 1.0;
					case 1:
						return prec * 5.0;
					case 2:
						return prec * 6.0;
					default:
						GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
					}
				case 1:
					return prec * (imin == 0 ? -2.0 : -4.0);
				case 2:
					return prec * 1.0;
				default:
					GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
				}
			default:
				return 0.0;
			}
		}
	}
	GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);

	return 0.0;
}

/*!
  \brief This function returns element Q(i,j), with i=node and j=nnode, of the precision matrix for
  the CRW defined in \c crwdef.
 
  This function allows also arguments ij that are not neighbours and can therefore also be used for
  pruning (using \ref GMRFLib_prune_graph()).
 
  \param[in]  node  First node
  \param[in] nnode  Second node
  \param[in] def The definition of the CRW1 or CRW2
 
  \sa \ref GMRFLib_crwdef_tp,  \ref GMRFLib_make_crw_graph,  \ref GMRFLib_prune_graph
*/
double GMRFLib_crw(int node, int nnode, void *def)
{
	/*
	 * this is the continous version for order=0, 1 or 2, which take into accont the positions consistently. the order=2
	 * uses the augmentation with the velocity.
	 * 
	 * TODO: make this also accept a cyclic argument. there should not be to much to change here. it is possible that it's
	 * easier if positions are present, then we can use the general expression only `defining' the work-entries correct. in 
	 * any case, then this function could be merged with GMRFLib_rw()? do we then need both? OOPS: if this is change to
	 * accept the cyclic case as well, then rememeber to change the make_graph routine as well! 
	 */
#define TP_POS 0
#define TP_VEL 1
#define SETUP_WORK_PTRS if (1){			\
	idelta  = &crwdef->work[2          ];	\
	idelta2 = &crwdef->work[2 +   (n+4)];	\
	idelta3 = &crwdef->work[2 + 2*(n+4)];	\
        isdelta = &crwdef->work[2 + 3*(n+4)];	\
        sidelta = &crwdef->work[2 + 4*(n+4)];}
#define SETUP_LOCAL_WORK_PTRS if (1){			\
	idelta  = &work[2          ];	\
	idelta2 = &work[2 +   (n+4)];	\
	idelta3 = &work[2 + 2*(n+4)];	\
        isdelta = &work[2 + 3*(n+4)];	\
        sidelta = &work[2 + 4*(n+4)];}

	int i, use_pos, n, idiff, imin, imax, node_i = -1, node_tp = -1, nnode_i = -1, nnode_tp = -1, order;
	double prec;

	double *idelta = NULL, *idelta2 = NULL, *idelta3 = NULL, *isdelta = NULL, *sidelta = NULL;

	GMRFLib_crwdef_tp *crwdef = (GMRFLib_crwdef_tp *) def;

	prec = GMRFLib_SET_PREC(crwdef) * (crwdef->prec_scale ? crwdef->prec_scale[0] : 1.0);
	n = crwdef->n;
	use_pos = (crwdef->position ? 1 : 0);
	order = crwdef->order;

	/*
	 * this is the easy case. Note that this case has an additional 'scale0' parameter
	 */
	if (order == 0) {
		return (node == nnode ? (crwdef->scale0 ? crwdef->scale0[node] : 1.0) * prec : 0.0);
	}
	assert(order > 0);

	if (crwdef->position) {
		/*
		 * check `workspace' and compute it if we don't have it. this is a bit tricky. since this can be called in the OMP
		 * mode, then we do this one thread at the time. any thread can compute the work-space, so we need a second check
		 * to avoid that this workspace is note allocated several times 
		 */
		if (!crwdef->work) {
#pragma omp critical
			{
				/*
				 */
				if (!crwdef->work) {

					/*
					 * we do this simple; we compute the workspace for the union of all cases. 
					 */
					double *delta = NULL, *work = NULL;

					delta = Calloc(n - 1, double);
					for (i = 0; i < n - 1; i++) {
						delta[i] = crwdef->position[i + 1] - crwdef->position[i];
					}
					work = Calloc(5 * (n + 4), double);
					SETUP_LOCAL_WORK_PTRS;
					for (i = 0; i < n - 1; i++) {
						idelta[i] = 1.0 / delta[i];
						idelta2[i] = SQR(idelta[i]);
						idelta3[i] = idelta[i] * idelta2[i];
					}
					for (i = 0; i < n - 2; i++) {
						sidelta[i] = idelta[i] + idelta[i + 1];
						isdelta[i] = 1.0 / (delta[i] + delta[i + 1]);
					}
					sidelta[n - 2] = idelta[n - 2];	/* yes */
					sidelta[-1] = idelta[0];	/* yes */

					Free(delta);
					crwdef->work = work;
				} else {
					SETUP_WORK_PTRS;
				}
			}
		} else {
			SETUP_WORK_PTRS;
		}
	}

	if (order == 1) {
		/*
		 * first order model, this is the simple case with no agumentation 
		 */

		imin = IMIN(node, nnode);
		imax = IMAX(node, nnode);
		idiff = imax - imin;
		if (idiff > order) {
			return 0.0;			       /* fast return, nothing here */
		}

		if (idiff == 0) {
			if (use_pos) {
				if (imin == 0) {
					return prec * idelta[0];
				}
				if (imin == n - 1) {
					return prec * idelta[n - 2];
				}
				return prec * (idelta[imin - 1] + idelta[imin]);
			} else {
				return prec * ((imin == 0 || imin == n - 1) ? 1.0 : 2.0);
			}
		} else {
			return prec * (use_pos ? -idelta[imin] : -1.0);
		}
	}

	/*
	 * the second order model is a bit more involved
	 * 
	 * first there is the switch beteen approximative/exact, then layout, the regular or irregular. 
	 */
	if (crwdef->layout == GMRFLib_CRW_LAYOUT_SIMPLE) {
		/*
		 * this is the approximative scheme with no augmentation.
		 * 
		 * do first the case with regular positions the irregular. the regular case is just a copy from GMRFLib_rw() 
		 */

		imax = IMAX(node, nnode);
		imin = IMIN(node, nnode);
		idiff = imax - imin;
		if (idiff > order) {
			return 0.0;			       /* fast return, nothing here */
		}

		if (!use_pos) {
			/*
			 * regular 
			 */
			if ((imax > order && imax < n - order - 1)) {
				/*
				 * internal node 
				 */
				return prec * (idiff == 0 ? 6.0 : (idiff == 1 ? -4.0 : 1.0));
			} else {
				if (imax > order) {
					/*
					 * map to the left egde 
					 */
					imax = n - 1 - IMIN(node, nnode);
					imin = n - 1 - IMAX(node, nnode);
				}
				switch (idiff) {
				case 0:
					switch (imax) {
					case 0:
						return prec * 1.0;
					case 1:
						return prec * 5.0;
					case 2:
						return prec * 6.0;
					default:
						GMRFLib_ASSERT(0, GMRFLib_ESNH);
					}
				case 1:
					return prec * (imin == 0 ? -2.0 : -4.0);
				case 2:
					return prec * 1.0;
				default:
					GMRFLib_ASSERT(0, GMRFLib_ESNH);
				}
			}
		} else {
			/*
			 * irregular case 
			 */
			switch (idiff) {
			case 0:
				return prec * 2.0 * (idelta2[imax - 1] * isdelta[imax - 2]
						     + idelta[imax - 1] * idelta[imax] * sidelta[imax - 1]
						     + idelta2[imax] * isdelta[imax]);
			case 1:
				return -prec * 2.0 * idelta2[imax - 1] * (idelta[imax - 2] + idelta[imax]);

			case 2:
				return prec * 2.0 * idelta[imax - 2] * idelta[imax - 1] * isdelta[imax - 2];

			default:
				GMRFLib_ASSERT(0, GMRFLib_ESNH);
			}
		}
	} else {
		/*
		 * this is the exact solution, which require augmentation 
		 */

		switch (crwdef->layout) {
			/*
			 * we need to make sure that node_i < nnode_i 
			 */

		case GMRFLib_CRW_LAYOUT_PAIRS:
			if (node < nnode) {
				node_i = node / 2;
				nnode_i = nnode / 2;
				node_tp = (GSL_IS_EVEN(node) ? TP_POS : TP_VEL);
				nnode_tp = (GSL_IS_EVEN(nnode) ? TP_POS : TP_VEL);
			} else {
				node_i = nnode / 2;
				nnode_i = node / 2;
				node_tp = (GSL_IS_EVEN(nnode) ? TP_POS : TP_VEL);
				nnode_tp = (GSL_IS_EVEN(node) ? TP_POS : TP_VEL);
			}
			break;
		case GMRFLib_CRW_LAYOUT_BLOCK:
			if (node % n < nnode % n) {
				node_i = node % n;
				nnode_i = nnode % n;
				node_tp = (node < n ? TP_POS : TP_VEL);
				nnode_tp = (nnode < n ? TP_POS : TP_VEL);
			} else {
				node_i = nnode % n;
				nnode_i = node % n;
				node_tp = (nnode < n ? TP_POS : TP_VEL);
				nnode_tp = (node < n ? TP_POS : TP_VEL);
			}
			break;
		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
		}

		idiff = nnode_i - node_i;
		assert(idiff >= 0);

		if (idiff > 1) {
			return 0.0;			       /* nothing here */
		}

		/*
		 * split in two cases to speedup 
		 */

		if (use_pos) {
			if (idiff == 0) {
				if (node_tp == TP_POS) {
					if (nnode_tp == TP_POS) {	/* TP_POS & TP_POS */
						if (nnode_i == 0) {
							return prec * 12.0 * idelta3[0];
						}
						if (nnode_i == n - 1) {
							return prec * 12.0 * idelta3[n - 2];
						}
						return prec * 12.0 * (idelta3[node_i] + idelta3[node_i - 1]);
					} else {	       /* TP_POS & TP_VEL */
						if (nnode_i == 0) {
							return prec * 6.0 * idelta2[0];
						}
						if (nnode_i == n - 1) {
							return prec * (-6.0) * idelta2[n - 2];
						}
						return prec * 6.0 * (idelta2[node_i] - idelta2[node_i - 1]);
					}
				} else {
					if (nnode_tp == TP_POS) {	/* TP_VEL & TP_POS */
						if (nnode_i == 0) {
							return prec * 6.0 * idelta2[0];
						}
						if (nnode_i == n - 1) {
							return prec * (-6.0 * idelta2[n - 2]);
						}
						return prec * 6.0 * (idelta2[node_i] - idelta2[node_i - 1]);
					} else {	       /* TP_VEL & TP_VEL */
						if (nnode_i == 0) {
							return prec * 4.0 * idelta[0];
						}
						if (nnode_i == n - 1) {
							return prec * 4.0 * idelta[n - 2];
						}
						return prec * 4.0 * (idelta[node_i] + idelta[node_i - 1]);
					}
				}
			} else {
				if (node_tp == TP_POS) {
					if (nnode_tp == TP_POS) {
						return prec * (-12.0) * idelta3[node_i];
					} else {
						return prec * 6.0 * idelta2[node_i];
					}
				} else {
					if (nnode_tp == TP_POS) {
						return prec * (-6.0) * idelta2[node_i];
					} else {
						return prec * 2.0 * idelta[node_i];
					}
				}
			}
		} else {
			if (idiff == 0) {
				if (node_tp == TP_POS) {
					if (nnode_tp == TP_POS) {	/* TP_POS & TP_POS */
						if (nnode_i == 0) {
							return prec * 12.0;
						}
						if (nnode_i == n - 1) {
							return prec * 12.0;
						}
						return prec * 24.0;
					} else {	       /* TP_POS & TP_VEL */
						if (nnode_i == 0) {
							return prec * 6.0;
						}
						if (nnode_i == n - 1) {
							return prec * (-6.0);
						}
						return 0.0;
					}
				} else {
					if (nnode_tp == TP_POS) {	/* TP_VEL & TP_POS */
						if (nnode_i == 0) {
							return prec * 6.0;
						}
						if (nnode_i == n - 1) {
							return prec * (-6.0);
						}
						return 0.0;
					} else {	       /* TP_VEL & TP_VEL */
						if (nnode_i == 0) {
							return prec * 4.0;
						}
						if (nnode_i == n - 1) {
							return prec * 4.0;
						}
						return prec * 8.0;
					}
				}
			} else {
				if (node_tp == TP_POS) {
					if (nnode_tp == TP_POS) {
						return prec * (-12.0);	/* POS & POS */
					} else {
						return prec * 6.0;	/* POS & VEL */
					}
				} else {
					if (nnode_tp == TP_POS) {
						return prec * (-6.0);	/* VEL & POS */
					} else {
						return prec * 2.0;	/* VEL & VEL */
					}
				}
			}
		}
		return 0.0;
	}

	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);

#undef SETUP_WORK_PTRS
#undef SETUP_LOCAL_WORK_PTRS
#undef TP_POS
#undef TP_VEL
}

/*!
  \brief This function returns element Q(i,j), with i=node and j=nnode, of the precision matrix for
  the RW on a 2D lattice defined in \c rw2ddef.

  Currently only order=2 is supported, however, it does incorporate correct boundary conditions for
  the non-cyclic case.
 
  \param[in]  node  First node
  \param[in] nnode  Second node
  \param[in] def The definition of the RW on the 2D lattice
 
  \sa \ref GMRFLib_rw2ddef_tp,  \ref GMRFLib_make_rw2d_graph
*/
double GMRFLib_rw2d(int node, int nnode, void *def)
{
	/*
	 * this function does work for non-neigbour arguments, hence it might be (and is) used for pruning. 
	 */
	int i, j, ii, jj, dx, dy, nrow, ncol, iref, jref, imin, jmin, imax, jmax;
	double prec;
	int map[] = { 0, 1, 2, 1, 0 };

	GMRFLib_rw2ddef_tp *rw2ddef = (GMRFLib_rw2ddef_tp *) def;

	nrow = rw2ddef->nrow;
	ncol = rw2ddef->ncol;
	prec = GMRFLib_SET_PREC(rw2ddef);
	prec *= (rw2ddef->prec_scale ? rw2ddef->prec_scale[0] : 1.0);

	GMRFLib_node2lattice(node, &i, &j, nrow, ncol);
	GMRFLib_node2lattice(nnode, &ii, &jj, nrow, ncol);

	if (rw2ddef->cyclic) {
		/*
		 * the cyclic case 
		 */

		dx = IMIN(IABS(i - ii), IABS(IMIN(i, ii) - IMAX(i, ii) + nrow));
		dy = IMIN(IABS(j - jj), IABS(IMIN(j, jj) - IMAX(j, jj) + ncol));
		dx = IABS(dx);
		dy = IABS(dy);

		switch (dx) {
		case 0:
			switch (dy) {
			case 0:
				return 20.0 * prec;
			case 1:
				return -8.0 * prec;
			case 2:
				return prec;
			default:
				return 0.0;
			}
		case 1:
			switch (dy) {
			case 0:
				return -8.0 * prec;
			case 1:
				return 2.0 * prec;
			default:
				return 0.0;
			}
		case 2:
			switch (dy) {
			case 0:
				return prec;
			default:
				return 0.0;
			}
		default:
			return 0.0;
		}
	} else if (rw2ddef->bvalue == GMRFLib_BVALUE_ZERO) {
		/*
		 * this is like the RW2D on a larger lattice, but here we condition on 0 outside 
		 */
		dx = IABS(i - ii);
		dy = IABS(j - jj);

		switch (dx) {
		case 0:
			switch (dy) {
			case 0:
				return 20.0 * prec;
			case 1:
				return -8.0 * prec;
			case 2:
				return prec;
			default:
				return 0.0;
			}
		case 1:
			switch (dy) {
			case 0:
				return -8.0 * prec;
			case 1:
				return 2.0 * prec;
			default:
				return 0.0;
			}
		case 2:
			switch (dy) {
			case 0:
				return prec;
			default:
				return 0.0;
			}
		default:
			return 0.0;
		}
	} else {
		/*
		 * this is with boundary effect, this is messy
		 * 
		 * first, rule out non-neighbours in the graph, so that the remaining terms are non-zero. this is similar to
		 * the prune function for the GMRFLib_make_rw2d_graph() 
		 */

		dx = IABS(i - ii);
		dy = IABS(j - jj);

		if ((IMAX(dx, dy) > 2) || (IMAX(dx, dy) == 2 && IMIN(dx, dy) >= 1)) {
			return 0.0;
		}

		/*
		 * these conditions are always true, no matter boundary or not 
		 */
		if (IMAX(dx, dy) == 2) {
			return prec;
		}
		if (dx == 1 && dy == 1) {
			return 2.0 * prec;
		}

		/*
		 * this is the internal case 
		 */
		if ((i > 1 && i < nrow - 2 && j > 1 && j < ncol - 2) || (ii > 1 && ii < nrow - 2 && jj > 1 && jj < ncol - 2)) {
			if (dx == 0) {
				return (dy == 0 ? 20.0 : -8.0) * prec;
			} else {
				return -8.0 * prec;
			}
		}

		/*
		 * boundaries! are both equal? 
		 */
		if ((i == ii) && (j == jj)) {
			int itmp, jtmp;

			itmp = (i > 1 ? map[IMAX(2, i - nrow + 5)] : i);
			jtmp = (j > 1 ? map[IMAX(2, j - ncol + 5)] : j);

			iref = IMAX(itmp, jtmp);
			jref = IMIN(itmp, jtmp);

			if (iref == 1) {
				return (jref == 1 ? 18.0 : 10.0) * prec;
			}
			if (iref == 2) {
				return (jref == 1 ? 19.0 : 11.0) * prec;
			}
			return 4.0 * prec;
		}

		/*
		 * ...then do the horizontal and vertical neighbours 
		 */
		imax = IMAX(i, ii);
		imin = IMIN(i, ii);
		jmax = IMAX(j, jj);
		jmin = IMIN(j, jj);

		/*
		 * corners? 
		 */
		if ((imin == 0 || imax == nrow - 1) && (jmin == 0 || jmax == ncol - 1)) {
			return -4.0 * prec;
		}

		/*
		 * at the edge? 
		 */
		if (imin == 0 || imax == nrow - 1 || jmin == 0 || jmax == ncol - 1) {
			return -6.0 * prec;
		}

		/*
		 * the remaining case 
		 */
		return -8.0 * prec;
	}

	return 0.0;
}

/*!
  \brief Make the graph suitable to the RW1 or RW2 model with regular locations defined in def.

  \param[out] graph  The graph for the RW1 or RW2 model with regular locations

  \param[in] def The definition of the RW1 or RW2 model with regular locations

  \remark \c GMRFLib_make_rw_graph() is to be used in connection with \c GMRFLib_rw() both using the
  RW1 or RW2 model defined using \c GMRFLib_rwdef_tp.

  \remark There is an alternative set of tools for the case where the locations are irregular, see
  \c GMRFLib_make_crw_graph(), \c GMRFLib_crw() and \c GMRFLib_crwdef_tp.

  \sa GMRFLib_rw, GMRFLib_rwdef_tp
*/
int GMRFLib_make_rw_graph(GMRFLib_graph_tp ** graph, GMRFLib_rwdef_tp * def)
{
	GMRFLib_make_linear_graph(graph, def->n, def->order, def->cyclic);
	return GMRFLib_SUCCESS;
}


/*!
  \brief Make the graph suitable to the (C)RW1 or CRW2 model with irregular locations defined in def.
  
  \param[out] graph  The graph for the (C)RW1 or CRW2 model with irregular locations
  
  \param[in] def The definition of the (C)RW1 or CRW2 model with irregular locations
  
  \remark \c GMRFLib_make_crw_graph() is to be used in connection with \c GMRFLib_crw() both using
  the (C)RW1 or CRW2 model defined using \c GMRFLib_crwdef_tp.
    
  \remark There is an alternative set of tools for the case where the locations are regular, see \c
  GMRFLib_make_rw_graph(), \c GMRFLib_rw() and \c GMRFLib_rwdef_tp.

  \sa GMRFLib_crw(), GMRFLib_crwdef_tp
*/
int GMRFLib_make_crw_graph(GMRFLib_graph_tp ** graph, GMRFLib_crwdef_tp * def)
{
	int i, *hold = NULL, n;
	GMRFLib_graph_tp *gg = NULL;

	n = def->n;

	if (def->order <= 1 || (def->order == 2 && def->layout == GMRFLib_CRW_LAYOUT_SIMPLE)) {
		// FIXME("\n\n!!!!modify the graph to a complete graph!!!");GMRFLib_make_linear_graph(graph, n, n, 0);
		GMRFLib_make_linear_graph(graph, n, def->order, 0);
		return GMRFLib_SUCCESS;
	}
	if (def->order == 2 && def->layout == GMRFLib_CRW_LAYOUT_PAIRS) {
		GMRFLib_make_linear_graph(graph, 2 * n, 3, 0);
		return GMRFLib_SUCCESS;
	}

	/*
	 * the rest is for the block'ed case for order 2. this is the tricky bit 
	 */
	GMRFLib_make_empty_graph(&gg);
	gg->n = 2 * n;
	gg->nnbs = Calloc(2 * n, int);
	gg->nbs = Calloc(2 * n, int *);
	hold = Calloc(2 * n * 5, int);

	for (i = 0; i < 2 * n; i++) {
		gg->nbs[i] = &hold[i * 5];		       /* this is ok */
	}

	for (i = 0; i < n; i++) {
		if (i == 0) {
			gg->nnbs[i] = 3;
			gg->nbs[i][0] = i + n;		       /* vel */
			gg->nbs[i][1] = i + 1;		       /* next pos */
			gg->nbs[i][2] = i + 1 + n;	       /* next vel */

			gg->nnbs[i + n] = 3;
			gg->nbs[i + n][0] = i;		       /* pos */
			gg->nbs[i + n][1] = i + 1;	       /* next pos */
			gg->nbs[i + n][2] = i + 1 + n;	       /* next vel */
		} else if (i == n - 1) {
			gg->nnbs[i] = 3;
			gg->nbs[i][0] = i + n;		       /* vel */
			gg->nbs[i][1] = i - 1;		       /* prev pos */
			gg->nbs[i][2] = i - 1 + n;	       /* prev vel */

			gg->nnbs[i + n] = 3;
			gg->nbs[i + n][0] = i;		       /* pos */
			gg->nbs[i + n][1] = i - 1;	       /* prev pos */
			gg->nbs[i + n][2] = i - 1 + n;	       /* prev vel */
		} else {
			gg->nnbs[i] = 5;
			gg->nbs[i][0] = i + n;		       /* vel */
			gg->nbs[i][1] = i + 1;		       /* next pos */
			gg->nbs[i][2] = i + 1 + n;	       /* next vel */
			gg->nbs[i][3] = i - 1;		       /* prev pos */
			gg->nbs[i][4] = i - 1 + n;	       /* prev vel */

			gg->nnbs[i + n] = 5;
			gg->nbs[i + n][0] = i;		       /* pos */
			gg->nbs[i + n][1] = i + 1;	       /* next pos */
			gg->nbs[i + n][2] = i + 1 + n;	       /* next vel */
			gg->nbs[i + n][3] = i - 1;	       /* prev pos */
			gg->nbs[i + n][4] = i - 1 + n;	       /* prev vel */
		}
	}
	GMRFLib_EWRAP0(GMRFLib_prepare_graph(gg));
	*graph = gg;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Makes the graph suitable to the RW model on a 2D lattice in \c def.
  
  This function creates a lattice graph for the RW model on a 2D lattice as defined in \c
  def. Currently only order=2 is supported. The dimention of the lattice must be at least, 5 x
  5. The mapping of the nodes to pixel indices, is defined similar to the \ref
  GMRFLib_make_lattice_graph(), ie using the functions \ref GMRFLib_node2lattice() and \ref
  GMRFLib_lattice2node().

  \param[out] graph  The graph for the RW model on a 2D lattice
  
  \param[in] def The definition of the RW model on a 2D lattice
  
  \sa GMRFLib_rw2d(), GMRFLib_rw2ddef_tp, GMRFLib_node2lattice(), GMRFLib_lattice2node()
*/
int GMRFLib_make_rw2d_graph(GMRFLib_graph_tp ** graph, GMRFLib_rw2ddef_tp * def)
{
	GMRFLib_graph_tp *g = NULL;

	GMRFLib_EWRAP0(GMRFLib_make_lattice_graph(&g, def->nrow, def->ncol, 2, 2, def->cyclic));
	GMRFLib_EWRAP0(GMRFLib_prune_graph(graph, g, (GMRFLib_Qfunc_tp *) GMRFLib_rw2d, (void *) def));

	GMRFLib_free_graph(g);

	return GMRFLib_SUCCESS;
}


int GMRFLib_crw_scale(void *def)
{
	/*
	 * This approach uses the constrained sampling approach, much faster
	 */
	GMRFLib_crwdef_tp *crwdef = Calloc(1, GMRFLib_crwdef_tp);
	GMRFLib_crwdef_tp *odef = (GMRFLib_crwdef_tp *) def;

	double *prec_scale_guess = Calloc(1, double);
	*prec_scale_guess = 1.0;

	crwdef->n = odef->n;
	assert(odef->order > 0);
	crwdef->order = odef->order;
	crwdef->prec = NULL;
	crwdef->log_prec = NULL;
	crwdef->log_prec_omp = NULL;
	crwdef->position = odef->position;
	assert(odef->layout == GMRFLib_CRW_LAYOUT_SIMPLE);
	crwdef->layout = odef->layout;
	crwdef->work = NULL;
	crwdef->scale0 = NULL;
	crwdef->prec_scale = prec_scale_guess;

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_make_crw_graph(&graph, crwdef);

	int i, free_position = 0;

	/*
	 * make sure we have defined the positions, as the code is easier with it
	 */
	if (!(crwdef->position)) {
		free_position = 1;
		crwdef->position = Calloc(graph->n, double);
		for (i = 0; i < graph->n; i++) {
			crwdef->position[i] = i;
		}
	}

	double *len = Calloc(graph->n, double);
	for (i = 0; i < graph->n; i++) {
		if (i == 0) {
			len[i] = (crwdef->position[i + 1] - crwdef->position[i]) / 2.0;
		} else if (i == graph->n - 1) {
			len[i] = (crwdef->position[i] - crwdef->position[i - 1]) / 2.0;
		} else {
			len[i] = (crwdef->position[i + 1] - crwdef->position[i - 1]) / 2.0;
		}
	}

	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);
	constr->nc = crwdef->order;
	constr->a_matrix = Calloc(constr->nc * graph->n, double);
	for (i = 0; i < graph->n; i++) {
		constr->a_matrix[i * constr->nc + 0] = len[i];
	}

	double len_acum = 0.0;
	for (i = 0; i < graph->n; i++) {
		len_acum += len[i];
	}
	if (crwdef->order == 1) {
		*prec_scale_guess = len_acum;
	} else {
		*prec_scale_guess = SQR(len_acum);
	}

	if (crwdef->order == 2) {
		len_acum = 0.0;
		for (i = 0; i < graph->n; i++) {
			len_acum += len[i];
			constr->a_matrix[i * constr->nc + 1] = len_acum;
		}
	}

	constr->e_vector = Calloc(constr->nc, double);
	GMRFLib_prepare_constr(constr, graph, GMRFLib_TRUE);

	double *c = Calloc(graph->n, double), eps = GMRFLib_eps(0.5);
	GMRFLib_problem_tp *problem;

	for (i = 0; i < graph->n; i++) {
		c[i] = eps;
	}

	int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
	GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

	while (!ok) {
		retval =
		    GMRFLib_init_problem(&problem, NULL, NULL, c, NULL, graph, GMRFLib_crw, (void *) crwdef, NULL, constr, GMRFLib_NEW_PROBLEM);
		switch (retval) {
		case GMRFLib_EPOSDEF:
			for (i = 0; i < graph->n; i++) {
				c[i] *= 10.0;
			}
			problem = NULL;
			break;
		case GMRFLib_SUCCESS:
			ok = 1;
			break;
		default:
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_ERROR(retval);
			abort();
			break;
		}

		if (++num_try >= num_try_max) {
			FIXME("This should not happen. Contact developers...");
			abort();
		}
	}
	GMRFLib_set_error_handler(old_handler);

	GMRFLib_Qinv(problem, GMRFLib_QINV_DIAG);

	double sum = 0.0;
	for (i = 0; i < graph->n; i++) {
		sum += log(*(GMRFLib_Qinv_get(problem, i, i))) * len[i];
	}

	odef->prec_scale = Calloc(1, double);
	odef->prec_scale[0] = exp(sum / (crwdef->position[graph->n - 1] - crwdef->position[0])) * *prec_scale_guess;

	Free(c);
	Free(crwdef);
	Free(len);
	GMRFLib_free_constr(constr);
	GMRFLib_free_graph(graph);
	GMRFLib_free_problem(problem);
	if (free_position) {
		Free(crwdef->position);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_rw_scale(void *def)
{
	GMRFLib_rwdef_tp *rwdef = Calloc(1, GMRFLib_rwdef_tp);
	GMRFLib_rwdef_tp *odef = (GMRFLib_rwdef_tp *) def;
	double *prec_scale_guess = Calloc(1, double);
	*prec_scale_guess = 1.0;

	rwdef->n = odef->n;
	assert(odef->order > 0);
	rwdef->order = odef->order;
	rwdef->cyclic = odef->cyclic;
	rwdef->prec = NULL;
	rwdef->log_prec = NULL;
	rwdef->log_prec_omp = NULL;
	rwdef->scale0 = odef->scale0;
	rwdef->prec_scale = prec_scale_guess;

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_make_rw_graph(&graph, rwdef);
	int i;
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_make_empty_constr(&constr);

	if (!rwdef->cyclic) {
		/*
		 * cyclic == FALSE
		 */
		if (rwdef->order == 0) {
			constr->nc = 0;
		} else if (rwdef->order == 1) {
			constr->nc = 1;
			constr->a_matrix = Calloc(constr->nc * graph->n, double);
			for (i = 0; i < graph->n; i++) {
				constr->a_matrix[i * constr->nc + 0] = 1.0;
			}
			*prec_scale_guess = (double) graph->n - 1;
		} else if (rwdef->order == 2) {
			constr->nc = 2;
			constr->a_matrix = Calloc(constr->nc * graph->n, double);
			for (i = 0; i < graph->n; i++) {
				constr->a_matrix[i * constr->nc + 0] = 1.0;
				constr->a_matrix[i * constr->nc + 1] = (i - graph->n / 2.0);
			}
			*prec_scale_guess = (double) ISQR(graph->n - 1);
		} else {
			assert(0 == 1);
		}
	} else {
		/*
		 * cyclic == TRUE
		 */
		if (rwdef->order == 0) {
			constr->nc = 0;
		} else if (rwdef->order == 1 || rwdef->order == 2) {
			constr->nc = 1;
			constr->a_matrix = Calloc(constr->nc * graph->n, double);
			for (i = 0; i < graph->n; i++) {
				constr->a_matrix[i * constr->nc + 0] = 1.0;
			}
			*prec_scale_guess = (double) (rwdef->order == 1 ? graph->n - 1 : ISQR(graph->n - 1));
		} else {
			assert(0 == 1);
		}
	}

	if (constr->nc) {
		constr->e_vector = Calloc(constr->nc, double);
		GMRFLib_prepare_constr(constr, graph, GMRFLib_TRUE);
	} else {
		GMRFLib_free_constr(constr);
		constr = NULL;
	}

	double *c = Calloc(graph->n, double), eps = GMRFLib_eps(0.5);
	GMRFLib_problem_tp *problem;

	for (i = 0; i < graph->n; i++) {
		c[i] = eps;
	}

	int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
	GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

	while (!ok) {
		retval = GMRFLib_init_problem(&problem, NULL, NULL, c, NULL, graph, GMRFLib_rw, (void *) rwdef, NULL, constr, GMRFLib_NEW_PROBLEM);
		switch (retval) {
		case GMRFLib_EPOSDEF:
			for (i = 0; i < graph->n; i++) {
				c[i] *= 10.0;
			}
			problem = NULL;
			break;
		case GMRFLib_SUCCESS:
			ok = 1;
			break;
		default:
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_ERROR(retval);
			abort();
			break;
		}

		if (++num_try >= num_try_max) {
			FIXME("This should not happen. Contact developers...");
			abort();
		}
	}
	GMRFLib_set_error_handler(old_handler);

	GMRFLib_Qinv(problem, GMRFLib_QINV_DIAG);

	double sum = 0.0;
	for (i = 0; i < graph->n; i++) {
		sum += log(*(GMRFLib_Qinv_get(problem, i, i)));
	}

	odef->prec_scale = Calloc(1, double);
	odef->prec_scale[0] = exp(sum / graph->n) * *prec_scale_guess;

	Free(c);
	Free(rwdef);
	GMRFLib_free_constr(constr);
	GMRFLib_free_graph(graph);
	GMRFLib_free_problem(problem);

	return GMRFLib_SUCCESS;
}

int GMRFLib_rw2d_scale(void *def)
{
	GMRFLib_rw2ddef_tp *rw2ddef = Calloc(1, GMRFLib_rw2ddef_tp);
	GMRFLib_rw2ddef_tp *odef = (GMRFLib_rw2ddef_tp *) def;

	rw2ddef->nrow = odef->nrow;
	rw2ddef->ncol = odef->ncol;
	rw2ddef->order = odef->order;
	assert(rw2ddef->order == 2);
	rw2ddef->cyclic = odef->cyclic;
	rw2ddef->bvalue = odef->bvalue;
	rw2ddef->log_prec = NULL;
	rw2ddef->log_prec_omp = NULL;
	rw2ddef->prec_scale = NULL;

	GMRFLib_graph_tp *graph = NULL;
	GMRFLib_make_rw2d_graph(&graph, rw2ddef);

	int i, j, k;
	double *c = NULL;
	GMRFLib_constr_tp *constr = NULL;
	GMRFLib_problem_tp *problem = NULL;

	if (rw2ddef->bvalue != GMRFLib_BVALUE_ZERO) {
		GMRFLib_make_empty_constr(&constr);
		constr->nc = (rw2ddef->cyclic ? 1 : 3);
		constr->a_matrix = Calloc(constr->nc * graph->n, double);
		if (constr->nc == 1) {
			for (i = 0; i < graph->n; i++) {
				constr->a_matrix[i * constr->nc + 0] = 1.0;
			}
		} else {
			for (j = 0; j < rw2ddef->ncol; j++) {
				for (i = 0; i < rw2ddef->nrow; i++) {
					GMRFLib_lattice2node(&k, i, j, rw2ddef->nrow, rw2ddef->ncol);
					constr->a_matrix[k * constr->nc + 0] = 1;
					constr->a_matrix[k * constr->nc + 1] = i - rw2ddef->nrow / 2.0;
					constr->a_matrix[k * constr->nc + 2] = j - rw2ddef->ncol / 2.0;
				}
			}
		}

		constr->e_vector = Calloc(constr->nc, double);
		GMRFLib_prepare_constr(constr, graph, GMRFLib_TRUE);

		double eps = GMRFLib_eps(0.75);
		c = Calloc(graph->n, double);
		for (i = 0; i < graph->n; i++) {
			c[i] = eps;
		}
	} else {
		/*
		 * The model is proper if BVALUE_ZERO is set, no need to add anything on the diagonal
		 */
		constr = NULL;
	}

	int retval = GMRFLib_SUCCESS, ok = 0, num_try = 0, num_try_max = 100;
	GMRFLib_error_handler_tp *old_handler = GMRFLib_set_error_handler_off();

	while (!ok) {
		retval =
		    GMRFLib_init_problem(&problem, NULL, NULL, c, NULL, graph, GMRFLib_rw2d, (void *) rw2ddef, NULL, constr, GMRFLib_NEW_PROBLEM);
		switch (retval) {
		case GMRFLib_EPOSDEF:
			for (i = 0; i < graph->n; i++) {
				c[i] *= 10.0;
			}
			problem = NULL;
			break;
		case GMRFLib_SUCCESS:
			ok = 1;
			break;
		default:
			GMRFLib_set_error_handler(old_handler);
			GMRFLib_ERROR(retval);
			abort();
			break;
		}

		if (++num_try >= num_try_max) {
			FIXME("This should not happen. Contact developers...");
			abort();
		}
	}
	GMRFLib_set_error_handler(old_handler);

	GMRFLib_Qinv(problem, GMRFLib_QINV_DIAG);

	double sum = 0.0;
	for (i = 0; i < graph->n; i++) {
		sum += log(*(GMRFLib_Qinv_get(problem, i, i)));
	}

	odef->prec_scale = Calloc(1, double);
	odef->prec_scale[0] = exp(sum / graph->n);

	Free(c);
	Free(rw2ddef);
	GMRFLib_free_constr(constr);
	GMRFLib_free_graph(graph);
	GMRFLib_free_problem(problem);

	return GMRFLib_SUCCESS;
}

/*
  Example for manual
 */

/*! \page ex_rw A worked out example smoothing a time-series data, using the routines in rw.c
  
Solve the same problem as in \ref ex_wa, now using the routines in rw.c

\par Program code:

\verbinclude example-doxygen-rw.txt

*/
