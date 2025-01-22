#ifndef __PARAM_CONSTR_H__
#define __PARAM_CONSTR_H__

int inla_parse_param_constraints(inla_tp * mb);
double inla_eval_param_constraint(int thread_id, Data_section_tp * ds);

#endif
