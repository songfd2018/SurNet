#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
# include <chrono> // time
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h>  //mkdir
#include <algorithm>    // sort
#include <vector> 


using namespace std;

extern "C" {

	int rand_cate(double* _prop, omprng _rng) {
		int res = 0;
		double u = _rng.runif();
		while (u > _prop[res]) {
			u = u - _prop[res];
			res++;
		}
		return res;
	}

	int rand_Ber(double _prob, omprng _rng) {
		int res = 1;
		double u = _rng.runif();
		if (u > _prob) {
			res = 0;
		}
		return res;
	}

	int rand_Poi(double _lam, omprng _rng) {
		int res = 1;
		double rat, p;
		rat = exp(-_lam);
		p = _rng.runif();
		while (p > rat) {
			p = p * _rng.runif();
			res++;
		}
		res--;
		return res;
	}

	int rand_discrete(int _dim, double* _prob, omprng _rng) {
		int res = 0;
		double u = _rng.runif();
		while (u > _prob[res]) {
			u = u - _prob[res];
			res++;
		}
		return res;
	}

	void rand_Dir(double *_xi, int _K, omprng _rng, double *_pi) {

		double *rn_gam = new double[_K];
		double sum_gam = 0.0;
		for (int k = 0; k < _K; k++) {
			rn_gam[k] = _rng.rgamma(_xi[k], 1.0);
			sum_gam += rn_gam[k];
		}
		for (int k = 0; k < _K; k++) {
			_pi[k] = rn_gam[k] / sum_gam;
		}

		delete[] rn_gam;
	}

	//Max value in a vector
	double vec_max(double* vec, int n) {
		double res = vec[0];
		for (int i = 1; i < n; i++) {
			if (res < vec[i]) {
				res = vec[i];
			}
		}
		return res;
	}

	void _update_N(int _P, int _n_eml, double* _knots,
		double** _X, double* _Y, int* _nu, int* _int_eml,
		double** _lambda, double** _beta,
		int _mem_pair, omprng _rng,
		int* _intention) {

		double exp_xbeta, fail_prob, lam_poi;
		int int_time;
		for (int k = 0; k < _n_eml; k++) {
			int_time = _int_eml[k] - 1; // the corresponding time interval for Y_{ijk}
			// cout << "int_time = " << int_time << endl;
			// exp(x_{ijk}^T beta_{lm}) 
			exp_xbeta = 1.0;
			for (int p = 0; p < _P; p++) {
				exp_xbeta = exp_xbeta * exp(_X[k][p] * _beta[p][_mem_pair]);
			}

			// cout << "exp_xbeta = " << exp_xbeta << endl;
			// ( 1- S_pch(y_ij | lambda_lm) )
			fail_prob = 1.0;
			for (int j = 0; j < int_time; j++) {
				fail_prob = fail_prob * exp(-_lambda[j][_mem_pair] * (_knots[j + 1] - _knots[j]));
			}
			fail_prob = fail_prob * exp(-_lambda[int_time][_mem_pair] * (_Y[k] - _knots[int_time]));
			fail_prob = 1.0 - fail_prob;
			// cout << "fail_prob = " << fail_prob << endl;

			lam_poi = exp_xbeta * (1 - fail_prob);
			
			//if (k == 0) {
			//	cout << "The " << k << " reply time of lambda  = " << lam_poi << endl;
			//	cout << "Membership pair indicator is " << _mem_pair << endl;
			//	cout << "exp_xbeta = " << exp_xbeta << ", and fail_prob = " << fail_prob << endl;
			//}
			_intention[k] = rand_Poi(lam_poi, _rng) + _nu[k];
		}
	}

	int _update_mp(int _C, int _J, int _P, int _n_eml, double* _knots,
		int* _fail_ind, double* _sur_interval, double* _sur_cov,
		double** _X, double* _Y, int* _nu, int* _int_eml,
		double** _lambda, double** _beta, double* _pi_sender, double* _pi_receiver,
		omprng _rng) {

		int res;
		double* log_p = new double[_C * _C];
		double* event_prob = new double[_C * _C];
		double exp_xbeta, fail_prob, max_log, sum_prob;
		int int_time, ind_mem;

		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				log_p[ind_mem] = log(_pi_sender[l]) + log(_pi_receiver[m]);

				for (int j = 0; j < _J; j++) {
					// Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
					log_p[ind_mem] = log_p[ind_mem] + _fail_ind[j] * log(_lambda[j][ind_mem]);
					// Corresponding to nu_ij (- lambda_p(y_ij - s_{p-1}) - sum_q = 1 ^ p-1 lambda_q (s_q-s_{q-1}))
					//                             = lambda_{lm}^T tilde A_{ij}
					log_p[ind_mem] = log_p[ind_mem] - _sur_interval[j] * _lambda[j][ind_mem];
				}

				for (int p = 0; p < _P; p++) {
					// Corresponding to nu_ij x_ijk^T beta_lm
					log_p[ind_mem] = log_p[ind_mem] + _sur_cov[p] * _beta[p][ind_mem];
				}

				// Corresponding to - exp(x_{ijk}^T beta_{lm}) * ( 1- S_pch(y_ij | lambda_lm) )
				for (int k = 0; k < _n_eml; k++) {
					int_time = _int_eml[k] - 1; // the corresponding time interval for Y_{ijk}
					// cout << "int_time = " << int_time << endl;
					// exp(x_{ijk}^T beta_{lm}) 
					exp_xbeta = 1.0;
					for (int p = 0; p < _P; p++) {
						exp_xbeta = exp_xbeta * exp(_X[k][p] * _beta[p][ind_mem]);
					}

					// cout << "exp_xbeta = " << exp_xbeta << endl;
				   // ( 1- S_pch(y_ij | lambda_lm) )
					fail_prob = 1.0;
					for (int j = 0; j < int_time; j++) {
						fail_prob = fail_prob * exp(-_lambda[j][ind_mem] * (_knots[j + 1] - _knots[j]));
					}
					fail_prob = fail_prob * exp(-_lambda[int_time][ind_mem] * (_Y[k] - _knots[int_time]));
					fail_prob = 1.0 - fail_prob;
					// cout << "fail_prob = " << fail_prob << endl;

					log_p[ind_mem] = log_p[ind_mem] - exp_xbeta * fail_prob;
					// cout << "log_prob["<< i << "]["<< ind_mem << "] = " << log_prob[i][ind_mem] << endl;

				}
			}
		}

		max_log = vec_max(log_p, _C * _C);
		sum_prob = 0.0;
		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				event_prob[ind_mem] = exp(log_p[ind_mem] - max_log);
				sum_prob = sum_prob + event_prob[ind_mem];
			}
		}

		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				event_prob[ind_mem] = event_prob[ind_mem] / sum_prob;
				// cout << "event_prob[" << l << "][" << m << "] = " << event_prob[ind_mem] << ", ";
			}
		}

		res = rand_discrete(_C * _C, event_prob, _rng);

		delete[] event_prob;
		delete[] log_p;

		return res;
	}

	int _update_mp_MH(int _C, int _J, int _P, int _n_eml, double* _knots,
		int* _fail_ind, double* _sur_interval, double* _sur_cov,
		double** _X, double* _Y, int* _nu, int* _int_eml,
		double** _lambda, double** _beta, double* _pi_sender, double* _pi_receiver,
		omprng _rng, int mp_cur) {

		double* event_prob = new double[_C * _C];
		double log_rho, exp_xbeta_cur, exp_xbeta_temp, fail_prob_cur, fail_prob_temp;
		int int_time, ind_mem, mp_temp;
		int l_cur, m_cur, l_temp, m_temp;
		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				event_prob[ind_mem] = 1.0 / _C / _C;
			}
		}
		mp_temp = rand_discrete(_C * _C, event_prob, _rng);

		l_cur = mp_cur / _C;
		m_cur = mp_cur - _C * l_cur;
		l_temp = mp_temp / _C;
		m_temp = mp_temp - _C * l_temp;

		log_rho = 0.0;

		log_rho = log(_pi_sender[l_temp]) + log(_pi_receiver[m_temp]) - log(_pi_sender[l_cur]) - log(_pi_receiver[m_cur]);

		for (int j = 0; j < _J; j++) {
			// Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
			log_rho = log_rho + _fail_ind[j] * (log(_lambda[j][mp_temp]) - log(_lambda[j][mp_cur]));
			// Corresponding to nu_ij (- lambda_p(y_ij - s_{p-1}) - sum_q = 1 ^ p-1 lambda_q (s_q-s_{q-1}))
			//                             = lambda_{lm}^T tilde A_{ij}
			log_rho = log_rho - _sur_interval[j] * (_lambda[j][mp_temp] - _lambda[j][mp_cur]);
		}

		for (int p = 0; p < _P; p++) {
			// Corresponding to nu_ij x_ijk^T beta_lm
			log_rho = log_rho + _sur_cov[p] * (_beta[p][mp_temp] - _beta[p][mp_cur]);
		}

		// Corresponding to - exp(x_{ijk}^T beta_{lm}) * ( 1- S_pch(y_ij | lambda_lm) )
		for (int k = 0; k < _n_eml; k++) {
			int_time = _int_eml[k] - 1; // the corresponding time interval for Y_{ijk}
										// cout << "int_time = " << int_time << endl;
										// exp(x_{ijk}^T beta_{lm}) 
			exp_xbeta_cur = 1.0;
			exp_xbeta_temp = 1.0;
			for (int p = 0; p < _P; p++) {
				exp_xbeta_cur = exp_xbeta_cur * exp(_X[k][p] * _beta[p][mp_cur]);
				exp_xbeta_temp = exp_xbeta_temp * exp(_X[k][p] * _beta[p][mp_temp]);
			}

			// cout << "exp_xbeta = " << exp_xbeta << endl;
			// ( 1- S_pch(y_ij | lambda_lm) )
			fail_prob_cur = 1.0;
			fail_prob_temp = 1.0;
			for (int j = 0; j < int_time; j++) {
				fail_prob_cur = fail_prob_cur * exp(-_lambda[j][mp_cur] * (_knots[j + 1] - _knots[j]));
				fail_prob_temp = fail_prob_temp * exp(-_lambda[j][mp_temp] * (_knots[j + 1] - _knots[j]));
			}

			fail_prob_cur = fail_prob_cur * exp(-_lambda[int_time][mp_cur] * (_Y[k] - _knots[int_time]));
			fail_prob_cur = 1.0 - fail_prob_cur;

			fail_prob_temp = fail_prob_temp * exp(-_lambda[int_time][mp_temp] * (_Y[k] - _knots[int_time]));
			fail_prob_temp = 1.0 - fail_prob_temp;
			// cout << "fail_prob = " << fail_prob << endl;

			log_rho = log_rho - exp_xbeta_temp * fail_prob_temp + exp_xbeta_cur * fail_prob_cur;
			// cout << "log_prob["<< i << "]["<< ind_mem << "] = " << log_prob[i][ind_mem] << endl;

		}

		delete[] event_prob;

		if (log_rho > log(_rng.runif())) {
			return mp_temp;
		}
		else {
			return mp_cur;
		}
	}

	void _update_beta(int _C, int _M, int _P, int _p, double* _knots,
		double** _sur_cov, int** _c_list,
		double** _X, double* _Y, int* _int_eml,
		double** _lambda,
		int* _mem_pair, omprng _rng, double* beta_prior, double beta_proposal,
		double** _beta) {

		int eml = 0;
		int ind_mem, n_eml, int_time;
		double exp_xbeta, fail_prob;
		double* log_rho = new double[_C * _C];
		double* beta_temp = new double[_C * _C];
		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				// sample from proposal distribution
				beta_temp[ind_mem] = _rng.rnorm(_beta[_p][ind_mem], beta_proposal);
				// prior
				if (_p == 0) {
					log_rho[ind_mem] = -pow(beta_temp[ind_mem] - beta_prior[0], 2) / 2 / beta_prior[1] + pow(_beta[_p][ind_mem] - beta_prior[0], 2) / 2 / beta_prior[1];
				}
				else {
					log_rho[ind_mem] = -pow(beta_temp[ind_mem], 2) / 2 / beta_prior[1] + pow(_beta[_p][ind_mem], 2) / 2 / beta_prior[1];
				}
			}
		}

		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;
				// if (_p == 0) {
				// 	cout << "The log_ratio for beta[" << l << "][" << m << "][" << _p << "] = " << log_rho[ind_mem] << endl;
				// }
			}
		}

		for (int i = 0; i < _M; i++) {
			n_eml = _c_list[i][2];
			ind_mem = _mem_pair[i];

			log_rho[ind_mem] = log_rho[ind_mem] + _sur_cov[i][_p] * (beta_temp[ind_mem] - _beta[_p][ind_mem]);

			for (int k = 0; k < n_eml; k++) {

				int_time = _int_eml[eml + k] - 1; // the corresponding time interval for Y_{ijk}
											// cout << "int_time = " << int_time << endl;
											// exp(x_{ijk}^T beta_{lm}) 
				exp_xbeta = 1.0;
				for (int p = 0; p < _P; p++) {
					exp_xbeta = exp_xbeta * exp(_X[eml + k][p] * _beta[p][ind_mem]);
				}

				// cout << "exp_xbeta = " << exp_xbeta << endl;
				// ( 1- S_pch(y_ij | lambda_lm) )
				fail_prob = 1.0;
				for (int j = 0; j < int_time; j++) {
					fail_prob = fail_prob * exp(-_lambda[j][ind_mem] * (_knots[j + 1] - _knots[j]));
				}
				fail_prob = fail_prob * exp(-_lambda[int_time][ind_mem] * (_Y[eml + k] - _knots[int_time]));
				fail_prob = 1.0 - fail_prob;
				// cout << "fail_prob = " << fail_prob << endl;

				log_rho[ind_mem] = log_rho[ind_mem] - exp_xbeta * fail_prob * (exp(_X[eml + k][_p] * (beta_temp[ind_mem] - _beta[_p][ind_mem])) - 1);

			}
			eml = eml + n_eml;
		}

		// update beta
		
		for (int l = 0; l < _C; l++) {
			for (int m = 0; m < _C; m++) {
				ind_mem = l * _C + m;

				// cout << "Before updating, beta[" << _p << "][" << ind_mem << "] = " << _beta[_p][ind_mem] << endl;
				// cout << "beta[" << _p << "][" << ind_mem << "] sampled from proposal distribution is " << beta_temp[ind_mem] << endl;
				if (log_rho[ind_mem] > log(_rng.runif())) {
					_beta[_p][ind_mem] = beta_temp[ind_mem];
				}
				// cout << "After updating, beta[" << _p << "][" << ind_mem << "] = " << _beta[_p][ind_mem] << endl;
			}
		}
		

		delete[] log_rho;
		delete[] beta_temp;
	}

	void _update_lambda(int _C, int _M, int _J, double* _knots,int** _c_list,
		int** _fail_ind, double* _Y, int* _int_eml,
		int* _mem_pair, int* _intention, omprng _rng, double* shape_prior, double* rate_prior,
		double** _lambda) {


		// auxiliary varaibles to update lambda
		// cout << "Assign a double pointer to lambda." << endl;
		double** shape_lam = new double*[_J];
		double** rate_lam = new double*[_J];

		int ind_mem;
		// prior
		for (int j = 0; j < _J; j++) {
			shape_lam[j] = new double[_C * _C];
			rate_lam[j] = new double[_C * _C];
			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;
					shape_lam[j][ind_mem] = shape_prior[j];
					rate_lam[j][ind_mem] = rate_prior[j];
				}
			}
		}

		// log_ratio
		int eml = 0;
		int n_eml, int_time;
		for (int i = 0; i < _M; i++) {
			n_eml = _c_list[i][2];
			ind_mem = _mem_pair[i];
			for (int j = 0; j < _J; j++) {
				shape_lam[j][ind_mem] = shape_lam[j][ind_mem] + _fail_ind[i][j];
			}
			for (int k = 0; k < n_eml; k++) {
				int_time = _int_eml[eml + k] - 1;
				for (int j = 0; j < int_time; j++) {
					rate_lam[j][ind_mem] = rate_lam[j][ind_mem] + _intention[eml + k] * (_knots[j + 1] - _knots[j]);
				}
				rate_lam[int_time][ind_mem] = rate_lam[int_time][ind_mem] + _intention[eml + k] * (_Y[eml + k] - _knots[int_time]);
			}
			
			eml = eml + n_eml;
		}

		// update lambda
		// double lam_next;
		for (int j = 0; j < _J; j++) {
			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;
					// cout << "Before updating, lambda[" << l << "][" << m << "][" << j << "] = " << _lambda[j][ind_mem] << endl;
					// cout << "The shape of lambda[" << l << "][" << m << "][" << j << "] is " << shape_lam[j][ind_mem] << ", and the rate of lambda[" << l << "][" << m << "][" << j << "] is " << rate_lam[j][ind_mem] << "." << endl;
					
					// lam_next = _rng.rgamma(shape_lam[j][ind_mem], 1.0 / rate_lam[j][ind_mem]);
					// cout << "lambda_temp = " << lam_next << endl;
					// _lambda[j][ind_mem] = lam_next;
					
					_lambda[j][ind_mem] = _rng.rgamma(shape_lam[j][ind_mem], 1.0 / rate_lam[j][ind_mem]);

					if (_lambda[j][ind_mem] < exp(-100)) {
						_lambda[j][ind_mem] = exp(-100);
					}

					// cout << "After updating, lambda[" << l << "][" << m << "][" << j << "] = " << _lambda[j][ind_mem] << endl;
				}
			}
		}

		for (int j = 0; j < _J; j++) {
			delete[] shape_lam[j];
			delete[] rate_lam[j];
		}
		delete[] shape_lam;
		delete[] rate_lam;
	}

	void _expectation_lambda(int _C, int _M, int _J, double* _knots, int** _c_list,
		int** _fail_ind, double* _Y, int* _int_eml,
		int* _mem_pair, int* _intention, double* shape_prior, double* rate_prior,
		double** _lambda) {


		// auxiliary varaibles to update lambda
		// cout << "Assign a double pointer to lambda." << endl;
		double** shape_lam = new double*[_J];
		double** rate_lam = new double*[_J];

		int ind_mem;
		// prior
		for (int j = 0; j < _J; j++) {
			shape_lam[j] = new double[_C * _C];
			rate_lam[j] = new double[_C * _C];
			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;
					shape_lam[j][ind_mem] = shape_prior[j];
					rate_lam[j][ind_mem] = rate_prior[j];
				}
			}
		}

		// log_ratio
		int eml = 0;
		int n_eml, int_time;
		for (int i = 0; i < _M; i++) {
			n_eml = _c_list[i][2];
			ind_mem = _mem_pair[i];
			for (int j = 0; j < _J; j++) {
				shape_lam[j][ind_mem] = shape_lam[j][ind_mem] + _fail_ind[i][j];
			}
			for (int k = 0; k < n_eml; k++) {
				int_time = _int_eml[eml + k] - 1;
				for (int j = 0; j < int_time; j++) {
					rate_lam[j][ind_mem] = rate_lam[j][ind_mem] + _intention[eml + k] * (_knots[j + 1] - _knots[j]);
				}
				rate_lam[int_time][ind_mem] = rate_lam[int_time][ind_mem] + _intention[eml + k] * (_Y[eml + k] - _knots[int_time]);
			}

			eml = eml + n_eml;
		}

		// calculate the posterior mean of lambda
		for (int j = 0; j < _J; j++) {
			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;
					_lambda[j][ind_mem] = shape_lam[j][ind_mem]/rate_lam[j][ind_mem];

					if (_lambda[j][ind_mem] < exp(-100)) {
						_lambda[j][ind_mem] = exp(-100);
					}
				}
			}
		}

		for (int j = 0; j < _J; j++) {
			delete[] shape_lam[j];
			delete[] rate_lam[j];
		}
		delete[] shape_lam;
		delete[] rate_lam;
	}

	void _update_pi(int _N, int _C, int _M,
		int** _c_list, double* _xi,
		int* _mem_pair, omprng _rng,
		double** _pi) {


		// auxiliary varaibles to update lambda
		// cout << "Assign a double pointer to lambda." << endl;
		double** sum_mem = new double*[_N];
		for (int s = 0; s < _N; s++) {
			sum_mem[s] = new double[_C];
			for (int l = 0; l < _C; l++) {
				sum_mem[s][l] = _xi[l];
			}
		}

		int ind_mem, l, m, sender, receiver;
		// log_ratio
		for (int i = 0; i < _M; i++) {
			sender = _c_list[i][0] - 1;
			receiver = _c_list[i][1] - 1;
			ind_mem = _mem_pair[i];
			l = ind_mem / _C;
			m = ind_mem - _C * l;

			sum_mem[sender][l] = sum_mem[sender][l] + 1;
			sum_mem[receiver][m] = sum_mem[receiver][m] + 1;
		}

		// update pi
		for (int s = 0; s < _N; s++) {
			//for (int l = 0; l < _C; l++) {
			// 	cout << "sum_mem[" << s << "][" << l << "] = " << sum_mem[s][l] << endl;
			//}
			rand_Dir(sum_mem[s], _C,
				_rng,
				_pi[s]);

			// set a lower bound of pi to avoid infinite value in updating xi
			for (l = 0; l < _C; l++) {
				if (_pi[s][l] < exp(-20.0)) {
					_pi[s][l] = exp(-20.0);
				}
			}
		}

		for (int s = 0; s < _N; s++) {
			delete[] sum_mem[s];
		}
		delete[] sum_mem;
	}
	
	void _update_xi(int _N, int _C,
		double** _pi, omprng _rng, double* xi_prior, double xi_proposal,
		double* _xi) {

		
		double alpha = 0.0;
		double alpha_temp = 0.0;
		double log_rho_alpha = 0.0;
		double* eta = new double[_C];
		double* eta_proposal = new double[_C];
		double* eta_temp = new double[_C];
		double log_rho_eta = 0.0;
		for (int l = 0; l < _C; l++) {
			alpha = alpha + _xi[l];
		}
		for (int l = 0; l < _C; l++) {
			eta[l] = _xi[l] / alpha;
			eta_proposal[l] = eta[l] * xi_proposal;
		}

		// update alpha
		alpha_temp = _rng.rbeta(alpha, _C - alpha) * _C;

		// prior
		log_rho_alpha = (xi_prior[0] - 1) * (log(alpha_temp) - log(alpha)) + (xi_prior[1] - 1) * (log(_C - alpha_temp) - log(_C - alpha));

		// likelihood
		log_rho_alpha = log_rho_alpha + _N * (lgamma(alpha_temp) - lgamma(alpha));
		for (int l = 0; l < _C; l++) {
			log_rho_alpha = log_rho_alpha + _N * (lgamma(alpha * eta[l]) - lgamma(alpha_temp * eta[l]));
			for (int s = 0; s < _N; s++) {
				log_rho_alpha = log_rho_alpha + (alpha_temp - alpha) * eta[l] * log(_pi[s][l]);
			}
		}

		// proposal
		log_rho_alpha = log_rho_alpha + (alpha_temp - 1.0) * log(alpha) + (_C - alpha_temp - 1.0) * log(_C - alpha);
		log_rho_alpha = log_rho_alpha - lgamma(alpha_temp) - lgamma(_C - alpha_temp);
		log_rho_alpha = log_rho_alpha - (alpha - 1.0) * log(alpha_temp) - (_C - alpha - 1.0) * log(_C - alpha_temp);
		log_rho_alpha = log_rho_alpha + lgamma(alpha) + lgamma(_C - alpha);

		if (log_rho_alpha > log(_rng.runif())) {
			alpha = alpha_temp;
		}

		// update eta
		rand_Dir(eta_proposal, _C, _rng,
			eta_temp);


		for (int l = 0; l < _C; l++) {
			log_rho_eta = log_rho_eta + _N * (lgamma(alpha * eta[l]) - lgamma(alpha * eta_temp[l]));
			for (int s = 0; s < _N; s++) {
				log_rho_eta = log_rho_eta + alpha * (eta_temp[l] - eta[l]) * log(_pi[s][l]);
			}
			log_rho_eta = log_rho_eta + lgamma(xi_proposal * eta[l]) + (xi_proposal * eta_temp[l] - 1.0) * log(eta[l]);
			log_rho_eta = log_rho_eta - lgamma(xi_proposal * eta_temp[l]) - (xi_proposal * eta[l] - 1.0) * log(eta_temp[l]);
		}

		if (log_rho_eta > log(_rng.runif())) {
			for (int l = 0; l < _C; l++) {
				eta[l] = eta_temp[l];
			}
		}

		for (int l = 0; l < _C; l++) {
			_xi[l] = alpha * eta[l];
		}

		delete[] eta;
		delete[] eta_proposal;
		delete[] eta_temp;
	}

	void calculate_loglike(int _C, int _N, int _M, int _J, int _P, int** _c_list, double* _knots,
		int** _fail_ind, double** _sur_interval, double** _sur_cov,
		double** _X, double* _Y, int* _nu, int* _int_eml,
		double** _lambda, double** _beta, double** _pi, double* _xi,
		double* _loglike) {

		double* log_p = new double[_C * _C];
		double* event_prob = new double[_C * _C];
		double exp_xbeta, fail_prob, max_log, sum_prob, sum_xi;
		int int_time, ind_mem, n_eml, eml, sender, receiver;

		eml = 0;
		for (int i = 0; i < _M; i++) {
			sender = _c_list[i][0] - 1;
			receiver = _c_list[i][1] - 1;
			n_eml = _c_list[i][2];

			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;

					log_p[ind_mem] = log(_pi[sender][l]) + log(_pi[receiver][m]);
					
					/*
					if (i == 4) {
						cout << "Pi part: log_prob[" << i << "][" << ind_mem << "] = " << log_p[ind_mem] << endl;
					}
					*/
					for (int j = 0; j < _J; j++) {
						// Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
						log_p[ind_mem] = log_p[ind_mem] + _fail_ind[i][j] * log(_lambda[j][ind_mem]);
						// Corresponding to nu_ij (- lambda_p(y_ij - s_{p-1}) - sum_q = 1 ^ p-1 lambda_q (s_q-s_{q-1}))
						//                             = lambda_{lm}^T tilde A_{ij}
						log_p[ind_mem] = log_p[ind_mem] - _sur_interval[i][j] * _lambda[j][ind_mem];
					}

					for (int p = 0; p < _P; p++) {
						// Corresponding to nu_ij x_ijk^T beta_lm
						log_p[ind_mem] = log_p[ind_mem] + _sur_cov[i][p] * _beta[p][ind_mem];
					}

					/*
					if (i == 4) {
						cout << "Failure time part: log_prob[" << i << "][" << ind_mem << "] = " << log_p[ind_mem] << endl;
					}
					*/

					// Corresponding to - exp(x_{ijk}^T beta_{lm}) * ( 1- S_pch(y_ij | lambda_lm) )
					for (int k = 0; k < n_eml; k++) {
						int_time = _int_eml[eml + k] - 1; // the corresponding time interval for Y_{ijk}
													// cout << "int_time = " << int_time << endl;
													// exp(x_{ijk}^T beta_{lm}) 
						exp_xbeta = 1.0;
						for (int p = 0; p < _P; p++) {
							exp_xbeta = exp_xbeta * exp(_X[eml + k][p] * _beta[p][ind_mem]);
						}

						
						// ( 1- S_pch(y_ij | lambda_lm) )
						fail_prob = 1.0;
						for (int j = 0; j < int_time; j++) {
							fail_prob = fail_prob * exp(-_lambda[j][ind_mem] * (_knots[j + 1] - _knots[j]));
						}
						fail_prob = fail_prob * exp(-_lambda[int_time][ind_mem] * (_Y[eml + k] - _knots[int_time]));
						fail_prob = 1.0 - fail_prob;
						
						/*
						if (i == 4) {
							cout << "For " << k << "-th response time, " << endl;
							for (int p = 0; p < _P; p++) {
								cout << "beta[ " << p << "][" << ind_mem << "] = " << _beta[p][ind_mem] << endl;
								cout << "X[" << k << "][" << p << "] = " << _X[eml + k][p] << endl;
							}
						
							for (int j = 0; j < int_time; j++) {
								cout << "lambda[ " << j << "][" << ind_mem << "] = " << _lambda[j][ind_mem] << endl;
							}

							cout << "exp_xbeta = " << exp_xbeta << endl;
							cout << "fail_prob = " << fail_prob << endl;
						}
						*/

						log_p[ind_mem] = log_p[ind_mem] - exp_xbeta * fail_prob;
					}
					//if (i == 4) {
					//	cout << "log_prob[" << i << "][" << ind_mem << "] = " << log_p[ind_mem] << endl;
					//}
				}
			}

			max_log = vec_max(log_p, _C * _C);
			sum_prob = 0.0;
			for (int l = 0; l < _C; l++) {
				for (int m = 0; m < _C; m++) {
					ind_mem = l * _C + m;
					event_prob[ind_mem] = exp(log_p[ind_mem] - max_log);
					sum_prob = sum_prob + event_prob[ind_mem];
				}
			}

			_loglike[0] = _loglike[0] + max_log + log(sum_prob);

			eml = eml + n_eml;
		}

		// cout << "After calculating the first part, loglike = " << _loglike[0] << endl;
		sum_xi = 0.0;
		for (int l = 0; l < _C; l++) {
			sum_xi = sum_xi + _xi[l];
		}
		_loglike[0] = _loglike[0] + _N * lgamma(sum_xi);

		for (int l = 0; l < _C; l++) {
			_loglike[0] = _loglike[0] - _N * lgamma(_xi[l]);
			for (int s = 0; s < _N; s++) {
				_loglike[0] = _loglike[0] + (_xi[l] - 1) * log(_pi[s][l]);
			}
		}
		//cout << "loglike = " << _loglike[0] << endl;
		delete[] event_prob;
		delete[] log_p;

	}


	double calculate_DIC(int _C, int _N, int _M, int _J, int _P, int** _c_list, double* _knots,
		int** _fail_ind, double** _X, double* _Y, int* _nu, int* _int_eml,
		int* _mem_pair, int* _intentions,
		double** _lambda, double** _beta, double** _pi, double* _xi) {

		double dic = 0.0;
		int int_time, ind_mem, l, m, n_eml, eml, sender, receiver;
		double xbeta, sur_prob, sum_xi;
		
		eml = 0;
		for (int i = 0; i < _M; i++) {
			sender = _c_list[i][0] - 1;
			receiver = _c_list[i][1] - 1;
			n_eml = _c_list[i][2];

			ind_mem = _mem_pair[i];
			l = ind_mem / _C;
			m = ind_mem - l * _C;

			dic = dic + log(_pi[sender][l]) + log(_pi[receiver][m]);

			/*
			if (i < 5) {
			cout << "Pi part: dic[" << i << "] = " << dic << endl;
			}
			*/

			for (int j = 0; j < _J; j++) {
				// Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
				dic = dic + _fail_ind[i][j] * log(_lambda[j][ind_mem]);
			}

			/*
			if (i < 5) {
			cout << "log lambda part: dic[" << i << "] = " << dic << endl;
			}
			*/

			for(int k = 0; k < n_eml; k++) {

				// nu_ijk * log(N_{ijk})
				if (_nu[eml + k] > 0) {
					dic = dic + _nu[eml + k] * log(_intentions[eml + k]);
				}

				// N_ijk * log(S_pch(y_ij | lambda_lm))
				int_time = _int_eml[eml + k] - 1; // the corresponding time interval for Y_{ijk}
												  // cout << "int_time = " << int_time << endl;
												  // exp(x_{ijk}^T beta_{lm}) 

				// S_pch(y_ij | lambda_lm)
				sur_prob = 1.0;
				for (int j = 0; j < int_time; j++) {
					sur_prob = sur_prob * exp(-_lambda[j][ind_mem] * (_knots[j + 1] - _knots[j]));
				}
				sur_prob = sur_prob * exp(-_lambda[int_time][ind_mem] * (_Y[eml + k] - _knots[int_time]));

				dic = dic + _intentions[eml + k] * log(sur_prob);

				// N_ijk * x_ijk * beta_lm - exp(x_ijk * beta_lm)
				xbeta = 0.0;
				for (int p = 0; p < _P; p++) {
					xbeta = xbeta + _X[eml + k][p] * _beta[p][ind_mem];
				}

				dic = dic + _intentions[eml + k] * xbeta - exp(xbeta);


				// - log(N_ijk!) = lgamma(N_ijk + 1)
				dic = dic - lgamma(_intentions[eml + k] + 1);
			}

			/*
			if (i < 5) {
				cout << "all: dic[" << i << "] = " << dic << endl;
			}
			*/

			eml = eml + n_eml;
		}

		// cout << "After calculating the first part, DIC = " << dic << endl;

		sum_xi = 0.0;
		for (int l = 0; l < _C; l++) {
			sum_xi = sum_xi + _xi[l];
		}
		dic = dic + _N * lgamma(sum_xi);

		for (int l = 0; l < _C; l++) {
			dic = dic - _N * lgamma(_xi[l]);
			for (int s = 0; s < _N; s++) {
				dic = dic + (_xi[l] - 1) * log(_pi[s][l]);
			}
		}

		return dic;
	}

	double calculate_DIC_conditional(int _C, int _N, int _M, int _J, int _P, int** _c_list, double* _knots,
		int** _fail_ind, double** _X, double* _Y, int* _nu, int* _int_eml,
		int* _mem_pair, int* _intentions,
		double** _lambda, double** _beta) {

		double dic = 0.0;
		int int_time, ind_mem, l, m, n_eml, eml, sender, receiver;
		double xbeta, log_sur_prob, sum_xi;

		eml = 0;
		for (int i = 0; i < _M; i++) {
			sender = _c_list[i][0] - 1;
			receiver = _c_list[i][1] - 1;
			n_eml = _c_list[i][2];

			ind_mem = _mem_pair[i];
			l = ind_mem / _C;
			m = ind_mem - l * _C;

			// dic = dic + log(_pi[sender][l]) + log(_pi[receiver][m]);

			/*
			if (i < 5) {
			cout << "Pi part: dic[" << i << "] = " << dic << endl;
			}
			*/

			for (int j = 0; j < _J; j++) {
				// Corresponding to nu_ijk log(lambda_{lmp_{ijk}}) = (log(lambda_{lm}))^T B_{ij}
				dic = dic + _fail_ind[i][j] * log(_lambda[j][ind_mem]);
			}

			/*
			if (std::isnan(dic)) {
				cout << "log lambda part induces nan DIC." << endl;
				cout << "At pair " << i << endl;
				return;
			}
			*/

			for (int k = 0; k < n_eml; k++) {

				// nu_ijk * log(N_{ijk})
				if (_nu[eml + k] == 1) {
					dic = dic + _nu[eml + k] * log(_intentions[eml + k]);
				}

				/*
				if (std::isnan(dic)) {
					cout << "nu log(N) part induces nan DIC." << endl;
					cout << "intentions[" << eml + k << "] = " << _intentions[eml + k] << endl;
					return;
				}
				*/

				// N_ijk * log(S_pch(y_ij | lambda_lm))
				int_time = _int_eml[eml + k] - 1; // the corresponding time interval for Y_{ijk}
												  // cout << "int_time = " << int_time << endl;
												  // exp(x_{ijk}^T beta_{lm}) 

												  // S_pch(y_ij | lambda_lm)
				log_sur_prob = 0.0;
				for (int j = 0; j < int_time; j++) {
					log_sur_prob = log_sur_prob -_lambda[j][ind_mem] * (_knots[j + 1] - _knots[j]);
				}
				log_sur_prob = log_sur_prob -_lambda[int_time][ind_mem] * (_Y[eml + k] - _knots[int_time]);

				dic = dic + _intentions[eml + k] * log_sur_prob;

				/*
				if (std::isnan(dic)) {
					cout << "N * S(y_ijk) part induces nan DIC." << endl;
					cout << "intentions[" << eml + k << "] = " << _intentions[eml + k] << endl;
					cout << "sur_prob = " << log_sur_prob << endl;
					return;
				}
				*/

				// N_ijk * x_ijk * beta_lm - exp(x_ijk * beta_lm)
				xbeta = 0.0;
				for (int p = 0; p < _P; p++) {
					xbeta = xbeta + _X[eml + k][p] * _beta[p][ind_mem];
				}

				dic = dic + _intentions[eml + k] * xbeta - exp(xbeta);

				/*
				if (std::isnan(dic)) {
					cout << "N * x beta part induces nan DIC." << endl;
					cout << "intentions[" << eml + k << "] = " << _intentions[eml + k] << endl;
					cout << "xbeta is " << xbeta << endl;
					return;
				}
				*/

				// - log(N_ijk!) = lgamma(N_ijk + 1)
				dic = dic - lgamma(_intentions[eml + k] + 1.0);

				/*
				if (std::isnan(dic)) {
					cout << "log(N!) part induces nan DIC." << endl;
					cout << "intentions[" << eml + k << "] = " << _intentions[eml + k] << endl;
					return;
				}
				*/
			}

			/*
			if (i < 5) {
			cout << "all: dic[" << i << "] = " << dic << endl;
			}
			*/

			eml = eml + n_eml;
		}

		// cout << "After calculating the first part, DIC = " << dic << endl;
		/*
		sum_xi = 0.0;
		for (int l = 0; l < _C; l++) {
			sum_xi = sum_xi + _xi[l];
		}
		dic = dic + _N * lgamma(sum_xi);

		for (int l = 0; l < _C; l++) {
			dic = dic - _N * lgamma(_xi[l]);
			for (int s = 0; s < _N; s++) {
				dic = dic + (_xi[l] - 1) * log(_pi[s][l]);
			}
		}
		*/

		return dic;
	}

	void MCMC_iteration(int* seed, int* nc, int* MCMC_iter,
		int* num_usr, int* num_mem, int* num_com, int* num_eml, int* num_int, int* num_cov,
		double* knots, int* com_list, double* sur_time, int* censor_ind, double* cov, int* int_eml,
		double* xi, double* lambda, double* beta,
		double* lambda_prior, double* lambda_prior_rate, double* beta_prior, double* xi_prior, int* update_pi, 
		double* sur_interval, int *fail_ind, double* sur_cov,
		int* mem_pair, double* pi, int* intentions, 
		double* log_like, 
		int* ind_after_burnin, int* re_iter, double* DIC) {

		// set random number generator
		// Set the seed for RNG
		omprng MCMC_Rng;
		// cout << "The seed is " << seed[0] << endl;
		MCMC_Rng.fixedSeed(seed[0]);
		// Set the number of cores for parallel
		omp_set_num_threads(nc[0]);

		// cout << "Load dim information" << endl;
		int N = num_usr[0];
		// cout << "The number of users is " << N << ";" << endl;
		int C = num_mem[0];
		// cout << "The number of membership is " << C << ";" << endl;
		int M = num_com[0];
		// cout << "The number of communicating pairs is " << M << ";" << endl;
		int N_eml = num_eml[0];
		// cout << "The number of emails is " << N_eml << ";" << endl;
		int J = num_int[0];
		// cout << "The number of intervals is " << J << ";" << endl;
		int P = num_cov[0];
		// cout << "The number of covariates is " << P << ";" << endl;

		int ind_mem, int_time, ind_res, eml;
		/////////////////////////////////////////////////////////////////////////////////////////////////
		// Allocate memory for parameters, latent variables and their auxiliary variables //
		////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////
		// Observed data //
		// Record the sender, receiver and email number between them
		int** c_list = new int*[M];
		for (int i = 0; i < M; i++) {
			c_list[i] = &com_list[3 * i];
		}

		// Record the covariates
		// cout << "Assign matrix a double pointer" << endl;
		double** X = new double*[N_eml];
		for (int eml = 0; eml < N_eml; eml++) {
			X[eml] = &cov[eml * P];
		}

		// Record survival time of each edge in J time intervals
		double** sur_interval_c = new double*[M];
		for (int i = 0; i < M; i++) {
			sur_interval_c[i] = &sur_interval[i * J];
		}

		// Record failure indicator of each edge in J time intervals
		int** fail_ind_c = new int*[M];
		for (int i = 0; i < M; i++) {
			fail_ind_c[i] = &fail_ind[i * J];
		}

		// Record the product of x_{ijk} and nu_{ijk} for each edge 
		double** sur_cov_c = new double*[M];
		for (int i = 0; i < M; i++) {
			sur_cov_c[i] = &sur_cov[i * P];
		}

		//////////////////
		// Parameter //
		 // Record the parameter  lambda
		double** lam_c = new double*[J];
		for (int j = 0; j < J; j++) {
			lam_c[j] = &lambda[j * C * C];
		}

		double* lam_rate_prior = new double[J];
		for (int j = 0; j < J; j++) {
			lam_rate_prior[j] = 1.0 / lambda_prior_rate[0];
		}

		// Record the parameter beta
		double** beta_c = new double*[P];
		for (int p = 0; p < P; p++) {
			beta_c[p] = &beta[p * C * C];
		}

		// Record beta in the Newton-Raphson algorithm
		double beta_temp;

		// Record the number of xi
		double sum_xi = 0.0;
		for (int l = 0; l < C; l++) {
			sum_xi = sum_xi + xi[l];
		}

		////////////////////////
		// Latent variables //
		// cout << "Assign a double pointer to gamma." << endl;
		double** pi_c = new double*[N];
		for (int s = 0; s < N; s++) {
			pi_c[s] = &pi[C * s];
		}

		///////////////////////////
		// 2. MCMC sampling //
		///////////////////////////
		// auxiliary variable
		// cout << "Start MCMC sampling." << endl;

		///////////////////////////////////////////////
		//  1) update N_{ijk}  and (l_{ij}, m_{ij}) //
		eml = 0;
		// cout << "Update latent variables." << endl;

		int sender, receiver, n_eml;
		double exp_xbeta, fail_prob, lam_poi;
#pragma omp parallel for
		for (int i = 0; i < M; i++) {
			
			sender = c_list[i][0] - 1;
			receiver = c_list[i][1] - 1;
			n_eml = c_list[i][2];

		
			_update_N(P, n_eml, knots,
				&(X[eml]), &(sur_time[eml]), &(censor_ind[eml]), &(int_eml[eml]),
				lam_c, beta_c,
				mem_pair[i], MCMC_Rng,
				&(intentions[eml]));

			mem_pair[i] = _update_mp(C, J, P, n_eml, knots,
				fail_ind_c[i], sur_interval_c[i], sur_cov_c[i],
				&(X[eml]), &(sur_time[eml]), &(censor_ind[eml]), &(int_eml[eml]),
				lam_c, beta_c, pi_c[sender], pi_c[receiver],
				MCMC_Rng);

			eml = eml + n_eml;
		}
		// cout << "The log-likelihood is equal to " << log_like[0] << endl;
		
		//////////////////////////////
		//  2) update beta_{lmr}  //
		// cout << "Update beta." << endl;
		for (int p = 0; p < P; p++) {
			_update_beta(C, M, P, p, knots,
				sur_cov_c, c_list,
				X, sur_time, int_eml,
				lam_c,
				mem_pair, MCMC_Rng, beta_prior, 0.1,
				beta_c);
		}
		
		///////////////////////////////////
		//  3) update lambda_{lmp}  //
		// cout << "Update lambda." << endl;
		_update_lambda(C, M, J, knots,c_list,
			fail_ind_c, sur_time, int_eml,
			mem_pair, intentions, MCMC_Rng, lambda_prior, lam_rate_prior,
			lam_c);

		/////////////////////////
		//  4) update pi_{il}  //
		// cout << "Update pi." << endl;
		_update_pi(N, C, M,
			c_list, xi,
			mem_pair, MCMC_Rng,
			pi_c);
		

		/////////////////////////
		//  5) update xi_{l}  //
		// cout << "update_pi = " << update_pi[0] << endl;
		if (update_pi[0] == 1 & C > 1) {
			// cout << "Update xi." << endl;
			_update_xi(N, C,
				pi_c, MCMC_Rng, xi_prior, 10.0,
				xi);
		}

		//////////////////////////////////
		//  6) calculate loglikelihood  //

		// cout << "Calculate log-likelihood function." << endl;
		calculate_loglike(C, N, M, J, P, c_list, knots,
			fail_ind_c, sur_interval_c, sur_cov_c,
			X, sur_time, censor_ind, int_eml,
			lam_c, beta_c, pi_c, xi,
			log_like);

		/*
		double DIC_test;
		DIC_test = calculate_DIC(C, N, M, J, P, c_list, knots,
			fail_ind_c, X, sur_time, censor_ind, int_eml,
			mem_pair, intentions,
			lam_c, beta_c, pi_c, xi);
		cout << "The DIC by posterior sampling is " << DIC_test << endl;
		*/

		////////////////////////
		//  7) calculate DIC  //
		/*
		// Output the extra posterior sampling to debug DIC
		// file variable to output the posterior sampling
		ofstream post_File;

		// result/simulation_v0_r1
		string output_dir = "./extra_MCMC_sampling/";
		int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		if (!check)
			cout << "Directory " << output_dir << " is created." << endl;
		else {
			cout << "Directory " << output_dir << " already exists." << endl;
		}
		
		string out_file;
		*/

		if (ind_after_burnin[0] == 1) {
			// The posterior expected value of the joint deivance
			double DIC_post;
			int* intentions_re = new int[N_eml];
			int* mp_pair_re = new int[M];
			double** pi_re = new double*[N];
			for (int s = 0; s < N; s++) {
				pi_re[s] = new double[C];
			}

			for (eml = 0; eml < N_eml; eml++) {
				intentions_re[eml] = intentions[eml];
			}

			for (int i = 0; i < M; i++) {
				mp_pair_re[i] = mem_pair[i];
			}

			for (int s = 0; s < N; s++) {
				for (int l = 0; l < C; l++) {
					pi_re[s][l] = pi_c[s][l];
				}
			}
			

			for (int iter = 0; iter < re_iter[0]; iter++) {
				///////////////////////////////////////////////
				//  1) update N_{ijk}  and (l_{ij}, m_{ij}) //
				eml = 0;
				// cout << "Update latent variables." << endl;
				for (int i = 0; i < M; i++) {

					sender = c_list[i][0] - 1;
					receiver = c_list[i][1] - 1;
					n_eml = c_list[i][2];


					_update_N(P, n_eml, knots,
						&(X[eml]), &(sur_time[eml]), &(censor_ind[eml]), &(int_eml[eml]),
						lam_c, beta_c,
						mp_pair_re[i], MCMC_Rng,
						&(intentions_re[eml]));

					mp_pair_re[i] = _update_mp(C, J, P, n_eml, knots,
						fail_ind_c[i], sur_interval_c[i], sur_cov_c[i],
						&(X[eml]), &(sur_time[eml]), &(censor_ind[eml]), &(int_eml[eml]),
						lam_c, beta_c, pi_re[sender], pi_re[receiver],
						MCMC_Rng);

					eml = eml + n_eml;
				}

				_update_pi(N, C, M,
					c_list, xi,
					mp_pair_re, MCMC_Rng,
					pi_re);

				DIC_post = calculate_DIC_conditional(C, N, M, J, P, c_list, knots,
					fail_ind_c, X, sur_time, censor_ind, int_eml,
					mp_pair_re, intentions_re,
					lam_c, beta_c);

				DIC[0] = DIC[0] + DIC_post;
				
				// cout << "Writing posterior sampling into the directory " << output_dir << endl;
				/*
				out_file = output_dir + "N_post_iter" + to_string(MCMC_iter[0]) + ".txt";
				// cout << "Writing alpha_post into " << out_file << endl;
				post_File.open(out_file.c_str(), ios::out | ios::app);
				for (eml = 0; eml < N_eml; eml++) {
					post_File << intentions_re[eml];
					post_File << " ";
				}
				post_File << endl;
				post_File.close();

				out_file = output_dir + "mp_post_iter" + to_string(MCMC_iter[0]) + ".txt";
				// cout << "Writing alpha_post into " << out_file << endl;
				post_File.open(out_file.c_str(), ios::out | ios::app);

				for (int i = 0; i < M; i++) {
					post_File << mp_pair_re[i];
					post_File << " ";
				}
				post_File << endl;
				post_File.close();

				out_file = output_dir + "pi_post_iter" + to_string(MCMC_iter[0]) + ".txt";
				// cout << "Writing alpha_post into " << out_file << endl;
				post_File.open(out_file.c_str(), ios::out | ios::app);

				for (int s = 0; s < N; s++) {
					for (int l = 0; l < C; l++) {
						post_File << pi_re[s][l];
						post_File << " ";
					}
				}
				post_File << endl;
				post_File.close();
				*/
			}

			DIC[0] = DIC[0] / re_iter[0];

			// cout << "At " << MCMC_iter[0] << "-th resampling, the posterior expected DIC is " << DIC[0] << endl;

			// calculate the fixed point deviance
			// Calculate E{lambda|D}
			double** lam_expect = new double*[J];
			for (int j = 0; j < J; j++) {
				lam_expect[j] = new double[C * C];
			}

			// calculate the posterior mean of lambda
			_expectation_lambda(C, M, J, knots, c_list,
				fail_ind_c, sur_time, int_eml,
				mem_pair, intentions, lambda_prior, lam_rate_prior,
				lam_expect);

			/*
			for (int j = 0; j < J; j++) {
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					cout << "lam_expect[" << j << "][" << ind_mem << "] = " << lam_expect[j][ind_mem] << endl;
				}
			}
			*/

			// estimate E{lambda|D} by sampling
			double** beta_expect = new double*[P];
			for (int p = 0; p < P; p++) {
				beta_expect[p] = new double[C * C];
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					beta_expect[p][ind_mem] = 0.0;
				}
			}

			double** beta_temp = new double*[P];
			for (int p = 0; p < P; p++) {
				beta_temp[p] = new double[C * C];
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					beta_temp[p][ind_mem] = beta_c[p][ind_mem];
				}
			}

			// double* xi_expect = new double[C];
			// double* xi_temp = new double[C];
			/*
			for (int l = 0; l < C; l++) {
				xi_expect[l] = 0.0;
				xi_temp[l] = xi[l];
			}
			*/

			for (int iter = 0; iter < re_iter[0]; iter++) {
				for (int p = 0; p < P; p++) {
					_update_beta(C, M, P, p, knots,
						sur_cov_c, c_list,
						X, sur_time, int_eml,
						lam_expect,
						mem_pair, MCMC_Rng, beta_prior, 0.1,
						beta_temp);
				}
				/*
				if (C > 1) {
					// cout << "Update xi." << endl;
					_update_xi(N, C,
						pi_c, MCMC_Rng, xi_prior, 10.0,
						xi_temp);
				}
				*/
				for (int p = 0; p < P; p++) {
					for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
						beta_expect[p][ind_mem] = beta_expect[p][ind_mem] + beta_temp[p][ind_mem];
					}
				}
				/*
				for (int l = 0; l < C; l++) {
					xi_expect[l] = xi_expect[l] + xi_temp[l];
				}
				*/
			}

			// posterior mean of beta
			for (int p = 0; p < P; p++) {
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					beta_expect[p][ind_mem] = beta_expect[p][ind_mem] / re_iter[0];
				}
			}

			//for (int l = 0; l < C; l++) {
			//	xi_expect[l] = xi_expect[l] / re_iter[0];
			//}

			DIC[1] = calculate_DIC_conditional(C, N, M, J, P, c_list, knots,
				fail_ind_c, X, sur_time, censor_ind, int_eml,
				mem_pair, intentions,
				lam_expect, beta_expect);

			/*
			cout << "At " << MCMC_iter[0] << "-th resampling, The fixed point DIC is " << DIC[1] << endl;

			out_file = output_dir + "lam_exp_iter" + to_string(MCMC_iter[0]) + ".txt";
			// cout << "Writing alpha_post into " << out_file << endl;
			post_File.open(out_file.c_str(), ios::out | ios::app);
			for (int j = 0; j < J; j++) {
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					post_File << lam_expect[j][ind_mem];
					post_File << " ";
				}
				post_File << endl;
			}
			post_File << endl;
			post_File.close();

			out_file = output_dir + "beta_exp_iter" + to_string(MCMC_iter[0]) + ".txt";
			// cout << "Writing alpha_post into " << out_file << endl;
			post_File.open(out_file.c_str(), ios::out | ios::app);
			for (int p = 0; p < P; p++) {
				for (ind_mem = 0; ind_mem < C * C; ind_mem++) {
					post_File << beta_expect[p][ind_mem];
					post_File << " ";
				}
				post_File << endl;
			}
			post_File << endl;
			post_File.close();
			*/
			

			delete[] intentions_re;
			delete[]  mp_pair_re;
			for (int s = 0; s < N; s++) {
				delete[] pi_re[s];
			}
			delete[] pi_re;

			for (int j = 0; j < J; j++) {
				delete[] lam_expect[j];
			}
			delete[] lam_expect;

			for (int p = 0; p < P; p++) {
				delete[] beta_expect[p];
			}
			delete[] beta_expect;

			for (int p = 0; p < P; p++) {
				delete[] beta_temp[p];
			}
			delete[] beta_temp;

			// delete[] xi_expect;
			// delete[] xi_temp;

		}
		
		delete[] c_list;
		delete[] X;
		delete[] sur_interval_c;
		delete[] fail_ind_c;
		delete[] sur_cov_c;
		delete[] lam_c;
		delete[] lam_rate_prior;
		delete[] beta_c;
		delete[] pi_c;
		
	}

}