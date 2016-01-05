#ifndef ENGINE_H
#define ENGINE_H

#include<stdio.h>
#include<stdlib.h>

#include"Arguments.h"
#include"Model.h"


#define MARK 9e99
// Substitutes la += log(1+exp(lb-la)) <=> a = a + b.
#define LOGADD(la,lb) if(la>8e99||lb-la>100) la=lb;else la+=log(1+exp(lb-la));

using namespace std;


//----------------
// Fast Möbius transform routines.
//
void sub_fumt(int j, int d, int S, int overlap, double *t, double *s, int n, int k){
	if (d < n && overlap < k){
		sub_fumt(j, d+1, S, overlap, t, s, n, k);
		S |= (1 << d);
		overlap += (d >= j+1);
		sub_fumt(j, d+1, S, overlap, t, s, n, k); 
	}
	else{
		s[S] = MARK;
		int jinS = ((S >> j) & 1);
		if (overlap + jinS <= k) s[S] = t[S];
		if (jinS) LOGADD(s[S], t[S - (1 << j)]);
	}
}

// Fast upward Möbius transform:
// s(S) := \sum_{T \subseteq S : |T| \leq k} t(T).
void fumt(double *t, double *s, int n, int k){
	double *tmp, *sprev = new double[1 << n];
	for (int T = 0; T < (1<<n); T ++) sprev[T] = t[T];
	for (int j = 0; j < n; j ++){
		sub_fumt(j, 0, 0, 0, sprev, s, n, k);
		tmp = sprev; sprev = s; s = tmp;
	}
	if (n % 2 == 0){ for (int S = 0; S < (1<<n); S ++) s[S] = sprev[S];}
	else{ tmp = sprev; sprev = s; s = tmp;}
	delete [] sprev;
}
void test_fumt(){
	int n = 4;
	double t[16], s[16];
	for (int S = 0; S < 16; S ++){
		t[S] = log(S+1);
	}
	fumt(t, s, n, 2);
	for (int S = 0; S < 16; S ++){
		cerr<<" s["<<S<<"] = "<<exp(s[S])<<endl;
	}
}



//
void sub_fdmt(int j, int d, int T, int overlap, double *s, double *t, int n, int k){
	//cerr<<" T = "<<T<<":"; print_set(cerr, T); 
	//cerr<<" overlap = "<<overlap<<endl; 
	if (d < n && overlap <= k){
		sub_fdmt(j, d+1, T, overlap, s, t, n, k);
		T |= (1 << d);
		overlap += (d <= j);
		sub_fdmt(j, d+1, T, overlap, s, t, n, k); 
	}
	else if (d == n){
		t[T] = s[T];
		int jinT = ((T >> j) & 1);
		if (jinT == 0) LOGADD(t[T], s[T + (1 << j)]);
		if (overlap > k) t[T] = MARK;
	}
}
// Fast downward Möbius transform:
// t(T) := \sum_{T \subseteq S} s(S).
void fdmt(double *s, double *t, int n, int k){
	double *tmp, *tprev = new double[1 << n];
	for (int S = 0; S < (1<<n); S ++) tprev[S] = s[S];
	for (int j = 0; j < n; j ++){
		
		sub_fdmt(j, 0, 0, 0, tprev, t, n, k);	
		tmp = t; t = tprev; tprev = tmp;
	}
	if (n % 2 == 0){ for (int T = 0; T < (1<<n); T ++) t[T] = tprev[T];}
	else{ tmp = t; t = tprev; tprev = tmp;}
	delete [] tprev;
}
void test_fdmt(){
	int n = 4;
	double t[16], s[16];
	for (int S = 0; S < 16; S ++){
		s[S] = log(S+1);
	}
	fdmt(s, t, n, 3);
	for (int T = 0; T < 16; T ++){
		cerr<<" t["<<T<<"] = "<<exp(t[T])<<endl;
	}
}

//----------------






// NOTE: Engine only calls a veru limited set of interface functions
// implemented in Model.
class Engine {
public:
	Engine(){}
	~Engine(){}	
	
	void init(Model *m){
		model = m;

		//test();
	
	}
	void test(){
		vector<int> T;
		T.push_back(1); T.push_back(2);
		
		cerr<<" log p(x_0 | x_1,2): "<<model->log_lcp(0, T)<<endl;
	
		for (int h = 0; h < model->num_layers(); h ++){
			int *Vh, *Vu, *Vl;
			int nh, nu, nl;
			model->layer(h, &Vh, &nh);	
			model->upper_layers(h, &Vu, &nu);	
			model->lower_layers(h, &Vl, &nl);
			cerr<<"Layer "<<h<<":";
			cerr<<endl<<" layer: ";
			print_nodes(cerr, Vh, nh);
			cerr<<endl<<" upper: ";
			print_nodes(cerr, Vu, nu);
			cerr<<endl<<" lower: ";
			print_nodes(cerr, Vl, nl);
			cerr<<endl;	
		}
		
		//test_fumt();
		test_fdmt(); 
		exit(1);

	}
	// Computes all edge probabilities.
	void compute_edge_probabilities(){
		int l = model->num_layers();
		// Separate computations for each layer. 
		for (int h = 0; h < l; h ++){
			if (model->edges_within(h))
				compute_edge_probabilities(h);
		}	
	}
	// Computes probabilities for edges pointing to the h-th layer.
	void compute_edge_probabilities(int h){
		cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;
		
		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);
		
		beta = new double*[nh];
		alpha = new double*[nh];
		
		k = model->max_indegree();
		
		cerr<<" . "<<nh<<" nodes"<<endl;
		
		for (int j = 0; j < nh; j ++){
			beta[j] = new double[1 << nh];
			//int i = Vh[j];
			compute_beta(j, Vh, nh, Vu, nu);
			alpha[j] = new double[1 << nh];	
			compute_alpha(j, Vh, nh);
		
			cerr<<" . Tables beta and alpha computed for node "<<j<<"."<<endl; 
		
		}
		
		cerr<<" . Tables beta and alpha are now ready."<<endl;

			
		// Step 3: Compute g_forward[];
		gf = new double[1 << nh];
		
		compute_g_forward(nh);
		
		//cerr<<" . gf[all] = "<<gf[(1<<nh)-1]<<endl;
		
		// Step 4: Compute g_backward[];
		
		gb = new double[1 << nh];
		
		compute_g_backward(nh);
		
		//cerr<<" . gb[all] = "<<gb[(1<<nh)-1]<<endl;
			
		// Check point: g_forward[] == g_backward[];
	
		// Step 5: Loop over end-point nodes j.
		double *gamma = new double[1 << nh];
		double *gfb = new double[1 << nh];
		for (int j = 0; j < nh; j ++){
			// Compute gfb[].
			for (int S = 0; S < (1<<nh); S ++){
				if ((S>>j) & 1){
					gfb[S] = -MARK;
				}
				else{
					int cS = (1 << nh) - 1 - S - (1 << j);
					gfb[S] = gf[S] + gb[cS];
					//cerr<<" j = "<<j<<": ";
					//print_set(cerr, S);
					//cerr<<"; ";
					//print_set(cerr, cS);
					//cerr<<": "<<gfb[S];
					//cerr<<endl;
				}
			}
			// Compute gamma_j() := gamma[].
			fdmt(gfb, gamma, nh, k);
			
			//cerr<<" gamma[0]: "<<gamma[0]<<endl;
			//cerr<<" gfb[0]: "<<gfb[0]<<endl;
			
			cerr<<" . Incoming edges for node "<<j<<":";
			// Loop over start-point nodes i.
			for (int i = 0; i < nh; i ++){
				if (i == j) continue;
				cerr<<" "<<i;		
				double log_prob = eval_edge(i,  j, gamma, beta[j], nh, k);
		
				model->print_edge_prob(
					cout, Vh[i], Vh[j], exp(log_prob - gb[(1<<nh)-1]));
			
			}
			cerr<<endl;
			
		} 
		delete [] gamma;
		delete [] gfb;
	
		delete [] gf; delete [] gb;
		for (int j = 0; j < nh; j ++){
			delete [] beta[j]; delete [] alpha[j];
		}
		delete [] beta; delete [] alpha;
		cerr<<" Edge probabilities now computed."<<endl; 
	}
	// beta[j][T] := sum_W gamma[j][T \cup W] .
	void compute_beta(int j, int *Vh, int nh, int *Vu, int nu){
		//cerr<<"Compute beta, j = "<<j<<", i = "<< Vh[j];
		//cerr<<", k = "<<k<<", nh = "<<nh<<":"<<endl;
		
		sub_init(0, 0, beta[j], MARK, 0, nh);
		//cerr<<" Initialization done."<<endl;

		int *S = new int[k];
		sub_beta(-1, nh, beta[j], 0, S, Vh[j], Vu, nu, 0);
		delete [] S;
	}
	void sub_beta(int jprev, int nh, double *b, int d, int *S, int i, int *Vu, int nu, int T){
		//cerr<<" S:"; print_nodes(cerr, S, d); 
		//cerr<<"; T:"; print_nodes(cerr, T, Vu, nu); cerr<<endl;
		double lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
		LOGADD(b[T], lb);
		if (d < k){
			for (int j = jprev + 1; j < nu; j ++){	
				if (Vu[j] != i){
					S[d] = Vu[j];
					if (j < nh) sub_beta(j, nh, b, d+1, S, i, Vu, nu, T | (1 << j));
					else sub_beta(j, nh, b, d+1, S, i, Vu, nu, T);
				}
			}
		}		
	}
	void sub_init(int d, int S, double *a, double value, int ones, int nh){
		//cerr<<" d: "<<d<<endl; 
		a[S] = value;
		if (d < nh){
			sub_init(d+1, S, a, value, ones, nh);
			if (ones < k) sub_init(d+1, S | (1 << d), a, value, ones+1, nh);
		}
	}
	// alpha[j][S] = sum_T beta[j][T]  (...times the prior).
	void compute_alpha(int j, int *Vh, int nh){
		// Use fast truncated upward Möbius transform.
		fumt(beta[j], alpha[j], nh, k); 
	}
	void compute_g_forward(int nh){
	
		sub_gf(0, nh, 0);
	
	}
	void sub_gf(int d, int nh, int S){
		if (d < nh){
			sub_gf(d+1, nh, S);
			sub_gf(d+1, nh, S | (1 << d)); 
		}
		else {
			double sum = MARK; 
			int T = S, J = 1;
			for (int j = 0; j < nh; j ++){
				if (T & 1){
					// Now S - J is a subset of S.
					double w = alpha[j][S - J] + gf[S - J];
					LOGADD(sum, w);
				} 		
				T >>= 1; J <<= 1;
			}
			if (S == 0) gf[0] = 0; else gf[S] = sum;
		}
	}
	void compute_g_backward(int nh){
	
		sub_gb(0, nh, 0);
	
	}
	void sub_gb(int d, int nh, int S){
		if (d < nh){
			sub_gb(d+1, nh, S);
			sub_gb(d+1, nh, S | (1 << d)); 
		}
		else {
			double sum = MARK; 
			int T = S, J = 1, complS = (1 << nh) - 1 - S;
			for (int j = 0; j < nh; j ++){
				if (T & 1){
					// Now S - J is a subset of S.
					double w = alpha[j][complS] + gb[S - J];
					LOGADD(sum, w);
				} 		
				T >>= 1; J <<= 1;
			}
			if (S == 0) gb[0] = 0; else gb[S] = sum;
		}
	}
	double eval_edge(int i, int j, double *a, double *b, int nh, int k){
		int T = 1 << i;  
		return sub_eval_edge(-1, 1, T, i, j, a, b, nh, k);
	}
	double sub_eval_edge(
		int tprev, int d, int T, int i, int j, double *a, double *b, int nh, int k){
		
		//cerr<<" T:"; print_set(cerr, T); 
		//cerr<<": "<<a[T]<<"; "<<b[T]<<endl;
			
		if (d == k){
			return a[T] + b[T];
		}
		double sum = a[T] + b[T];
		for (int t = tprev + 1; t < nh; t ++){
			if (t == i || t == j) continue;
			int Tnext = T | (1 << t);
			double w = sub_eval_edge(t, d+1, Tnext, i, j, a, b, nh, k);
			LOGADD(sum, w);
		}
		return sum;
	}
	
	Model *model;
	
	double **alpha, **beta;
	double *gf, *gb;
	int k;
};



#endif
