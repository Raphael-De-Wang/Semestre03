#ifndef MODEL_H
#define MODEL_H

#include<stdio.h>
#include<stdlib.h>

#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<math.h>

#include"Arguments.h"


using namespace std;

#define MAX_COUNT 2*8192

static double lr[16][MAX_COUNT]; 
static double lr2[16][MAX_COUNT]; 


// log Gamma(1/a + z) / Gamma(1/a) if a > 0;
// log Gamma(-a + z) / Gamma(-a) otherwise.
double log_gammaratio(int a, int z){
	if (z < MAX_COUNT && a > 0) return lr[a][z];
	if (z < MAX_COUNT && a < 0) return lr2[-a][z];
	cerr<<" Counts larger than "<<MAX_COUNT<<
		" occured. Increase MAX_COUNT. Exit now.\n"; 
	exit(1);  
}
double log_gammaratio(double z){
	return 0.50 * log(M_PI) + (z - 0.50) * log(z) - z + 1.0 / (12 * z);  
}

template<class T> void printvec(ostream& f, const vector<T>& v){
	int size= v.size();
	int i;
	for(i=0; i<size; i++){
		f<<v.at(i)<<" ";
	}
	f<<endl;
}
template<class T> void printvecs(ostream& f, const vector<vector<T> >& v){
	int size= v.size();
	int i;
	for(i=0; i<size; i++){
		printvec(f, v.at(i));
	}
}
void print_nodes(ostream& f, int *V, int n){
	for (int j = 0; j < n; j ++){
		f<<" "<<V[j];
	}
}
void print_nodes(ostream& f, int T, int *V, int n){
	int j = 0; 
	while (T > 0){
		if (T & 1){
			f<<" "<<V[j];
		}
		j ++;
		T >>= 1;
	}
}
void print_set(ostream& f, int T){
	int j = 0; 
	while (T > 0){
		if (T & 1){
			f<<" "<<j;
		}
		j ++;
		T >>= 1;
	}
}


// Reads and stores data.
class Data{
public:
	Data(){}
	~Data(){}
	void init(){
		read_data();
		downcode();
	}
	void read_data(){
		ifstream ifs(Arguments::datafile, ios::in);
		if (!ifs){
			fprintf(stderr, " Cannot read file %s.\n", Arguments::datafile);
			exit(1);
		}
		fprintf(stderr, " Reading file %s...\n", Arguments::datafile);
		
		char* buffer= new char[10000];
		ifs.getline(buffer, 10000);
		char* pch= strtok(buffer," \t");
		// http://www.cplusplus.com/ref/cstring/strtok.html
		while (pch != NULL){
			string tempstring(pch);
			heads.push_back(tempstring);
			pch = strtok(NULL, " \t");
        }
		numattributes = heads.size();
	
		fprintf(stderr, " Heading read: %d attributes.\n", numattributes);
	
		dm.clear();
		vector<int> temp;
		while (true){
			temp.clear();
			ifs.getline(buffer, 10000);
			pch = strtok(buffer," \t");
			while (pch != NULL){
				temp.push_back(atoi(pch));
				pch = strtok(NULL, " \t");
			}
			if ((int)temp.size() != numattributes) break;
			
			dm.push_back(temp);
			
			if ((int)dm.size() >= atoi(Arguments::maxnumrecords)) break;
		}
		numrecords = dm.size();
		fprintf(stderr, " Data read: %d lines.\n", numrecords); 
		delete [] buffer;
	}
	void downcode(){
		int inuse[4096][256];
		for (int v = 0; v < 255; v ++){
			for (int i = 0; i < numattributes; i ++){
				inuse[v][i] = -1;
			}
		}
		for (int t = 0; t < numrecords; t ++){
			for (int i = 0; i < numattributes; i ++){
				inuse[dm[t][i]][i] = 1;
			}
		}
		arities.clear();
		maxarity = 0;
		for (int i = 0; i < numattributes; i ++){
			int count = 0;
			for (int v = 0; v < 255; v ++){
				if (inuse[v][i] == 1){
					inuse[v][i] = count ++;
				}
			}
			arities.push_back(count);	
			if (count > maxarity) maxarity = count;
		}
		for (int t = 0; t < numrecords; t ++){
			for (int i = 0; i < numattributes; i ++){
				dm[t][i] = inuse[dm[t][i]][i];
			}
		}
		//print_data(cerr);
	}
	int get_index(string s){
		int i = 0;
		while (i < (int)heads.size() && heads[i] != s){ i ++; }
		return i;
	}
	void print_data(ostream & f){
		for (int i = 0; i < numattributes; i ++){
			f << " " << heads[i];
		}
		f << endl; 
		for (int t = 0; t < numrecords; t ++){
			for (int i = 0; i < numattributes; i ++){
				f << " " << dm[t][i];
			}
			f << endl; 
		}
	}
		
	vector<string> heads;
	vector< vector<int> > dm;
	vector<int> arities;
	int numattributes;
	int numrecords;
	int maxarity;
};


// NOTE: We assume that children come first, that is, 
// V[0] is the grandest child and V[n-1] is the grandest parent.
// Layers respectively from 0 to numlayers.
// Note that this is reverse to the input ordering.
class Layering{
public:
	Layering(){}
	~Layering(){ delete [] V;}

	void init(Data & data){
		set_layers(data);
		
		print_layers();
	
	}
	void set_layers(Data & data){
	
		V = new int[data.numattributes];
		int j = 0;
		cnh.clear(); nh.clear(); cnh.clear();
		
		if (Arguments::layeringfile[0] == '%'){
			numlayers = 1;
			int layersize = data.numattributes;
			nh.push_back(layersize);
			cnh.push_back(j);
			for (int l = 0; l < layersize; l ++){
				int i = j;
				V[j++] = i;
			}
			
			edgeswithin.clear();
			edgeswithin.push_back(true);
			return;
		}
	
		parselayeringfile(Arguments::layeringfile, true);
		numlayers = nodesinlayers.size();
	
		for (int h = numlayers-1; h >= 0; h --){
			int layersize = nodesinlayers[h].size();
			nh.push_back(layersize);
			cnh.push_back(j);
			for (int l = 0; l < layersize; l ++){
				int i = data.get_index(nodesinlayers[h][l]);
				V[j++] = i;
			}
		}
	
	}
	void print_layers(){
		for (int h = 0; h < numlayers; h ++){
			cerr<<"Layer "<<h<<":";
			for (int j = 0; j < nh[h]; j ++){
				cerr<<" "<<V[cnh[h] + j];
			}
			cerr<<endl;
		}
	}
	bool edges_within(int h){
		return edgeswithin[numlayers - h - 1];
	}	
	void parselayeringfile(char *layerfilename, bool verbose){
	
		ifstream ifs(layerfilename, ios::in);
		if( !ifs){
			cerr<<"ERROR: The file '"<<layerfilename<<"' could not be opened\n";
			cerr<<"Exiting\n";
			exit(1);
		}
		if(verbose){
			cerr<<"Reading layer file: '"<<layerfilename<<"'\n";
		}

		char* buffer= new char[10000];
		char* pch;
		int i;
		char tempchar;
		while(true){
			ifs.getline(buffer,10000);
			pch = strtok(buffer, " \t");
			vector<string> tempstringvec;	//will hold the tokenized line
			while(pch != NULL){
				string tempstring(pch);
				tempstringvec.push_back(tempstring);
				pch = strtok(NULL, " \t");
			}

			if(tempstringvec.size() < 3){
				cerr<<"   ERROR: in reading file '"<<layerfilename<<"'"<<endl;
				cerr<<"          each line must be of the form:"<<endl;
				cerr<<"          layername\t1/0\tn1 ... nn"<<endl;
				cerr<<"          EXITING"<<endl;
				exit(1);
			}

			layernames.push_back(tempstringvec.at(0));
			if(tempstringvec.at(1) == "0"){
				edgeswithin.push_back(0);
			}
			else if( tempstringvec.at(1) == "1"){
				edgeswithin.push_back(1);
			}
			else{
				cerr<<"   ERROR: in reading file '"<<layerfilename<<"'"<<endl;
				cerr<<"          the values in the second column can either be '1' or '0'"<<endl;
				cerr<<"          EXITING"<<endl;
				exit(1);
			}

			vector<string> tempnodesinlayer;
			for(i = 2; i < (int)tempstringvec.size();i ++){
				tempnodesinlayer.push_back(tempstringvec.at(i));
			}
			nodesinlayers.push_back(tempnodesinlayer);

			ifs>>tempchar;
			ifs.putback(tempchar);
			if(ifs.eof()){
				break;
			}
		}
		delete [] buffer;
		buffer = 0;
		pch = 0;
		cerr.width(13);
		cerr<<"Layer Name";
		cerr.width(24);
		cerr<<"Edges Within Allowed";
		cerr.width(14);
		cerr<<"Node Names"<<endl;
		for(i = 0; i < (int)layernames.size(); i ++){
			cerr.width(13);
			cerr<<layernames.at(i);
			cerr.width(15);
			if(edgeswithin.at(i)){
				cerr<<"true";
			}
			else{
				cerr<<"false";
			}
			cerr<<"             ";
			printvec(cerr, nodesinlayers.at(i));
		}
	}


	vector<string> layernames; 
	vector<bool> edgeswithin;
	vector<vector<string> > nodesinlayers;
	int *V;
	int numlayers;
	vector<int> nh;
	vector<int> cnh;
};







// Blocks for computing sufficient statistics.
class Tnode{
public:
	double pseudocount;
	short int arity;
	short int count;
	short int depth;
	Tnode** children;
	Tnode(const int a, const int d){
		pseudocount = 1;
		arity = a;
		count = 0;
		depth = d;
		children = new Tnode*[arity];
		for(int v = 0; v < arity; v ++){
			children[v] = NULL;
		}
	}
	~Tnode(){
		for(int v = 0; v < arity; v ++){
			if (children[v] != NULL){
				delete children[v];
			}
		}
		delete children;
		//free_memory();
	}
	void free_memory(){
		//cerr<<" Deleting d="<<depth<<" c="<<count<<endl;
		
		for(int v = 0; v < arity; v ++){
			if (children[v] != NULL){
				children[v]->free_memory();
				delete children[v];
			}
		}
		//delete [] children;
	}
	// Note: Actually S is not needed.
	void update(vector<int> & S, vector<int> & u){
		if (u.size() > 0){	
			//int j = S[u.size()-1];
			int v = u[u.size()-1]; u.pop_back();
			if (children[v] == NULL){
				//cerr<<" Creating d="<<depth<<" c="<<count<<endl;
				children[v] = new Tnode(arity, depth+1);
			}
			children[v]->update(S, u);
		}
		count ++;
	}
	void print(ostream & f){
		f << "(" << count << ")" << endl;
		for (int v = 0; v < arity; v ++){
			if (children[v] != NULL){
				for (int d = 0; d < depth; d ++) f << " "; 
				f << " " << v << " ";
				children[v]->print(f);
			}
		}
	}
	double evaluate(int d, int a, int r){
		double sum = 0;
		if (depth == d){
				sum = log_gammaratio(a, count);
				//cerr << "count1: " << count << endl;
		}
		else if (depth <= d - 1){
			for (int v = 0; v < arity; v ++){
				if (children[v] != NULL){
					double value = children[v]->evaluate(d, a, r);
					sum += value;
				}
			}
			if (depth == d - 1){
				sum -= log_gammaratio(r, count);
				//cerr << "count2: " << count << endl;
			}
		}
		return sum;
	}
};


class Model{
public:
	Model(){}
	~Model(){}
	
	void init(){
		data.init();
		m = data.numrecords;
		n = data.numattributes;

		layering.init(data);

		for (int a = 1; a < 16; a ++){
			double aa = 1.0 / a;
			lr[a][0] = 0; lr2[a][0] = 0;
			for (int k = 1; k < MAX_COUNT; k ++){
				lr[a][k] = lr[a][k-1] + log(aa + k - 1);
				lr2[a][k] = lr2[a][k-1] + log(a + k - 1);
			}
		}
		
		//test();		
	}
	void test(){
		for (int i = 0; i < 1000000; i ++){
			vector<int> T;
			T.clear();
			T.push_back(1);T.push_back(2);
			T.push_back(3);T.push_back(4);
			cerr<<" Test ("<<i<<") : "<<log_lcp(5 + i % 8, T)<<endl;
		}
		exit(1);
	}

	//--- Interface functions --- below ----------
	//
	// Log of local conditional probability, log p(xi | xS).
	double log_lcp(int i, vector<int> & T){
		
		return log_lcp_mult(i, T);
	}
	double log_lcp(int i, int *T, int d){
		vector<int> TT;
		for (int j = 0; j < d; j ++) TT.push_back(T[j]);
		return log_lcp_mult(i, TT);
	}
	// Prio probability of T as the parents of i.
	double log_prior(int i, int *T, int d){
		// - Log {n-1 choose d}.
		return lr[1][d] + lr[1][n-1-d] - lr[1][n-1];
	}
	int num_layers(){
		return layering.numlayers;
	}
	void layer(int h, int** Vh, int* nh){
		*Vh = &(layering.V[layering.cnh[h]]);
		*nh = layering.nh[h];
	}
	void upper_layers(int h, int** Vu, int* nu){//similar to above but upper means parent
		if (h < layering.numlayers){
			*Vu = &(layering.V[layering.cnh[h+1]]);
			*nu = n - layering.cnh[h+1];
		}
		else{ *Vu = NULL; *nu = 0;}	
	}
	void lower_layers(int h, int **Vl, int *nl){//similar to above but lower means children
		if (h > 0){
			*Vl = layering.V;
			*nl = layering.cnh[h];
		}
		else{ *Vl = NULL; *nl = 0;}
	}
	int max_indegree(){
		return atoi(Arguments::maxindegree);
	}
	bool edges_within(int h){
		return layering.edges_within(h);
	}
	int num_nodes(){
		return n;
	}
	void print_edge_prob(ostream& f, int i, int j, double p){
		f<<data.heads[i]<<" -> "<<data.heads[j]<<" "<<p<<endl; 
	} 
	
	//
	//--- Interface functions --- above ----------	
	
	
private:	
	// Dirichlet-multinomial model
	double log_lcp_mult(int i, vector<int> & T){
		vector<int> S;
		S.push_back(i);
		for (int k = 0; k < (int)T.size(); k ++){
			S.push_back(T[k]);
		}
		
		//cerr << " Evaluation yields: " << evaluate(S) << " , log_gammaratio(1, 2) = " 
		//	<< log_gammaratio(1, 2) << " , a = " <<  data.arities[S[0]] << endl;
		
		return evaluate(S);	
	}
	double evaluate(vector<int> & S){
		Tnode root(data.maxarity+1, 0);
		vector<int> u;
		for (int t = 0; t < m; t ++){
			u.clear();
			for (int k = 0; k < (int)S.size(); k ++){
				int j = S[k];
				u.push_back(data.dm[t][j]);
			}
			//cerr << t+1 << "th update: ";
			root.update(S, u);
			//cerr << endl;
		}
		//root.print(cerr);
	
		return root.evaluate(S.size(), 1, -data.arities[S[0]]);
	}
		
	Data data;
	int n;
	int m;
	Layering layering;	
	
};





#endif
