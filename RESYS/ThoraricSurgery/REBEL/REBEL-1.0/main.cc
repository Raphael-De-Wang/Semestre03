#include<stdio.h>
#include<stdlib.h>
#include"math.h"

#include"Arguments.h"
#include"Model.h"
#include"Engine.h"

void welcome(){
	fprintf(stderr,
	" ~~~ h e l l o ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); 
	fprintf(stderr, 
	" REBEL       Rapid Exact Bayesian Edge Learning             v 0.1\n");
	fprintf(stderr,
	" (c) 2005    Mikko Koivisto      HIIT BRU, University of Helsinki\n"); 
	fprintf(stderr,
	" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); 

}
void goodbye(){
	fprintf(stderr,
	" ~~~ g o o d b y e ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}



int main(int argc, char **argv){
	welcome();

	Arguments::init(argc, argv);
		
	Model model;

	model.init();
	
	Engine engine;

	engine.init(&model);

	engine.compute_edge_probabilities();

	goodbye();
	return 1;
}
