#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include<stdio.h>
#include<stdlib.h>

using namespace std;

class Arguments {
public:
	static char* datafile;
	static char* layeringfile;
	static char* maxindegree;
	static char* model;
	static char* task;
	static char* maxnumrecords;
	

	static void init(int argc, char **args){
		for(int i = 1; i < argc; i ++){
			
			if(args[i][0]=='-' && args[i][1]=='d' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::datafile = args[j];
					j ++;
				} 
			}
			else if(args[i][0]=='-' && args[i][1]=='l' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::layeringfile = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='m' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::maxindegree = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='u' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::maxnumrecords = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='M' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::model = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='T' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::task = args[j]; 
					j ++;
				}
			}
		}
	
		print_arguments(stderr);
	}
	static void print_arguments(FILE *f){
		fprintf(f, " -d Data file:\n");
		fprintf(f, "   %62s\n", Arguments::datafile);	
		//fprintf(f, " -l Layering file:\n");
		//fprintf(f, "   %62s\n", Arguments::layeringfile);
		fprintf(f, " -m Maximum indegree:\n");
		fprintf(f, "   %62s\n", Arguments::maxindegree);
		fprintf(f, " -u Maximum number of data records read:\n");
		fprintf(f, "   %62s\n", Arguments::maxnumrecords);
		//fprintf(f, " -M Model:\n");
		//fprintf(f, "   %62s\n", Arguments::model);
		//fprintf(f, " -T Task (Infer=I, Generate=G):\n");
		//fprintf(f, "   %62s\n", Arguments::task);
	}
};

char* Arguments::datafile = "testdata.dat";
char* Arguments::layeringfile = "%";
char* Arguments::maxindegree = "3";
char* Arguments::model = "M";
char* Arguments::task = "I";
char* Arguments::maxnumrecords = "999999";


#endif
