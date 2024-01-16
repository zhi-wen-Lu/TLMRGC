
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <list>
#include <iomanip>
#include <string.h>
#include "sais.h"
#include "libbsc/bsc.h"
using namespace std;
int Ref_len,compress_len,tar_len,Ver_len,Dec_len;
//char *Ref_seq;
char *Dec_seq;
char *Ver_seq;
char *Comp_seq;
char *temp;
const int min_size = 1<<23;
int pos_vec_len;
static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
inline std::string generateString(const std::string &chr, int length = 5) {
	srand(time(0));
	string res = chr + "_";
	for (int i = 0; i < length; ++i) {
		res += alphanum[rand() % 62];
	}
	return res;
}
struct POSITION_RANGE{
	int begin, length;
};
POSITION_RANGE *pos_vec = new POSITION_RANGE[min_size];


void readTar(FILE *Tstream){
	vector<char> taxaBuffer;
	char c = getc(Tstream);
	int _ver_len = 0;
	while(c != EOF){
		if(c == '>'|| c == '@'){
			while ((c = getc(Tstream))!= '\n')
			{
				taxaBuffer.push_back(c);
			}
		}
		else{
			while ((c = getc(Tstream)) != '>' && c != EOF)
			{
				if(!isalpha(c)){
					continue;
				}
				Ver_seq[_ver_len++] = c;
				if(_ver_len ==13570){
					int ww =0;
				}
			}
			
		}
	}
	Ver_len =_ver_len;
}

void readDec(FILE *Tstream){
	vector<char> taxaBuffer;
	char c = getc(Tstream);
	int _Dec_len = 0;
	while(c != EOF){
		if(c == '>'|| c == '@'){
			while ((c = getc(Tstream))!= '\n')
			{
				taxaBuffer.push_back(c);
			}
		}
		else{
			while ((c = getc(Tstream)) != '>' && c != EOF)
			{
				if(!isalpha(c)){
					continue;
				}
				Dec_seq[_Dec_len++] = c;
				if(_Dec_len ==13570){
					int ww =0;
				}
			}
			
		}
	}
	Dec_len =_Dec_len;
}

void Verification(char *Ver_seq,char *Tar_seq){
	int tot = 0;
	if(Dec_len == Ver_len){
		for(int y = 0; y < Ver_len; y++){
			if(Dec_seq[y] == Ver_seq[y]){
				tot ++;
			}else if(Dec_seq[y] != Ver_seq[y]){
				int m = y;
				cout << m << endl;
			}
		}
		cout << "success" << endl;
	}else if(Dec_len != Ver_len){
		for(int y = 0; y < Ver_len; y++){
			if(Dec_seq[y] == Ver_seq[y]){
				tot++;
			}else if(Dec_seq[y] != Ver_seq[y]){
				int m = y;
				cout << "error" << endl;
				cout << m << endl;
			}
		}
	}
	int ww = tot;
}



int main(int argc, char **argv){
	int c;
	int rmq =0;
	bool flag = false;;
	while ((c = getopt (argc, argv, "r:")) != -1)
		switch (c){
		case 'r':
	        rmq = 1;
        	break;
		if (isprint (optopt))
                	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
             	else
               		fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
		return 1;
		default:
		abort ();
	}
	char defile[150];

   vector<string> filename = {"chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", //default chr name list
                    "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa", "chr10.fa", 
                    "chr11.fa", "chr12.fa", "chr13.fa", "chr14.fa", "chr15.fa", "chr16.fa", "chr17.fa", 
                    "chr18.fa", "chr19.fa", "chr20.fa", "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"};
   for(int d = 0 ; d < 24; d ++){
	   char ref[150];char tar[150];char code[150];char decode[150];
		FILE *Dstream;
		FILE *Tstream;
		sprintf(ref,"%s%s",argv[argc-1],filename[d].c_str());
		sprintf(decode,"%s%s",argv[argc-2],filename[d].c_str());
		Dstream = fopen(decode, "r");
		Tstream = fopen(tar, "r");
		int max = 1 << 28;
		Comp_seq = new char[max];
		temp= new char[max];		
		Dec_seq = new char[max];	
		Ver_seq = new char[max];
		readTar(Tstream);
		readDec(Dstream);
		Verification(Ver_seq,Dec_seq);
		

		
   }

	return EXIT_SUCCESS;
}
