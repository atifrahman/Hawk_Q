//
//  main.cpp
//  kmer
//

#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
using namespace std;

#include <pthread.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "kmer.h"

#include "specialfunctions.h"


#define MAX_REC_LEN 10240

int noSamples;

int kmerLength=31;

int NUM_THREADS=16;

pthread_mutex_t readFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t outFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t hashTable_mutex = PTHREAD_MUTEX_INITIALIZER;

#pragma pack(push, 1)
class Kmer
{
public:
	long long int kmer;
	unsigned short * counts;
    	double pVal;
    
    Kmer(int noSamples);
    ~Kmer();
    void show();
	void freeMemory();
};
#pragma pack(pop)

class KeyVal
{
public:
	long long int val;
	int count;
};

class HashTable
{
public:
    long long int totalKmers;
    long long int totalTests;
    vector <Kmer *> kmers[HASH_TABLE_LENGTH];
	pthread_mutex_t hashTableBuckets_mutex[HASH_TABLE_LENGTH];

	
	long long int * totalKmerCounts;
	double * phenotypeValues;
	double likelihoodNULL;

	
    HashTable();
    void insertKmer(long long int val, int count, int sampleNo);
    void show();
	void computeNullLikelihood();		
    void computePvalues();
    void computeLikelihoodRatios();	
    void dumpKmers(double sigLevel);
    
};

HashTable *ht;

class Factorials
{
public:
    double *factorials;
    Factorials(int max);
    double getFactorial(int number);
};

Factorials *f;

long long int getInt(char *s)
{
	long long int val=0;
	int i=0;
	char ch;
	while(s[i])
	{
		val=val<<2;
		
		ch=s[i];
		if(ch=='A')
		{
			val=val|0;
		}
		else if(ch=='C')
		{
			val=val|1;
		}
		else if(ch=='G')
		{
			val=val|2;
		}
		else
		{
			val=val|3;
		}
		i++;
        
	}
	return val;
}

char * getKmer(long long int val, char *kmer, int kmerLength)
{
    
	int temp=0;
	for(int i=kmerLength-1;i>=0;i--)
	{
		temp=val&3;
		val=val>>2;
		if(temp==0)
			kmer[i]='A';
		else if(temp==1)
			kmer[i]='C';
		else if(temp==2)
			kmer[i]='G';
		else if(temp==3)
			kmer[i]='T';
	}
	kmer[kmerLength]='\0';
	return kmer;
}


unsigned long int getHash(unsigned long long int key)
{
    /*
     key = (~key) + (key << 18); // key = (key << 18) - key - 1;
     key = key ^ (key >> 31);
     key = key * 21; // key = (key + (key << 2)) + (key << 4);
     key = key ^ (key >> 11);
     key = key + (key << 6);
     key = key ^ (key >> 22);
     */
	return (unsigned long int) key;
}


Kmer::Kmer (int noSamples)
{
	counts=new unsigned short[noSamples];
	for(int i=0;i<noSamples;i++)
        	counts[i]=0;
//	caseCounts=(char *)calloc(noCases,sizeof(char));
    
	
    pVal=0;

}

Kmer::~Kmer()
{
//    delete [] caseCounts;
//    delete [] controlCounts;
}

void Kmer::freeMemory()
{
	delete [] counts;
//	free(caseCounts);
//	free(controlCounts);
}


void Kmer::show()
{
    char kmerString[100];
    
    cout<<getKmer(kmer,kmerString,31)<<" ";
    
    for(int k=0;k<noSamples;k++)
    {
        cout<<(int)counts[k]<<" ";
        
    }
    cout<<pVal<<endl;
    
}


HashTable::HashTable()
{
    totalKmerCounts=new long long int[noSamples];
    phenotypeValues=new double[noSamples];	
    for(int i=0;i<noSamples;i++)
        totalKmerCounts[i]=0;
    
	for(int i=0;i<HASH_TABLE_LENGTH;i++)
	{
		pthread_mutex_init(&hashTableBuckets_mutex[i],NULL);
	}
    
    totalKmers=0;
    totalTests=0;
}


void HashTable::insertKmer(long long int val, int count, int sampleNo)
{
	unsigned long int index=getHash(val) % HASH_TABLE_LENGTH;
	int found=0;

	pthread_mutex_lock(&hashTableBuckets_mutex[index]);
    
	for(int i=0;i<kmers[index].size();i++)
	{
	 if(val==kmers[index][i]->kmer)
        {
            
              kmers[index][i]->counts[sampleNo]=count;
            	found=1;
		break;
	  }
	}
	if(found==0)
	{
		Kmer* km=new Kmer(noSamples);
		km->kmer=val;
		
              km->counts[sampleNo]=count;
        	kmers[index].push_back(km);
       
		totalKmers++;
        	

	}
	pthread_mutex_unlock(&hashTableBuckets_mutex[index]);
 
}

vector<double> regress(int noSamples, double *x, double *y)
{
	vector<double> result;
	double a,b;
	double x_mean=0, y_mean=0;
	double s_x=0,s_y=0,s_xy=0;

	for(int i=0;i<noSamples;i++)
	{
		x_mean+=x[i];
		y_mean+=y[i];
	}
	x_mean=x_mean/noSamples;
	y_mean=y_mean/noSamples;

	for(int i=0;i<noSamples;i++)
	{
		s_x+=(x[i]-x_mean)*(x[i]-x_mean);
		s_y+=(y[i]-y_mean)*(y[i]-y_mean);
		s_xy+=(x[i]-x_mean)*(y[i]-y_mean);
	}

	b=s_xy/s_x;
	a=y_mean-b*x_mean;
	result.push_back(a);
	result.push_back(b);
	return result;

}

double normalProb(double x, double m, double s)
{
	double logpi=0.9189385;
	  
  	return	-log(s)-logpi-(x-m)*(x-m)/(2*s*s);     
}

double logPI=0.9189385;

void HashTable::computeNullLikelihood()
{
	double *y=phenotypeValues;
	double y_mean=0;

	for(int k=0;k<noSamples;k++)
	{
		y_mean+=y[k];
	}
	y_mean=y_mean/noSamples;

	double e_null=0;

	for(int k=0;k<noSamples;k++)
	{
		e_null+=(y[k]-y_mean)*(y[k]-y_mean);
	}
	e_null=sqrt(e_null/noSamples);

	likelihoodNULL=0;

	for(int k=0;k<noSamples;k++)
	{
		likelihoodNULL+=(-log(e_null)-logPI-(y[k]-y_mean)*(y[k]-y_mean)/(2*e_null*e_null));
	}
}


void * likelihoodRatio_thread(void *threadid)
{
	long tid;
   	tid = (long)threadid;
    
	double likelihoodNull, likelihoodAlt, likelihoodRatio;

	double *y;

	double *x=new double[noSamples];

	double e_alt,y_p;

    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
            
		for(int k=0;k<noSamples;k++)
		{
			x[k]=ht->kmers[i][j]->counts[k]/(double)ht->totalKmerCounts[k];
		}          
    		y = ht->phenotypeValues;

		vector<double> result=regress(noSamples,x,y);
		
		e_alt=0;

		for(int k=0;k<noSamples;k++)
		{
			y_p=result[0]+result[1]*x[k];
			e_alt+=(y[k]-y_p)*(y[k]-y_p);
		}

		e_alt=sqrt(e_alt/noSamples);

		likelihoodNull=ht->likelihoodNULL, likelihoodAlt=0;

		for(int k=0;k<noSamples;k++)
		{
			likelihoodAlt+=(-log(e_alt)-logPI-(y[k]-(result[0]+result[1]*x[k]))*(y[k]-(result[0]+result[1]*x[k]))/(2*e_alt*e_alt));
		}

		likelihoodRatio=likelihoodAlt-likelihoodNull;
            
		if(likelihoodRatio<0)
		{
			likelihoodRatio=0;
		}		

	       double pVal=alglib::chisquarecdistribution(1, 2*likelihoodRatio);
	
     
            ht->kmers[i][j]->pVal=pVal;
	
 		

        }
    }

	delete []x;

	pthread_exit(NULL);


}

void HashTable::computeLikelihoodRatios()
{
    


	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
		  rc = pthread_create(&threads[t], NULL, likelihoodRatio_thread, (void *)t);
		  if (rc){
			 exit(-1);
		  }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
         		exit(-1);
      		}
      	}

}

FILE *outFile;
double pValThreshold;


void * dump_thread(void *threadid)
{
	long tid;
   	tid = (long)threadid;

	char kmerString[100];
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
            if(ht->kmers[i][j]->pVal<=pValThreshold)
            {
 			pthread_mutex_lock(&outFile_mutex);

                      fprintf(outFile,"%s\t%e\t",getKmer(ht->kmers[i][j]->kmer, kmerString, 31),ht->kmers[i][j]->pVal);

			for(int k=0;k<noSamples;k++)
            		{
                		fprintf(outFile,"%d\t",ht->kmers[i][j]->counts[k]);
                
            		}
			fprintf(outFile,"\n");
			pthread_mutex_unlock(&outFile_mutex);
                
            }
        }
        
    }
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
		ht->kmers[i][j]->freeMemory();
		delete ht->kmers[i][j];
        }
        ht->kmers[i].clear();
    }


	pthread_exit(NULL);

}


void HashTable::dumpKmers(double sigLevel)
{
    outFile=fopen("out_wo_bonf.kmerDiff","a");
    
	pValThreshold=sigLevel/(double)CUTOFF;    


	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
		  rc = pthread_create(&threads[t], NULL, dump_thread, (void *)t);
		  if (rc){
			 exit(-1);
		  }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
         		exit(-1);
      		}
      	}

	fclose(outFile);
    
}

void HashTable::show()
{
/*    
    for(int i=0;i<HASH_TABLE_LENGTH;i++)
    {
        for(int j=0;j<kmers[i].size();j++)
        {
            kmers[i][j]->show();
        }
        
    }
*/
    cout<<totalKmers<<endl;
    cout<<totalTests<<endl;
}


void getKeyVal(char *s, KeyVal* kv)
{
	long long int	val=0;
	int countVal=0;
	int i=0;
	char ch;
	while(i<KMER_LENGTH)
	{
		
		ch=s[i++];
		
		val=val<<2|bases[ch];		
        
	}
	i++;
	kv->val=val;
	while(1)
	{
		ch=s[i];
		if(ch=='\0'||ch=='\n')
		{	
			break;
		}
		countVal=countVal*10+ch-'0';
		i++;
	        
	}
	kv->count=countVal;

}

FILE ** kmerFiles;
long long int *vals;
unsigned short int * counts;

struct ThreadArg
{
	long long int valBar;
	int threadID;
};

void * readSamples(void *threadid)
{
	ThreadArg *ta=(ThreadArg *)threadid;
	long long int valBar=ta->valBar;
   	int threadNo=ta->threadID;
	long long int val;
	int count;

	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

	for(int i=threadNo;i<noSamples;i+=NUM_THREADS)
    		{
			if(vals[i]!=-1 && vals[i]<valBar)
			{
				ht->insertKmer(vals[i], counts[i], i);
				vals[i]=-1;
				counts[i]=-1;
			}
			if(vals[i]==-1)
			{
       		while(fscanf(kmerFiles[i],"%lld %d\n",&val,&count)!=EOF)
        		{
			
							
				if(val<valBar)
				{
 			           	ht->insertKmer(val, count, i);
            			}
				else
				{
					vals[i]=val;
					counts[i]=count;
					break;
				}
        		}
 			}
    		}


	delete []line;
	pthread_exit(NULL);

}




int main(int argc, const char * argv[])
{
    
    noSamples=atoi(argv[1]);
     
    ht=new HashTable();

	FILE *countsFile=fopen("sample_total_kmers.txt","r");
	FILE *phenotypeFile=fopen("phenotype_values.txt","r");

	 
	for(int i=0;i<noSamples;i++)
    	{
		fscanf(countsFile,"%lld\n",&ht->totalKmerCounts[i]);
	}
	for(int i=0;i<noSamples;i++)
    	{
		fscanf(phenotypeFile,"%lf\n",&ht->phenotypeValues[i]);
	}

	fclose(countsFile);
	fclose(phenotypeFile);


	ht->computeNullLikelihood();

	FILE *outFile=fopen("out_wo_bonf.kmerDiff","w");
    	
	fclose(outFile);


    
	kmerFiles=new FILE*[noSamples];
	char *kmerFilename;
    	kmerFilename=new char[5000];
    
	vals=new long long int[noSamples];
	counts=new unsigned short int[noSamples];
	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);



    
    char *temp;
    char kmerString[100];
    unsigned short count;
	
	long long int valMax=0x3FFFFFFFFFFFFFFF;
	long long int valBar=0;
	long long int val=0;
	long long int valInc=0x0010000000000000;
//	long long int valInc=0x0000000000100000;

	FILE *sortedFile=fopen("sorted_files.txt","r");
	 
	for(int i=0;i<noSamples;i++)
    	{
		fscanf(sortedFile,"%s\n",kmerFilename);
       //	sprintf(kmerFilename,"%d_case_kmers_sorted.txt",(i+1));
       	kmerFiles[i]=fopen(kmerFilename,"r");

		if(kmerFiles[i]==NULL)
		{
			cout<<kmerFilename<<" file doesn't exist"<<endl;
		}

		vals[i]=-1;
		counts[i]=-1;
		
	}
	ThreadArg *thArgs[NUM_THREADS];

	while(valBar<valMax)
	{
		valBar+=valInc;


		pthread_t threads[NUM_THREADS];
		int rc;
		void *status;
		for(int i=0;i<NUM_THREADS;i++)
		{
			thArgs[i]=new ThreadArg;
			thArgs[i]->valBar=valBar;
			thArgs[i]->threadID=i;
			rc = pthread_create(&threads[i], NULL, readSamples, (void *)thArgs[i]);
			if (rc)
			{
				 exit(-1);
			}
		}

		for(int i=0;i<NUM_THREADS;i++)
		{
			rc = pthread_join(threads[i], &status);
			delete thArgs[i];
			if (rc) 
	  		{
         			exit(-1);
      			}
		}

		
	    	ht->computeLikelihoodRatios();
    
    	
    		ht->dumpKmers(SIG_LEVEL);
 
		cout<<valBar<<endl;

	}

    	ht->show();




	outFile=fopen("out_wo_bonf.kmerDiff","r");

    FILE *fileOut=fopen("out_unique.kmerDiff","w");
  
	FILE *totalKmerFile=fopen("total_kmers.txt","w");
	fprintf(totalKmerFile,"%lld\n",ht->totalKmers);

	int count1, count2;
	double pVal;



    
    return 0;
}




