#include "Tools.h"

#if DEBUG
	#ifndef CK
		#define CK 1
	#endif

	#ifndef CK_2
		#define CK_2 0
	#endif
	
	#ifndef CK_3
		#define CK_3 1
	#endif
#endif

#ifndef STORE_ON_DISK
	#define STORE_ON_DISK 0
#endif

#define N_TAXA 128

//Input parameters
dataTypelenSeq k_min=16;

dataTypeNChar n=0; //eBWT size
string output_new; //output name
vector<string> NAME;

              //A B C D E F G H I J K L M N O P Q R S T U V W X Y
int ORD[25] = {0,0,1,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,3,0,0,0,0,0};
#define ord(c) (ORD[c-65])
//terminator character at the end of the reads
char END_STRING = '#';

dataTypeNSeq numData=0;

double cosine_dist_loop_array(bool *arr_a, bool *arr_b, vector<dataTypelenSeq> &sizes)
{
    double a_dot_b = 0.0;
    double a_mag = 0;
    double b_mag = 0;
	double dist = 1.0;
	
    for (size_t i = 0; i < sizes.size(); ++i)
    {
        a_dot_b += (arr_a[i] and arr_b[i])?(double)sizes[i]:0.0;
        a_mag += (arr_a[i])?(double)sizes[i]:0.0;
        b_mag += (arr_b[i])?(double)sizes[i]:0.0;
    }
	
	#if CK
	fprintf(stderr,"\na_dot_b %f, a_mag %f, b_mag %f\n", a_dot_b, a_mag, b_mag);
	#endif
	
	a_mag *=b_mag;
	if(a_mag > 0.0 && a_dot_b > 0.0)
		dist -= (a_dot_b/sqrt(a_mag));
	
	#if CK
	fprintf(stderr,"dist %f\n", dist);
	#endif
	
    return dist;
}

bool Close(dataTypeNChar ind_start, dataTypeNChar ind_fin, FILE *inCDA,FILE *inEBWT,std::bitset<N_TAXA> &isIn, dataTypeNChar *n_clust){
	
	//(*n_clust)++;
	bool res = false;
	//ind_start: starting index of the eBWT positional cluster
	//ind_fin: final index (NOT included)
	//isIn stores colors in the cluster
	//B_Symb stores symbols in the cluster
	std::bitset<ALF> B_Symb;
	
	//To upload in memory positional cluster in eBWT and cda
	
	dataTypeNChar c_len=ind_fin-ind_start;
	dataTypeNSeq *bufferCDA = new dataTypeNSeq[c_len];
	dataTypedimAlpha *buffereBWT = new dataTypedimAlpha[c_len];
	
	fseek(inCDA,ind_start*sizeof(dataTypeNSeq), SEEK_SET);
	dataTypeNChar numchar=fread(bufferCDA,sizeof(dataTypeNSeq),c_len,inCDA);
	assert(numchar==c_len);
	fseek(inEBWT,ind_start*sizeof(dataTypedimAlpha), SEEK_SET);
	numchar=fread(buffereBWT,sizeof(dataTypedimAlpha),c_len,inEBWT);
	assert(numchar==c_len);
	
	for(dataTypeNChar index=0; index<c_len; index++){
		isIn[bufferCDA[index]]=1;
		if(buffereBWT[index]==END_STRING){
			B_Symb[ALF-1]=1;
		}
		else{
			B_Symb[ord(buffereBWT[index])]=1;
		}
	}
	
	//Check positional cluster maximal
	if(B_Symb.count()>1){
		#if CK_3 
			fprintf(stderr,"Close returns 1 --> [%d,%d)\n",ind_start,ind_fin);
			fprintf(stderr,"CDA:\t");				
			for(dataTypeNChar index=0; index<c_len; index++)
				fprintf(stderr,"%d", (int)bufferCDA[index]);
			fprintf(stderr,"\neBWT:\t");
			for(dataTypeNChar index=0; index<c_len; index++)
				fprintf(stderr,"%c",buffereBWT[index]);
			fprintf(stderr,"\n");
		#endif
		res=true;
		(*n_clust)++;
	}
	#if CK_3 
	else{
		fprintf(stderr,"Close returns 0 --> [%d,%d)\n",ind_start,ind_fin);
		fprintf(stderr,"CDA:\t");						
		for(dataTypeNChar index=0; index<c_len; index++)
			fprintf(stderr,"%d", (int)bufferCDA[index]);
		fprintf(stderr,"\neBWT:\t");
		for(dataTypeNChar index=0; index<c_len; index++)
			fprintf(stderr,"%c",buffereBWT[index]);
		fprintf(stderr,"\n");
	}
	#endif
		
	delete [] bufferCDA;
	delete [] buffereBWT;
	
	return res;
}


#if STORE_ON_DISK
void WriteToFile(std::bitset<N_TAXA> &isIn, FILE* fileOutM){

	#if CK
		fprintf(stderr,"\nWriteToFile --> isIn: ");				
		for(dataTypeNSeq s=0; s<numData; s++)
			fprintf(stderr,"%d,", (int)isIn[s]);
		fprintf(stderr,"\n\n");
	#endif
	
	dataTypeNSeq *bufferM = new dataTypeNSeq[numData];
	for(dataTypeNSeq i=0; i<numData; i++)
		bufferM[i] = isIn[i]+48;
	
	fwrite(bufferM, sizeof(dataTypeNSeq),numData,fileOutM);
	
	delete [] bufferM;
	
	return;
	
}
#endif

#if STORE_ON_DISK
dataTypeNChar ComputeCosineSimMatrix(dataTypeNChar n_clusters, FILE *fileBinVect, FILE *fileOutMatrix, vector<dataTypelenSeq> &sizes){
#else
dataTypeNChar ComputeCosineSimMatrix(dataTypeNChar n_clusters, FILE *fileOutMatrix, vector<dataTypelenSeq> &sizes, vector<std::bitset<N_TAXA>> &summaryBinVect){
#endif
	
	//Allocate matrix to fill reading binary vectors
	bool** mat = new bool*[numData];
	for(dataTypeNSeq i=0; i<numData; i++)
		mat[i] = new bool[n_clusters];
	
	dataTypeNSeq row, row2;
	
	#if STORE_ON_DISK
	//Read file fileBinVect
	dataTypeNSeq *bufferM = new dataTypeNSeq[numData];
	dataTypeNSeq numRead = fread(bufferM, sizeof(dataTypeNSeq),numData,fileBinVect);
	//counter for lines
	dataTypeNChar j=0;
	while(numRead>0){
		for(row=0; row<numData; row++){
			mat[row][j] = (bufferM[row]=='0'?false:true);
		}
		j++;
		numRead = fread(bufferM, sizeof(dataTypeNSeq),numData,fileBinVect);
	}
	assert(j==sizes.size());
	delete [] bufferM;
	#else
	//Transpose summaryBinVect
	for(dataTypeNChar j=0; j<summaryBinVect.size(); j++){
		for(row=0; row<numData; row++){
			mat[row][j] = (summaryBinVect[j]).test(row);
		}
	}
	//Clear summaryBinVect
	vector<std::bitset<N_TAXA>>().swap(summaryBinVect);
	#endif
	
	fprintf(fileOutMatrix,"\t%d\n",numData);
	//Compute OutMatrix values
	for(row=0; row<numData; row++){ 
		fprintf(fileOutMatrix,"%s\t",&NAME[row][0]);
		for(row2=0; row2<row; row2++){
			double cosSim = cosine_dist_loop_array(mat[row],mat[row2],sizes);
			fprintf(fileOutMatrix,"%f\t",cosSim);
		}
		fprintf(fileOutMatrix,"%f\t",0.0);
		for(row2=row+1; row2<numData; row2++){
			double cosSim = cosine_dist_loop_array(mat[row],mat[row2],sizes);
			fprintf(fileOutMatrix,"%f\t",cosSim);
		}
		fprintf(fileOutMatrix,"\n");
	}
	
	for(dataTypeNSeq i=0; i<numData; i++)
		delete[] mat[i];
	delete[] mat;
	
	return 1;
}
#if STORE_ON_DISK
dataTypeNChar ClusterAnalysis(FILE* fileOut,FILE* fileInLCP,FILE *fileInCDA,FILE *fileInEBWT, dataTypeNChar *n_clust, vector<dataTypelenSeq> &sizes){
#else
dataTypeNChar ClusterAnalysis(FILE* fileInLCP,FILE *fileInCDA,FILE *fileInEBWT, dataTypeNChar *n_clust, vector<dataTypelenSeq> &sizes, vector<std::bitset<N_TAXA>> &summ){
#endif
	dataTypeNChar ClustersLength=0;
	
	//colors is a vector containing the IDs of CDA elements to analyze
	
	bool started=false;
	dataTypelenSeq pred = -1;
	bool grow = true;
	dataTypeNChar tail_len = 0;
	dataTypeNChar pos_init = 0, pos_end = 0;
	dataTypelenSeq c_weight=1;

	//Scan LCP file to detect maximal eBWT positional clusters
	//To read LCP and CDA
	dataTypeNChar numcharLCP;
	dataTypeNChar index_file=0;
	dataTypelenSeq *bufferLCP = new dataTypelenSeq[BUFFERSIZE];
	
	fseek(fileInLCP, 0, SEEK_SET);
	
	//Read BUFFERSIZE elements from LCP-array
	numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERSIZE,fileInLCP);
	
	while(numcharLCP>0){
		for(dataTypeNChar indexbuffer=0; indexbuffer<numcharLCP; indexbuffer++){	
			if(bufferLCP[indexbuffer]>=k_min) {
				if (not started){
					started=true, grow=true;
					pos_init=index_file+indexbuffer-1;
					c_weight=bufferLCP[indexbuffer];
					tail_len=0;
				}
				else //cluster is open
				{
					//Case1: increasing sequence of LCP values
					if (grow){
						if(bufferLCP[indexbuffer]<pred) grow = false;
					}
					//Case2:non-increasing sequence of LCP values
					else{ 
						if (bufferLCP[indexbuffer]<pred) tail_len=0; 
						else if (bufferLCP[indexbuffer]==pred) tail_len++;
						else { //close cluster and start a new one
							std::bitset<N_TAXA> isIn;
							pos_end = index_file + indexbuffer-1-tail_len;
							if(Close(pos_init,pos_end,fileInCDA,fileInEBWT,isIn,n_clust)==1){	
								//Close returns 1 if the eBWT positional cluster is maximal		
								c_weight=(c_weight>bufferLCP[indexbuffer-tail_len-2])?bufferLCP[indexbuffer-tail_len-2]:c_weight; //set c_weight								
								#if CK_3
								fprintf(stderr,"c_weight=%d\n",c_weight);
								#endif
								//Update votes
								//isIn: all colors in the cluster (not divided)
								ClustersLength+=c_weight;
								sizes.push_back(c_weight);
								#if STORE_ON_DISK
								WriteToFile(isIn,fileOut); 
								#else
								summ.push_back(isIn);
								#endif
							}
							//Start a new cluster
							started=true, grow=true;
							pos_init=pos_end;
							c_weight=bufferLCP[indexbuffer-tail_len];
							tail_len=0;
						}
					}
				}//end-else if(not started)		
			}
			else if (started){ //close cluster
				std::bitset<N_TAXA> isIn;
				pos_end = index_file + indexbuffer;
				if(Close(pos_init,pos_end,fileInCDA,fileInEBWT,isIn,n_clust)==1){
					c_weight=(c_weight>bufferLCP[indexbuffer-1])?bufferLCP[indexbuffer-1]:c_weight; //set c_weight
					#if CK_3
					fprintf(stderr,"c_weight=%d\n",c_weight);
					#endif
					ClustersLength+=c_weight;
					sizes.push_back(c_weight);
					#if STORE_ON_DISK
					WriteToFile(isIn,fileOut); 
					#else
					summ.push_back(isIn);
					#endif
				}
				started=false;
			}//end-else if(b_lcp[indexbuffer]>=k_min)
			pred=bufferLCP[indexbuffer];
		}//end-for
		
		index_file+=numcharLCP;
		//Read BUFFERSIZE elements from LCP-array and CDA-array
		numcharLCP=fread(bufferLCP,sizeof(dataTypelenSeq),BUFFERSIZE,fileInLCP);
		sizes.shrink_to_fit();
	}//end-while
			
	//Close a possibly last cluster 
	if(started){
		std::bitset<N_TAXA> isIn;
		if(Close(pos_init,index_file,fileInCDA,fileInEBWT,isIn,n_clust)==1){
			//Read last LCP value
			fseek(fileInLCP, (n-1)*sizeof(dataTypelenSeq), SEEK_SET);
			fread(bufferLCP,sizeof(dataTypelenSeq),1,fileInLCP);
			c_weight=(c_weight>bufferLCP[0])?bufferLCP[0]:c_weight; //set c_weight
			#if CK_3
			fprintf(stderr,"c_weight=%d\n",c_weight);
			#endif
			ClustersLength+=c_weight;
			sizes.push_back(c_weight);
			#if STORE_ON_DISK
			WriteToFile(isIn,fileOut); 
			#else
			summ.push_back(isIn);
			#endif
		}
	}	
	assert(index_file==n);
	
	#if CK && STORE_ON_DISK==0
	assert(*n_clust==summ.size());
	fprintf(stderr,"\nsumm size = %d\nsumm -->\n",*n_clust);			
	for(dataTypeNChar j=0; j<*n_clust; j++){
		for(dataTypeNSeq id=0; id<numData; id++){
			fprintf(stderr,"%d,", (int)summ[j][id]);
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	#endif
	
	delete [] bufferLCP;
	
	return ClustersLength;
}

#if STORE_ON_DISK
dataTypeNChar ClusterVotes(FILE* InFileLCP,FILE *InFileCDA, FILE *OutFile,FILE *InFileEBWT, FILE *f_out)
#else
dataTypeNChar ClusterVotes(FILE* InFileLCP,FILE *InFileCDA, FILE *OutFile,FILE *InFileEBWT)
#endif
{
	time_t t_total=0;
	clock_t c_total=0;

	time_start(&t_total, &c_total); //start time
	
	//Detect positional clusters and write Table partial_sum
	dataTypeNChar n_clust=0;
	
	vector<dataTypelenSeq> ClusterSize;
	
	#if STORE_ON_DISK
	dataTypeNChar len_clust_vote=ClusterAnalysis(f_out,InFileLCP,InFileCDA,InFileEBWT,&n_clust, ClusterSize);
	#else
	vector<std::bitset<N_TAXA>> summary;
	dataTypeNChar len_clust_vote=ClusterAnalysis(InFileLCP,InFileCDA,InFileEBWT,&n_clust, ClusterSize,summary);
	#endif
	
	fprintf(stderr,"ClusterVotes -> Total length of positional clusters %lu\n",len_clust_vote);
	fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	
	#if STORE_ON_DISK
	fseek(f_out, 0, SEEK_SET);
	assert(ComputeCosineSimMatrix(n_clust,f_out,OutFile,ClusterSize));
	#else
	assert(ComputeCosineSimMatrix(n_clust,OutFile,ClusterSize,summary));
	#endif
	
	return n_clust;
}



int main(int argc, char **argv) {
	
	time_t t_total=0;
   
    clock_t c_total=0;
	
	//INPUT
	if( argc < 5) {
		fprintf(stderr,"Error usage %s fileFasta fileInfo fileOutput k_min\n",argv[0]); 
		exit(1);
	}

	string fileFasta=argv[1];
	string fileInfo=argv[2];
	output_new=argv[3];

	sscanf(argv[4], "%d", &k_min);
	
	//print input parameters
	fprintf(stdout,"k_min: %d\n",k_min);
	
	fprintf(stderr,"sizeof(dataTypeNChar)=%lu bytes.\n",sizeof(dataTypeNChar));
	fprintf(stderr,"sizeof(dataTypelenSeq)=%lu bytes.\n",sizeof(dataTypelenSeq));
	fprintf(stderr,"sizeof(dataTypeNSeq)=%lu bytes.\n",sizeof(dataTypeNSeq));
	
	//Read fileInfo
	std::ifstream in_list;
	in_list.open(fileInfo.c_str(), std::ifstream::in);
	if (!in_list.is_open()){
		fprintf(stderr,"Error opening file %s\n",fileInfo.c_str());
		exit (EXIT_FAILURE);
	}
	
	while(getline(in_list,fileInfo,'\t')){
		NAME.push_back(fileInfo);
		getline(in_list,fileInfo,'\n');
		numData++;
	}
	
	fprintf(stdout,"numData: %u\n",numData);
	
	if(numData > N_TAXA){
		fprintf(stderr,"ERROR numData is %u and N_TAXA is %u\nPlease increase N_TAXA and compile again",numData,N_TAXA);
		exit(EXIT_FAILURE);
	}
	
    //Files .lcp, .ebwt and .cda 
	string fnLCP=fileFasta+".lcp";
	string fnCDA=fileFasta+".cda";
	string fnEBWT=fileFasta+".ebwt";
	
	FILE *InLCP= fopen(fnLCP.c_str(), "rb");
	FILE *InCDA= fopen(fnCDA.c_str(), "rb");
	FILE *InEBWT= fopen(fnEBWT.c_str(), "rb");
	
    if(InLCP==NULL){
		fprintf(stderr,"Error opening file %s\n",fnLCP.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if(InCDA==NULL){
        fprintf(stderr,"Error opening file %s\n",fnCDA.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
    }
    if (InEBWT==NULL){
        fprintf(stderr,"Error opening file %s\n",fnEBWT.c_str());
        printf("fopen failed, errno = %d\n", errno);
        exit (EXIT_FAILURE);
	}
	
	//Set the eBWT length
    fseek(InEBWT, 0, SEEK_END);
	n = ftell(InEBWT)/sizeof(dataTypedimAlpha);
	
	//Output file(s)
	//filebuf out_fb;
	//out_fb.open (output_new.c_str(),std::ios::out);
	//ostream OutNew(&out_fb);
	FILE *OutNew= fopen(output_new.c_str(), "w");
	
    #if STORE_ON_DISK
	std::stringstream ssout;
	ssout << "BinaryVectors_k_" << (int)k_min << "_" << output_new << ".txt\0";
	
	string fnOUT=ssout.str();
	FILE *OutP= fopen(fnOUT.c_str(), "w+");
    if ((OutP==NULL)){
		std::cerr << "Error opening " << fnOUT << "." << std::endl;
		printf("fopen failed, errno = %d\n", errno);
		exit (EXIT_FAILURE);
	}
    #endif
    
	time_start(&t_total, &c_total); //start time

	//Call compute binary matrix
	#if STORE_ON_DISK
		dataTypeNChar n_col = ClusterVotes(InLCP,InCDA, OutNew,InEBWT,OutP);
	#else
		dataTypeNChar n_col = ClusterVotes(InLCP,InCDA, OutNew,InEBWT);
	#endif
	
	#if STORE_ON_DISK
	fprintf(stdout,"\nOutput a binary matrix of size n_clusters x n_taxa = %ld x %u \n", n_col, numData);
	#endif
	fprintf(stdout,"\nOutput an upper triangular matrix of size n_taxa x n_taxa = %u x %u \n", numData, numData);

	//Close files
	fclose(InLCP);
    fclose(InCDA);
    fclose(InEBWT);

	#if STORE_ON_DISK
	fclose(OutP);
	#endif
	fclose(OutNew);
	
    fprintf(stdout,"Time: %.6lf\n", time_stop(t_total, c_total));
	return 0;
}

