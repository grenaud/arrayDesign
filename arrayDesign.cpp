/*
 * arrayDesign 
 * Date: Apr-17-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com (based on http://plindenbaum.blogspot.de/2011/08/memory-mapping-human-genome-with-mmap.html)
 *
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <stdint.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <zlib.h>

#include "utils.h"

using namespace std;



typedef struct{
    string chrName;
    unsigned int coordinate;
    char ref;
    char alt;
} SNVvariation;




/* it is the structure created by samtools faidx */
typedef struct faidx1_t
	 {
	 int32_t line_len, line_blen;
	 int64_t len;
	 uint64_t offset;
	 }faidx1_t,*FaidxPtr;
/**
 * wrapper for a mmap, a fileid and some faidx indexes
 */
class IndexedGenome
	{
	private:
	    /* used to get the size of the file */
	    struct stat buf;
	    /* genome fasta file file descriptor */
	    int fd;
	    /* the mmap (memory mapped) pointer */
	    char *mapptr;
	    /** reads an fill a string */
	    bool readline(gzFile in,string& line)
	    {
		if(gzeof(in)) return false;
		line.clear();
		int c=-1;
		while((c=gzgetc(in))!=EOF && c!='\n') line+=(char)c;
		return true;
	    }
	    
	public:
	    /* maps a chromosome to the samtools faidx index */
	    map<string,faidx1_t> name2index;

	    /** constructor 
	     * @param fasta: the path to the genomic fasta file indexed with samtools faidx
	     */
	    IndexedGenome(const char* fasta):fd(-1),mapptr(NULL)
	    {
		string faidx(fasta);
		//cout<<fasta<<endl;
			string line;
			faidx+=".fai";
			/* open *.fai file */
			//cout<<faidx<<endl;
			ifstream in(faidx.c_str(),ios::in);
			if(!in.is_open()){
			    cerr << "cannot open " << faidx << endl;
			    exit(EXIT_FAILURE);
			}
			/* read indexes in fai file */
			while(getline(in,line,'\n'))
	{
				if(line.empty()) continue;
				const char* p=line.c_str();
				char* tab=(char*)strchr(p,'\t');
				if(tab==NULL) continue;
				string chrom(p,tab-p);
				++tab;
				faidx1_t index;
				if(sscanf(tab,"%ld\t%ld\t%d\t%d",
					&index.len, &index.offset, &index.line_blen,&index.line_len
					)!=4)
					{
					cerr << "Cannot read index in "<< line << endl;
					exit(EXIT_FAILURE);
					}
				/* insert in the map(chrom,faidx) */
				name2index.insert(make_pair(chrom,index));
				}
			/* close index file */
			in.close();

			/* get the whole size of the fasta file */
			if(stat(fasta, &buf)!=0)
				{
				perror("Cannot stat");
				exit(EXIT_FAILURE);
				}
			
 			/* open the fasta file */
			fd = open(fasta, O_RDONLY);
   			 if (fd == -1)
				{
				perror("Error opening file for reading");
				exit(EXIT_FAILURE);
    				}
    			/* open a memory mapped file associated to this fasta file descriptor */
			mapptr = (char*)mmap(0, buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
    			if (mapptr == MAP_FAILED)
				{
				close(fd);
				perror("Error mmapping the file");
				exit(EXIT_FAILURE);
    				}
			}
		/* destructor */
		~IndexedGenome()
			{
	/* close memory mapped map */
			if(mapptr!=NULL && munmap(mapptr,buf.st_size) == -1)
				{
				perror("Error un-mmapping the file");
   			 	}
   			 /* dispose fasta file descriptor */
			if(fd!=-1) close(fd);
			}
			
		/* return the base at position 'index' for the chromosome indexed by faidx */
	    void print(const FaidxPtr faidx,int64_t index, int length,unsigned int sizeProbes,unsigned int tiling, const vector<SNVvariation> & vectorOfSNV,unsigned int maxreps){

		unsigned int indexVecVar=0;
		// cout<<"print "<<index<<" "<<length<<endl;
		for(long i=0;i<long(length);i+=long(tiling)){ //iterate over region
		    // cout<<"region "<<i<<endl;
		    int64_t index2=index;
		    //detecting variants
		    vector<SNVvariation>  snvToUse;
		    int64_t st=index2;
		    int64_t en=index2+sizeProbes;

		    for(unsigned int idxVecVar=indexVecVar;idxVecVar<vectorOfSNV.size();idxVecVar++){
			if(vectorOfSNV[idxVecVar].coordinate>en){//coordinate is still too far
			    break;
			}

			if(vectorOfSNV[idxVecVar].coordinate<st){//keep going, increase other start because we can discard those
			    indexVecVar++;
			}
		    
			if(vectorOfSNV[idxVecVar].coordinate >= st &&
			   vectorOfSNV[idxVecVar].coordinate <= en){//in range
			    snvToUse.push_back(vectorOfSNV[idxVecVar]);
			}
		    }

		    
		    
		    // cout<<"region "<<i<<endl;
		    if(snvToUse.empty()){ //no variants
			// cout<<"test "<<index2<<endl;
			string strToPrint="";
	    
			for(unsigned int j=0;j<sizeProbes;j++){ //for each char
			    long pos= faidx->offset +
				index2 / faidx->line_blen * faidx->line_len +
				index2 % faidx->line_blen
				;
			    //cout<<char(toupper(mapptr[pos]));
			    strToPrint+=char(toupper(mapptr[pos]));
			    index2++;
			}

			for(unsigned int rep=0;rep<maxreps;rep++){ //for each repetition
			    cout<<strToPrint<<endl;
			}

		    }else{ //there are variants in sequence we have to compute every 2^#variants possible sequences
			// cout<<"has variants"<<endl;
			vector<bool> hasVariant (sizeProbes,false);
			for(unsigned int varind=0;varind<snvToUse.size();varind++){
			    hasVariant[ snvToUse[varind].coordinate-index2 ]=true;//setting flag for having a variant at that position
			}
			// cout<<vectorToString(hasVariant)<<endl;
			// cout<<index2<<endl;
			unsigned int numberOfrepToUse=maxreps/pow(2, snvToUse.size() );
			
			if(pow(2, snvToUse.size()) > maxreps) {
			    cerr<<"major logic issue: the number of reps is greater than the max previously found "<<maxreps<<"\t"<<pow(2, snvToUse.size() )<<endl;
			    exit(1);
			}
			    
			if(numberOfrepToUse!=1 && //if equal to 1, that means that snvTouse is maximal already
			   (numberOfrepToUse%2) !=0){
			    cerr<<"major logic issue: the number of reps is not congruent to 0 mod 2 major "<<maxreps<<"\t"<<pow(2, snvToUse.size() )<<"\t"<<numberOfrepToUse<<"\t"<<(numberOfrepToUse%2)<<endl;
			    exit(1);
			}
			
			for(unsigned int varind=0;varind<pow(2, snvToUse.size() );varind++){//for each the 2^#var sequences
			    index2=index;
			    // cout<<endl<<"varind "<<varind<<endl;
			    vector<bool> touse;
			    unsigned int indexInsnvToUse=0;
			    for(unsigned int k=0;k<snvToUse.size();k++){//generating possible booleans	    
				unsigned int mask=pow(2,k);
				touse.push_back(bool((varind & mask)>>k));
			    }
			    //cout<<vectorToString(touse)<<endl;	
			    if( touse.size() != snvToUse.size() ){
				cerr<<"Internal error"<<endl;
				exit(1);
			    }
			    // cout<<"to use "<<vectorToString(touse)<<endl;
			    string strToPrint="";
			    for(unsigned int j=0;j<sizeProbes;j++){
				long pos= faidx->offset +
				    index2 / faidx->line_blen * faidx->line_len +
				    index2 % faidx->line_blen
				    ;
				if(hasVariant[j]){
				    if(snvToUse[indexInsnvToUse].ref != char(toupper(mapptr[pos]))){
					cerr<<"ERROR: Discrepancy between reference allele for variant at "<<snvToUse[indexInsnvToUse].coordinate<<" and "<<index2<<endl;
					exit(1);
				    }

				    if(touse[indexInsnvToUse]){ //use reference variant
					//cout<<char(toupper( snvToUse[indexInsnvToUse].alt ));
					strToPrint+=char(toupper( snvToUse[indexInsnvToUse].alt ));
				    }else{ //use alternative variant					
					//cout<<char(toupper( snvToUse[indexInsnvToUse].ref ));
					strToPrint+=char(toupper( snvToUse[indexInsnvToUse].ref ));
				    }

				    indexInsnvToUse++;
				}else{
				    //cout<<char(toupper(mapptr[pos]));
				    strToPrint+=char(toupper(mapptr[pos]));
				}
				index2++;
			    }
			    for(unsigned int rep=0;rep<numberOfrepToUse;rep++){ //for each repetition
				cout<<strToPrint<<endl;
			    }

			
			}//for varind
			//cout<<endl<<"END "<<endl;

		    }//end else if no variant
		    index+=tiling;
		}//end  for(long i=0;i<long(length);i+=long(tiling)){ //iterate over region
	    }//end print
};

int main (int argc, char *argv[]) {
    unsigned int sizeProbes=35;
    int64_t tiling=1;
    vector<SNVvariation> vectorOfSNV;
    string snvFilename="";

    const string usage = "  "+string(argv[0])+" (options) <path to faidx chromosome> <bed with non-overlapping sorted list of region>\n\n"+
	"This program list the sequences for a given list of regions from a bed file for a single chromosome\n"+
	"Options:\n"+
	"\t"+"--size"+"\t"+"Size of the probes (default: "+stringify(sizeProbes)+")\n"+
	"\t"+"--tile"+"\t"+"Size of the tiling (default: "+stringify(tiling)+")\n"+
	"\t"+"--snv [file]"+"\t"+"File containing the single nucleotide variation in the following format\n\n\n"+
	"\t\tchr[tab]pos[tab]reference allele[tab]alternative allele\n\n"+
	"\t\tThey must contain a single chromosome and be sorted on the coordinate\n"+
	"\n";

   
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cerr<<"Usage:"<<usage<<endl;
	return 1;
    }

    
    for(int i=1;i<(argc-2);i++){


        if(string(argv[i]) == "--size" ){
            sizeProbes =destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "--tile" ){
            tiling =destringify<int>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "--snv" ){
            snvFilename = string(argv[i+1]);
	    i++;
            continue;
        }

    }

    string line;
    bool  foundVariants=false;
    string chrFound="";
    unsigned int indexVecVar=0;
    unsigned int coordFound=0;
    unsigned int maxVarInProbe=0;

    if(snvFilename != "" ){
	ifstream snvFile;
	snvFile.open(snvFilename.c_str(), ios::in);

	foundVariants=true;
	if (snvFile.is_open()){

	    if( getline (snvFile,line)){
		//cout<<"test1 "<<line<<endl;
		vector<string> fields=allTokens(line,'\t');
		SNVvariation toadd;
		toadd.chrName    =                           fields[0];
		toadd.coordinate = destringify<unsigned int>(fields[1]);
		toadd.ref        =                           fields[2][0];
		toadd.alt        =                           fields[3][0];
		vectorOfSNV.push_back(toadd);
		// cout<<vectorOfSNV.size()<<endl ;
		chrFound   = toadd.chrName;
		coordFound = toadd.coordinate;
		if(toadd.ref == toadd.alt){
		    cerr<<"The following line has the same ref. as alt. allele "<<line<<endl;
		    return 1;
		}
	    }

	    while ( getline (snvFile,line)){
		//cout<<"test2 "<<line<<endl;
		if(line.empty())
		    continue;
		vector<string> fields=allTokens(line,'\t');
		SNVvariation toadd;
		toadd.chrName    =                           fields[0];
		toadd.coordinate = destringify<unsigned int>(fields[1]);
		toadd.ref        =                           fields[2][0];
		toadd.alt        =                           fields[3][0];
		
		if(chrFound != toadd.chrName){
		    cerr<<"Cannot modify chromosome halfway in variant file from "<<chrFound<<" to "<<toadd.chrName<<endl;
		    return 1;
		}

		if(coordFound > toadd.coordinate){
		    cerr<<"Coordinate must be sorted last found =  "<<coordFound<<" now found = "<<toadd.coordinate<<endl;
		    return 1;
		}

		if(toadd.ref == toadd.alt){
		    cerr<<"The following line has the same ref. as alt. allele "<<line<<endl;
		    return 1;
		}

		coordFound=toadd.coordinate;
		vectorOfSNV.push_back(toadd);
	    }
	}else{
	    cerr << "Unable to open file "<<snvFilename<<endl;
	    return 1;
	}
	snvFile.close();

	////////////////////////////////////////////
	//computing max variations per probe size
	////////////////////////////////
	cerr<<"First pass, reading file for the first time to determine max # of variants"<<endl;
	ifstream myFile;
	unsigned int stFound=0;
	unsigned int enFound=0;
	string filename = string(argv[argc-1]);
	myFile.open(filename.c_str(), ios::in);
	bool firstLine=true;
	if (myFile.is_open()){
	    while ( getline (myFile,line)){
		vector<string> fields=allTokens(line,'\t');
		string       ch=fields[0];
		unsigned int st=destringify<unsigned int>(fields[1]);
		unsigned int en=destringify<unsigned int>(fields[2]);

		if(en<st){
		    cerr<<"Cannot have start greater than end for line   "<<line<<endl;
		    return 1;
		}

		if(firstLine){
		    chrFound=ch;
		    firstLine=false;
		    stFound=st;
		    enFound=en;
		}else{
		    if(chrFound != ch){
			cerr<<"Cannot modify chromosome halfway in bed file from "<<chrFound<<" to "<<ch<<endl;
			return 1;
		    }
		
		    if(en<enFound){
			cerr<<"Cannot have end coordinate lesser than the previous one for "<<line<<endl;
			return 1;
		    }

		    if(st<stFound){
			cerr<<"Cannot have start coordinate lesser than the previous one for "<<line<<endl;
			return 1;
		    }

		}

		// if(genome->name2index.find(ch) == genome->name2index.end()){
		//     cerr<<"Cannot find chr "<<ch<<endl;
		// }else{
		// }
		// faidx1_t & findx=genome->name2index[ch];

		//finding variants in that region
		vector<SNVvariation> tosend;
		if(foundVariants){
		    for(unsigned int i=indexVecVar;i<vectorOfSNV.size();i++){
			if(vectorOfSNV[i].coordinate>en){//coordinate is still too far
			    break;
			}

			//plus 2 due to the negative and 
			if(vectorOfSNV[i].coordinate<(st+2)){//keep going, increase other start because we can discard those 
			    indexVecVar++;
			}
		    
			if(vectorOfSNV[i].coordinate >= (st+2) &&
			   vectorOfSNV[i].coordinate <= en){//in range
			    tosend.push_back(vectorOfSNV[i]);
			}
		    }
		}
		
		//checking region 
		//  void print(const FaidxPtr faidx,int64_t index, int length,unsigned int sizeProbes,unsigned int tiling, const vector<SNVvariation> & vectorOfSNV){
		unsigned int indexVecVar=0;
		// cout<<"print "<<index<<" "<<length<<endl;
		int length=en-st+1-sizeProbes;
		int64_t index=st;
		for(long i=0;i<long(length);i+=long(tiling)){ //iterate over region
		    int64_t index2=index;
				    
		    vector<SNVvariation>  snvToUse;
		    int64_t st=index2;
		    int64_t en=index2+sizeProbes;

		    for(unsigned int idxVecVar=indexVecVar;idxVecVar<vectorOfSNV.size();idxVecVar++){
			if(vectorOfSNV[idxVecVar].coordinate>en){ //coordinate is still too far
			    break;
			}

			if(vectorOfSNV[idxVecVar].coordinate<st){//keep going, increase other start because we can discard those
			    indexVecVar++;
			}
		    
			if(vectorOfSNV[idxVecVar].coordinate >= st &&
			   vectorOfSNV[idxVecVar].coordinate <= en){//in range
			    snvToUse.push_back(vectorOfSNV[idxVecVar]);
			}
		    }

		    if(snvToUse.size()>maxVarInProbe){
			maxVarInProbe=snvToUse.size();
		    }

		    index+=tiling;

		}//end for(long i=0;i<long(length);i+=long(tiling)){ //iterate over region
		// genome->print(&findx,st,en-st+1-sizeProbes,sizeProbes,tiling,tosend);
		//cout<<endl;
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<filename<<endl;
	    return 1;
	}


    }

    maxVarInProbe=pow(2,maxVarInProbe);//we need 2^(#variant reads)

    
    ifstream myFile;
    string genomeFile = string(argv[argc-2]);
    //cerr<<genomeFile<<endl;
    IndexedGenome* genome=new IndexedGenome(genomeFile.c_str());
    cerr<<"mapped into memory"<<endl;
    unsigned int stFound=0;
    unsigned int enFound=0;
	
    string filename = string(argv[argc-1]);
    myFile.open(filename.c_str(), ios::in);
    bool firstLine=true;
    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    vector<string> fields=allTokens(line,'\t');
	    //vector<string> fields2=allTokens(fields[1],'-');
	    //cerr<<"line "<<line<<endl;
	    string       ch=fields[0];
	    unsigned int st=destringify<unsigned int>(fields[1]);
	    unsigned int en=destringify<unsigned int>(fields[2]);

	    if(en<st){
		cerr<<"Cannot have start greater than end for line   "<<line<<endl;
		return 1;
	    }

	    // if(st!=0)
	    // 	st--;

	    if(firstLine){
		chrFound=ch;
		firstLine=false;
		stFound=st;
		enFound=en;
	    }else{
		if(chrFound != ch){
		    cerr<<"Cannot modify chromosome halfway in bed file from "<<chrFound<<" to "<<ch<<endl;
		    return 1;
		}
		
		if(en<enFound){
		    cerr<<"Cannot have end coordinate lesser than the previous one for "<<line<<endl;
		    return 1;
		}

		if(st<stFound){
		    cerr<<"Cannot have start coordinate lesser than the previous one for "<<line<<endl;
		    return 1;
		}

	    }

	    if(genome->name2index.find(ch) == genome->name2index.end()){
		cerr<<"Cannot find chr "<<ch<<endl;
	    }else{
		//cout<<"found"<<endl;
	    }
	    faidx1_t & findx=genome->name2index[ch];

	    //finding variants in that region
	    vector<SNVvariation> tosend;
	    if(foundVariants){
		// cout<<"index "<<indexVecVar<<"\t"<<vectorOfSNV.size()<<endl;
		for(unsigned int i=indexVecVar;i<vectorOfSNV.size();i++){
		    // cout<<"index "<<indexVecVar<<"\t"<<vectorOfSNV[i].coordinate<<endl;
		    if(vectorOfSNV[i].coordinate>en){//coordinate is still too far
			break;
		    }

		    //plus 2 due to the negative and 
		    if(vectorOfSNV[i].coordinate<(st+2)){//keep going, increase other start because we can discard those 
			indexVecVar++;
		    }
		    
		    if(vectorOfSNV[i].coordinate >= (st+2) &&
		       vectorOfSNV[i].coordinate <= en){//in range
			tosend.push_back(vectorOfSNV[i]);
		    }
		}
	    }


	    genome->print(&findx,st,en-st+1-sizeProbes,sizeProbes,tiling,tosend,maxVarInProbe);
	    //cout<<endl;
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	return 1;
    }

    cerr<<"finished succesfully"<<endl;
    return 0;
}

