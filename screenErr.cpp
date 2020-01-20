//test file: /home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w4000/h1.reads.to_cons.bam
// /home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w2000e02/h1.reads.to_cons.bam

#include "htslib/hts.h"
#include "htslib/sam.h"
#include <iostream>
#include <string>
#include <assert.h>
#include <algorithm>
#include <map>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;

#include <sstream>
//to measure time
#include <chrono>

class AlignedRead {
public:
	int index;
	string chrom;
	//alignment starting position on read
	int alignStart, alignEnd;
	//read starting position on reference
	int startPos;
	int endPos;
	int length;
	char* seq, *seqrc;
	string name;
	int n;
	bam1_t *aln;	

	~AlignedRead() {
		if (seq != NULL) {
			delete seq;
		}
	}
	int LookupRefPositionInRead(int refPos) {
		int curReadPos = alignStart;
		int curRefPos = startPos;
    	int k, l;
		uint32_t *cigar = bam_get_cigar(aln);
		if (refPos < startPos) {
			return -1;
		}
		if (refPos > endPos ) {
			return -1;
		}
    	for (k = 0; k < aln->core.n_cigar; ++k) {
			int opLen = bam_cigar_oplen(cigar[k]);
			int op = bam_cigar_op(cigar[k]);
			if (op == BAM_CSOFT_CLIP) {
				curReadPos = opLen;
				continue;
			}
			if ((op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) and 
					curRefPos <= refPos and curRefPos + opLen > refPos) {
				return curReadPos + (refPos-curRefPos);
			}
			else if (op == BAM_CDEL and curRefPos <= refPos and curRefPos + opLen > refPos) {
				return curReadPos;
			}
			if (bam_cigar_type(cigar[k]) & 1 ) {
				curReadPos += opLen;
			}
			if (bam_cigar_type(cigar[k]) & 2 ){
				curRefPos += opLen;
			}
		}
		return -1;
	}
};

class Valid_indel {
public:
	//size of the valid indel
	int indel_size;
	//the portion of reads span the indel that have the same indel
	float portion;
	//position of the valid indel on reference
	uint32_t indel_pos;
};

//first parameter: a bam1_t type pter variable
//second parameter: a reference to an AlignedRead variable, pass by reference
//function call: IndexBamRead(alignedRead, *curRead);
void IndexBamRead(bam1_t *alignedRead, AlignedRead &read) {
	//
	// Convert from bam to read.
	//
		 read.startPos=alignedRead->core.pos;
		 read.endPos=read.startPos;

		 
		 read.length = alignedRead->core.l_qseq;
		 read.seq=new char[read.length];
		 read.name = string(bam_get_qname(alignedRead));


		 uint8_t *q = bam_get_seq(alignedRead);
		 for (int i=0; i < read.length; i++) {
			 read.seq[i]=seq_nt16_str[bam_seqi(q,i)];	
		 }
		 uint32_t *cigarPtr = bam_get_cigar(alignedRead);

		 int leftClip = 0;
		 read.alignStart = 0;
		 if (alignedRead->core.n_cigar > 0 and bam_cigar_op(cigarPtr[0]) == BAM_CSOFT_CLIP) {
			 read.alignStart = bam_cigar_oplen(cigarPtr[0]);
		 }
		 else if (alignedRead->core.n_cigar > 1 and bam_cigar_op(cigarPtr[2]) == BAM_CSOFT_CLIP) {
			 read.alignStart = bam_cigar_oplen(cigarPtr[1]);
		 }

		 read.alignEnd = read.alignStart;
		 for (int i = 0; i < alignedRead->core.n_cigar; i++) {
			 /* bam_cigar_type returns a bit flag with:
			  *   bit 1 set if the cigar operation consumes the query
			  *   bit 2 set if the cigar operation consumes the reference
			 */
			//if consume reference
			 if (bam_cigar_type(cigarPtr[i]) & 0x2 ) {
				 read.endPos += bam_cigar_oplen(cigarPtr[i]);
			 }
			//if consume query
			 if (bam_cigar_type(cigarPtr[i]) & 1 ) {
				 read.alignEnd += bam_cigar_oplen(cigarPtr[i]);
			 }
		 }
		 read.aln=bam_dup1(alignedRead);
}

//check if current discovered indel is a new one
//if it is new, return ture
bool Check_new_indel(vector <Valid_indel> all_indel_info,
					 int32_t indel_pos,
					 int range_btw_indels) {
	bool new_indel = true;
	//loop through all_indel_info, to check if there is the same indel before
	for (vector <Valid_indel>::iterator it = all_indel_info.begin() ; 
		 it != all_indel_info.end();
		 it++) {
		int32_t pos = it -> indel_pos;
		if (indel_pos > pos - range_btw_indels &&
			indel_pos < pos + range_btw_indels) {
			new_indel = false;
			break;
		}
	}
    return new_indel;
}

int main(int argc, char* argv[]) {
	htsFile *htsfp;
	bam_hdr_t *samHeader;	
	
	if (argc <= 1) {
		exit(0);
	}
		
	htsfp = hts_open(argv[1],"r");
    samHeader = sam_hdr_read(htsfp);
	bam1_t *alignedRead = bam_init1();
	int bottom=0;
	int top=0;

    AlignedRead* curRead= new AlignedRead;
	int res=sam_read1(htsfp, samHeader, alignedRead);
	while ( res >= 0 and alignedRead->core.flag & 2304 ) {
		res=sam_read1(htsfp, samHeader, alignedRead);
	}

    if (res < 0) {
		cerr << "NO reads" << endl;
		exit(0);
	}
	
    IndexBamRead(alignedRead, *curRead);
    /*testing
    cout << curRead->startPos << endl;
    cout << curRead->endPos << endl;
    cout << curRead->alignStart << endl;
    cout << curRead->alignEnd << endl;
    */

    curRead -> index = 0;
    //
	// Read queue is a first in, first out data structure. 
	// For now the map implementation allows arbitrary indexing.
	// The key value is the order in which the read is stored in the alignments.
	//
	map<int, AlignedRead*> readQueue;
	//
	// Bottom is the read that has been in the queue the longest time,
    //e.g. smallest start position.
	//
    readQueue[bottom] = curRead;

    //This parameter define the region to check large indels
    //region_len times indel_size = the number of bases should be looked ahead/behind
    int region_len = 5;
	//lower bound of large indels
	int indel_min_size = 100;
	//position of the valid indel on reference
	int32_t indel_pos = 0;
	//size of the valid indel
	int indel_size = 0;
	//if the valid indel is the first in current bam file
	bool if_first_indel = true;
	//the position of previous valid indel on reference
	int32_t pre_indel_pos = 0;
	//if current indel is within the range of range_btw_indels,
	//then it is not counted as a new one
	int range_btw_indels = 1000;
	//store all the valid indel information: length, positions, portion
	vector<Valid_indel> all_indel_info;

	//search large indels
	//cigar string of alignedRead: a pter
	//curRead -> aln: duplication of original bam1_t read of curRead
	
	string dir = argv[1];
	//want dir to be the directory to store the text file below
	//e.g.: string dir = "/home/cmb-16/mjc/quentin/testksw/racon/build/bin/test/w300/";
	//20 is the length of h1.reads.to_cons.bam;
	dir = dir.substr(0, dir.length() - 20);
	//cout << dir << endl;

	//write screening results in a text document
	ofstream text;
	bool write_error_region = false;
	if (write_error_region) {
		text.open(dir + "racon_para_screenErr.txt");
		text << dir << "\n";
	}

	//loop while there is a new read
	int cur = top;
	while (curRead != NULL) {
		curRead = readQueue[cur];

		//Read in all sequences that may overlap with this read
		//Note that now overlaps should include both directions

		//forward direction: adding more overlapping reads
		//startPos: read starting position on reference
		if (readQueue[top] -> startPos < curRead -> endPos) {
			//while there is still a read
			while ( sam_read1(htsfp, samHeader, alignedRead) > 0) {
				//TODO: check this line of code
				while (alignedRead -> core.flag & 2304 && sam_read1(htsfp, samHeader, alignedRead) > 0) {}
				AlignedRead* nextRead = new AlignedRead;
				top++;
				nextRead -> index = top;
				readQueue[top] = nextRead;
				IndexBamRead(alignedRead, *nextRead);
				//cout << "storing " << nextRead->name << "\t" << top << "\t" << nextRead->startPos << endl;
				//if (nextRead->name == "") {
				//	cerr << "here!" << endl;
				//}
				if (alignedRead -> core.pos > curRead -> endPos) {
					break;
				}
			}
		}

		//delete bottom reads that not overlapping with curRead
		//note that the non-overlapping reads may not be at the bottom: lengths are different
		//thus, get rid of the non-overlapping reads at the bottom first
		//then check overlapping (with indels) every time before indels detection
		for (int j = bottom; j <= top; j++) {
			if (readQueue[j] -> endPos < curRead -> startPos) {
				readQueue.erase(bottom);
				bottom++;
			}
			else {
				break;
			}
		}

		//test if readQueue contains all the overlapping reads with curRead
		/*
		cerr << "bottom: " << bottom << " cur: " << cur << " top: " << top << endl;
		cerr << "bottom end: " << readQueue[bottom] -> endPos 
			 << " cur start: " << readQueue[cur] -> startPos
			 << " cur end: " << readQueue[cur] -> endPos
			 << " top start: " << readQueue[top] -> startPos << endl;
		*/

		//loop through all reads overlapping with curRead
		//to detect common large indels

		//curRead -> aln: duplication of original bam1_t read of curRead
		uint32_t* curCigarPtr = bam_get_cigar(curRead -> aln);
		//aligned bases of current read on reference
		int curAlignedRefBase = 0;
		//current read's start postion on reference 
		int curStartPos = curRead -> aln -> core.pos;

		for (int i = 0; i < curRead -> aln -> core.n_cigar; i++) {
			//the number of reads span the valid indel
			int read_span = 0;
			//the number of reads having the same valid indel
			int same_indel = 0;
			//the portion of reads span the indel that have the same indel
			float portion = 0;

			//if consume reference
			if (bam_cigar_type(curCigarPtr[i]) & 0x2 ) {
				curAlignedRefBase += bam_cigar_oplen(curCigarPtr[i]);
			}
			//if there is a large indel
			if (bam_cigar_oplen(curCigarPtr[i]) >= indel_min_size &&
				(bam_cigar_op(curCigarPtr[i]) == BAM_CINS ||
				bam_cigar_op(curCigarPtr[i]) == BAM_CDEL)){
				//indel's position on reference
				indel_pos = curStartPos + curAlignedRefBase;
				//testing
				//cout << indel_pos << endl;
				//indel's size
				indel_size = bam_cigar_oplen(curCigarPtr[i]);

				for (int j = bottom; j <= top; j++) {
					if (j != cur){
						//if this read span the indel
						if (readQueue[j] -> startPos < indel_pos && readQueue[j] -> endPos > indel_pos) {
							read_span ++;
							uint32_t* jCigarPtr = bam_get_cigar(readQueue[j] -> aln);
							int jAlignedRefBase = 0;
							int jStartPos = readQueue[j] -> aln -> core.pos;
							int j_indel_pos = 0;
							for (int k = 0; k < readQueue[j] -> aln -> core.n_cigar; k++) {
								//if consume reference
								if (bam_cigar_type(jCigarPtr[k]) & 0x2 ) {
									jAlignedRefBase += bam_cigar_oplen(jCigarPtr[k]);
								}
								//if there is a same size indel
								if (bam_cigar_oplen(jCigarPtr[k]) >= 0.5 * indel_size &&
									bam_cigar_oplen(jCigarPtr[k]) <= 1.5 * indel_size &&
									(bam_cigar_op(jCigarPtr[k]) == BAM_CINS ||
									bam_cigar_op(jCigarPtr[k]) == BAM_CDEL)){
									//indel's position on reference
									j_indel_pos = jStartPos + jAlignedRefBase;
									if (j_indel_pos > (indel_pos - indel_size*region_len) &&
										j_indel_pos < (indel_pos + indel_size*region_len)) {
										same_indel ++;
										break;
									}
								}
							}	
						}
					}
				}
				portion = float(same_indel) / float(read_span);
				if (portion >= 0.5) {
					//if this is the first valid indel
					if (if_first_indel) {
						//set pre_indel_pos to current position
						pre_indel_pos = indel_pos;
						//output in terminal
						cout << "portion of indels: " << portion 
							<< " position: " << indel_pos
							<< " indel size: " << indel_size << endl;
						//write to the text file
						if (write_error_region) {
							text << "portion of indels: " << portion 
							 	<< " position: " << indel_pos
							 	<< " indel size: " << indel_size << "\n";
						}
						//create a valid_indel object
						Valid_indel curIndel;
						curIndel.portion = portion;
						curIndel.indel_pos = indel_pos;
						curIndel.indel_size = indel_size;
						//add curIndel to the vector of all valid indels
						all_indel_info.push_back(curIndel);
						//after the first valid indel, set if_first_indel be false
						if_first_indel = false;
					}
					else {
						//Try to omit those indels that are too close: should be counted once
						//after checking, if the indel is a new one, add it to all_indel_info
						//Check_new_indel(indel_pos) returns true if it's a new one
						if (Check_new_indel(all_indel_info, indel_pos, range_btw_indels)) {
							//create a valid_indel object
							Valid_indel curIndel;
							curIndel.portion = portion;
							curIndel.indel_pos = indel_pos;
							curIndel.indel_size = indel_size;
							//add curIndel to the vector of all valid indels
							all_indel_info.push_back(curIndel);
						} 
						//TODO: however, the indels may not come in the order of position
						//for simplicity, assume they come in order; fix this later
						if (indel_pos < pre_indel_pos - range_btw_indels ||
							indel_pos > pre_indel_pos + range_btw_indels) {
							//set pre_indel_pos to current position
							pre_indel_pos = indel_pos;
							//output in terminal
							cout << "portion of indels: " << portion 
								<< " position: " << indel_pos
								<< " indel size: " << indel_size << endl;
							//write to the text file
							if (write_error_region) {
								text << "portion of indels: " << portion 
							 		<< " position: " << indel_pos
							 		<< " indel size: " << indel_size << "\n";
							}
						}
					}
					
				}
			}
		}
		
		//increase cur index by 1: moving forward
		cur = cur + 1;

		if (readQueue.find(cur) == readQueue.end()) {
			//cerr << "bottom: " << bottom << " cur: " << cur << " top: " << top << endl;
      		//assert(readQueue.size() == 0);
			//cerr << "At the end of the queue!!! " << cur - 1 << "\t" 
			//	 << curRead -> startPos << "\t" << curRead -> name << endl;
			cerr << "bottom: " << bottom << " cur: " << cur << " top: " << top << endl;
			AlignedRead* curRead= new AlignedRead;
			curRead->index = cur;
			int res=sam_read1(htsfp, samHeader, alignedRead);
			while (alignedRead->core.flag & 2304 && sam_read1(htsfp, samHeader, alignedRead) > 0) {}
			if (res < 0) {
				curRead=NULL;
				break;
			}
		}
	}
	if (write_error_region) {
		//close text
		text.close();
	}

	//testing
	cout << endl << endl << endl;
	cout << "testing vector indels" << endl;
	if (all_indel_info.size() > 0) {
		for (vector <Valid_indel>::iterator it = all_indel_info.begin() ; 
		 	it != all_indel_info.end();
		 	it++) {
			cout << it -> indel_pos << endl;
		}
	}
}



