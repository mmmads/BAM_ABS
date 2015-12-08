// ==========================================================================
//                             BAM_ABS
// ==========================================================================

//
// ==========================================================================
// Author: Saima Sultana Tithi <saima5@vt.edu>, Hong Tran <hongt1@vt.edu>
// ==========================================================================

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/modifier.h>
#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <math.h>
#include <map>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

using namespace std;

typedef multimap<string, string> AmbiguousMap;
typedef AmbiguousMap::iterator mapAmbIter;

typedef multimap<string, string> UniqueMap;
typedef UniqueMap::iterator mapUniIter;

double pAC = 0.00025;
double pAT = 0.00025;
double pAG = 0.0005;
double pCA = 0.00025;
double pCT = 0.0005;
double pCG = 0.00025;
double pTC = 0.0005;
double pTA = 0.00025;
double pTG = 0.00025;
double pGA = 0.0005;
double pGT = 0.00025;
double pGC = 0.00025;
double pSNP = 0.001;
double CG = 0.85;
double CH = 0.02;

char getRevComp(char c)
{
    switch(c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        default:
            return c;
    }
}

double convertPhred33(char c)
{
    int x = c - 33;
    return pow(10.0, (-x/10.0));
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

void prior(string seq, string genomicSeq, double phredArray[], double priorArr[])
{
	for (int i = 0; i < seq.length(); i++)
    {
		priorArr[i] = 0.0;
	}

    for (int i = 0; i < seq.length(); i++)
    {
        if (phredArray[i] > 0)
        {
            if (genomicSeq[i] == 'A')
            {
                if (seq[i] == 'A')
                    priorArr[i] = 1 - pSNP;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pAC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pAC * CG;
                    else
                        priorArr[i] = pAC * CH;
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pAG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pAT + pAC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pAT + pAC * (1-CG);
                    else
                        priorArr[i] = pAT + pAC * (1-CH);
                }
            }
            else if (genomicSeq[i] == 'C')
            {
                if (seq[i] == 'A')
                    priorArr[i] = pCA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) * CG;
                    else
                        priorArr[i] = (1-pSNP) * CH;
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pCG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pCT + (1-pSNP) * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pCT + (1-pSNP) * (1-CG);
                    else
                        priorArr[i] = pCT + (1-pSNP) * (1-CH);
                }
            }
            else if (genomicSeq[i] == 'G')
            {
                if (seq[i] == 'A')
                    priorArr[i] = pGA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pGC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pGC * CG;
                    else
                        priorArr[i] = pGC * CH;
                }
                else if (seq[i] == 'G')
                    priorArr[i] = 1 - pSNP;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pGT + pGC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pGT + pGC * (1-CG);
                    else
                        priorArr[i] = pGT + pGC * (1-CH);
                }
            }
            else //T at genomic location
            {
                if (seq[i] == 'A')
                    priorArr[i] = pTA;
                else if (seq[i] == 'C')
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = pTC * CG;
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = pTC * CG;
                    else
                        priorArr[i] = pTC * CH;
                }
                else if (seq[i] == 'G')
                    priorArr[i] = pTG;
                else
                {
                    if ( (i+1) < seq.length() && seq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) + pTC * (1-CG);
                    else if ( (i+1) == seq.length() && genomicSeq[i+1] == 'G')
                        priorArr[i] = (1-pSNP) + pTC * (1-CG);
                    else
                        priorArr[i] = (1-pSNP) + pTC * (1-CH);
                }
            }
        }
        else
		{
            priorArr[i] = 1;
		}
    }
}

double sum(double observedProb[], int size)
{
    double sum = 0.0;
    for(int i = 0; i < size; i++)
    {
        sum += observedProb[i];
    }
    return sum;
}

int main(int argc, char const ** argv)
{
	cout << "Start processing...\n";
	if (argc != 4)
    {
        cout << "USAGE: ./main file.fa ambiguous_read_file unique_overlap_read_file\n";
        return 1;
    }

	// Try to load index and create on the fly if necessary.
	seqan::FaiIndex faiIndex;
	if (seqan::read(faiIndex, argv[1]) != 0)
	{
		if (build(faiIndex, argv[1]) != 0)
		{
		    cout << "ERROR: Fai index file could not be loaded or built.\n";
		    return 1;
		}
		if (write(faiIndex) != 0)  // Name is stored from when reading.
		{
		    cout << "ERROR: Fai index file could not be written to disk.\n";
		    return 1;
		}
	}

    string line;
	double rate[15];
	int lineNum = 0;
	ifstream myfile("snp_methylation_rate.txt");
    if (myfile.is_open()) {
        while (getline(myfile, line)  && lineNum <= 14) {
            vector<string> tokens = split(line, '=');
			double rateVal = atof(tokens.at(1).c_str());
			rate[lineNum] = rateVal;
			lineNum++;
        }
        myfile.close();
    }
	else {
		cout << "ERROR: Snp and methylation rate file could not be loaded.\n";
	    return 1;
	}

	if(rate[0] <= 0 || rate[0] >=1) {
		pAC = 0.00025;
		cout << "ERROR: pAC could not be read. The default value is used.\n";
	}
	else
		pAC = rate[0];
	if(rate[1] <= 0 || rate[1] >=1) {
		pAT = 0.00025;
		cout << "ERROR: pAT could not be read. The default value is used.\n";
	}
	else
		pAT = rate[1];
	if(rate[2] <= 0 || rate[2] >=1) {
		pAG = 0.0005;
		cout << "ERROR: pAG could not be read. The default value is used.\n";
	}
	else
		pAG = rate[2];
	if(rate[3] <= 0 || rate[3] >=1) {
		pCA = 0.00025;
		cout << "ERROR: pCA could not be read. The default value is used.\n";
	}
	else
		pCA = rate[3];
	if(rate[4] <= 0 || rate[4] >=1) {
		pCT = 0.0005;
		cout << "ERROR: pCT could not be read. The default value is used.\n";
	}
	else
		pCT = rate[4];
	if(rate[5] <= 0 || rate[5] >=1) {
		pCG = 0.00025;
		cout << "ERROR: pCG could not be read. The default value is used.\n";
	}
	else
		pCG = rate[5];
	if(rate[6] <= 0 || rate[6] >=1) {
		pTC = 0.0005;
		cout << "ERROR: pTC could not be read. The default value is used.\n";
	}
	else
		pTC = rate[6];
	if(rate[7] <= 0 || rate[7] >=1) {
		pTA = 0.00025;
		cout << "ERROR: pTA could not be read. The default value is used.\n";
	}
	else
		pTA = rate[7];
	if(rate[8] <= 0 || rate[8] >=1) {
		pTG = 0.00025;
		cout << "ERROR: pTG could not be read. The default value is used.\n";
	}
	else
		pTG = rate[8];
	if(rate[9] <= 0 || rate[9] >=1) {
		pGA = 0.0005;
		cout << "ERROR: pGA could not be read. The default value is used.\n";
	}
	else
		pGA = rate[9];
	if(rate[10] <= 0 || rate[10] >=1) {
		pGT = 0.00025;
		cout << "ERROR: pGT could not be read. The default value is used.\n";
	}
	else
		pGT = rate[10];
	if(rate[11] <= 0 || rate[11] >=1) {
		pGC = 0.00025;
		cout << "ERROR: pGC could not be read. The default value is used.\n";
	}
	else
		pGC = rate[11];
	if(rate[12] <= 0 || rate[12] >=1) {
		pSNP = 0.001;
		cout << "ERROR: pSNP could not be read. The default value is used.\n";
	}
	else
		pSNP = rate[12];
	if(rate[13] <= 0 || rate[13] >=1) {
		CG = 0.85;
		cout << "ERROR: pCG could not be read. The default value is used.\n";
	}
	else
		CG = rate[13];
	if(rate[14] <= 0 || rate[14] >=1) {
		CH = 0.02;
		cout << "ERROR: pCH could not be read. The default value is used.\n";
	}
	else
		CH = rate[14];
	//cout << pAC << "," << pAT << "," << pAG << "," << pCA << "," << pCT << "," << pCG << "," << pTC << "," << pTA << "," << pTG << "," << pGA << "," << pGT << "," << pGC << "," << pSNP << "," << CG << "," << CH << endl;

	AmbiguousMap myAmbMap;
    ifstream myfileAmb(argv[2]);
    if (myfileAmb.is_open()) {
        while (getline(myfileAmb, line)) {
            string::size_type pos = line.find_first_of('\t');
            string key = line.substr(0, pos);
            string value = line.substr(pos);
            value.erase(0, value.find_first_not_of('\t'));
            myAmbMap.insert(make_pair(key, value));
        }
        myfileAmb.close();
    }
	else {
		cout << "ERROR: Ambiguous read file could not be loaded.\n";
	    return 1;
	}

	/*
	int MAPQ = atoi(argv[4]);
	if(MAPQ < 0 || MAPQ > 255) {
		cout << "MAPQ value is not between 0 and 255. The default value 30 is used.\n";
		MAPQ = 30;
	}*/

	UniqueMap myUniMap;
    ifstream myfileUni(argv[3]);
    if (myfileUni.is_open()) {
        while (getline(myfileUni, line)) {
            vector<string> tokens = split(line, '\t');
			if(tokens.at(4).compare(tokens.at(9))) {
		        string key = tokens.at(3) + "_" + tokens.at(0) + "_" + tokens.at(1);
		        string value = tokens.at(4) + "\t" + tokens.at(6) + "\t" + tokens.at(10)
		                + "\t" + tokens.at(11) + "\t" + tokens.at(12) + "\t" + tokens.at(13)
		                + "\t" + tokens.at(14) + "\t" + tokens.at(15) + "\t" + tokens.at(16);
		        myUniMap.insert(make_pair(key, value));
			}
        }
        myfileUni.close();
    }
	else {
		cout << "ERROR: Unique overlap read file could not be loaded.\n";
	    return 1;
	}

	//ofstream outputfile;
  	//outputfile.open ("All_reads_with_probability.txt");
	ofstream outputfile1;
  	outputfile1.open ("Reads_with_highest_probable_location.sam");

	//string header = "#All multi-reads in sam format with probability(999 means probability cannot be calculated)\n";
	//string header1 = "#Multi-reads having the highest probability\n";

	//outputfile << header;
	//outputfile1 << header1;

    mapAmbIter m_it, s_it;
    mapUniIter mUni_it, sUni_it;

    for (m_it = myAmbMap.begin();  m_it != myAmbMap.end();  m_it = s_it)
    {
        string theKey = (*m_it).first;
        pair<mapAmbIter, mapAmbIter> keyRange = myAmbMap.equal_range(theKey);

        // Iterate over all map elements with key == theKey
		vector<double> likelihoodArr;
		std::map<double, std::string> likelihoodMap;
		string output = "";
        for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it)
        {
            string theValue = (*s_it).second;
            vector<string> tokens = split(theValue, '\t');
            unsigned flag = atoi(tokens.at(0).c_str());
            string seq = tokens.at(8);
            string phred33 = tokens.at(9);

			double multiSeqErr[phred33.length()];

            for (int i = 0; i < phred33.length(); i++)
            {
               multiSeqErr[i] = convertPhred33(phred33[i]);
            }

			//get the genome seq using seqan library
            string genomicSeq;

			// Translate sequence name to index.
			unsigned idx = 0;
			if (!getIdByName(faiIndex, tokens.at(1), idx))
			{
				if (!getIdByName(faiIndex, "chr" + tokens.at(1), idx))
                {
                    std::cerr << "ERROR: FAI index file has no entry for chromosome " << tokens.at(1) << "\n";
                    return 1;
                }
			}
			//unsigned seqLength = sequenceLength(faiIndex, idx);
			//cout << seqLength << "\n";

			// Load characters of chromosome
			seqan::CharString seqChrPrefix;
			unsigned startPos = atoi(tokens.at(2).c_str()) - 2;
			unsigned endPos = startPos + seq.length() + 2;

			if (readRegion(seqChrPrefix, faiIndex, idx, startPos, endPos) != 0)
			{
				std::cerr << "ERROR: Could not load reference sequence.\n";
				return 1;
			}
			//cout << "Seq:" << seqChrPrefix << "\n";

			std::stringstream stream;
			stream<<seqChrPrefix;

			genomicSeq = stream.str();
			//std::cout << genomicSeq.length() << "\n";

			if(flag == 16 || flag == 272) {
				string revGenomicSeq(genomicSeq);
				for (int i = 0; i < genomicSeq.length(); i++)
                {
                   revGenomicSeq[genomicSeq.length() - 1 -i] = getRevComp(genomicSeq[i]);
                }
				genomicSeq = revGenomicSeq;
			}
			genomicSeq = genomicSeq.substr(1);

			double priorArr[seq.length()];
            prior(seq, genomicSeq, multiSeqErr, priorArr);

            if(flag == 16 || flag == 272) {
    			double revPriorArr[seq.length()];
				for (int i = 0; i < seq.length(); i++)
                {
                   revPriorArr[seq.length() - 1 -i] = priorArr[i];
                }
				for (int i = 0; i < seq.length(); i++)
                {
                   priorArr[i] = revPriorArr[i];
                }

				double revMultiSeqErr[phred33.length()];
				for (int i = 0; i < phred33.length(); i++)
                {
                   revMultiSeqErr[phred33.length() - 1 -i] = multiSeqErr[i];
                }
				for (int i = 0; i < phred33.length(); i++)
                {
                   multiSeqErr[i] = revMultiSeqErr[i];
                }

                for (int i = 0; i < seq.length(); i++)
                {
                   seq[tokens.at(8).length() - 1 -i] = getRevComp(tokens.at(8)[i]);
                }
            }

            double observed[seq.length()];
			//initialize observed array
			for(int a = 0; a < seq.length(); a++)
			{
				observed[a] = 0.0;
			}

            double posterior[seq.length()];
			//initialize posterior array
			for(int a = 0; a < seq.length(); a++)
			{
				posterior[a] = 0.0;
			}

            string theUniKey = theKey + "_" + tokens.at(1) + "_" + tokens.at(2);
            if (myUniMap.count(theUniKey) != 0)
            {
                pair<mapUniIter, mapUniIter> keyRangeUni = myUniMap.equal_range(theUniKey);

                for (int i = 0; i < seq.length(); i++)
                {
                    if (multiSeqErr[i] > 0)
                    {
                        double observedProb[seq.length()];
						//initialize observedProb array
						for(int a = 0; a < seq.length(); a++)
						{
							observedProb[a] = 0.0;
						}

                        int baseTotal = 0;
                        // Iterate over all map elements with key == theUniKey
                        for (sUni_it = keyRangeUni.first;  sUni_it != keyRangeUni.second;  ++sUni_it)
                        {
                            string theUniValue = (*sUni_it).second;
                            vector<string> uniTokens = split(theUniValue, '\t');

                            double uniSeqErr[uniTokens.at(8).length()];

                            for (unsigned x = 0; x < uniTokens.at(8).length(); x++)
                            {
                               uniSeqErr[x] = convertPhred33(uniTokens.at(8)[x]);
                            }
                            unsigned ambPos = atoi(tokens.at(2).c_str());
                            unsigned uniPos = atoi(uniTokens.at(1).c_str());
                            if(uniPos <= (ambPos+i) && (uniPos+uniTokens.at(7).length()-1) >= (ambPos+i))
                            {
                                if(uniTokens.at(7).at(ambPos+i-uniPos) == seq[i])
                                {
                                    observedProb[baseTotal] = 1 - uniSeqErr[ambPos+i-uniPos]
                                            - multiSeqErr[i] + (uniSeqErr[ambPos+i-uniPos] * multiSeqErr[i]);
                                }
                                else
                                {
                                    observedProb[baseTotal] = uniSeqErr[ambPos+i-uniPos]
                                            + multiSeqErr[i] - (uniSeqErr[ambPos+i-uniPos] * multiSeqErr[i]);
                                }
                                baseTotal++;
                            }
                        }
                        if (baseTotal != 0)
                        {
                            observed[i] = sum(observedProb, seq.length()) / baseTotal;
                        }
                        else
                        {
                            observed[i] = 0.999999999999999;
                        }
						double value = (priorArr[i] * observed[i]) / (priorArr[i] * observed[i] + (1-priorArr[i]) * (1-observed[i]));
						//std::cout << std::setprecision(15) << value << "\n";

						if (value > 0)
							posterior[i] = log10(value);
                    }
                    else
                    {
                        posterior[i] = 0.0;
                    }
                }//end for
				double likelihood = sum(posterior, seq.length());
				//output = theKey + "\t" + tokens.at(1) + "_" + tokens.at(2) + "\t";
				output = theKey + "\t" + theValue + "\t";
				//outputfile << output;
				//outputfile << std::setprecision(15) << likelihood << "\n";
				likelihoodArr.push_back(likelihood);
				likelihoodMap[likelihood] = output;
            }
            else
            {
				//output = theKey + "\t" + tokens.at(1) + "_" + tokens.at(2) + "\t";
				output = theKey + "\t" + theValue + "\t";
				//outputfile << output;
				//outputfile << 999 << "\n";
            }

        }//end for
		if(!likelihoodArr.empty()) {
			vector<double>::const_iterator it;
			it = std::max_element(likelihoodArr.begin(), likelihoodArr.end());
			std::map<double, std::string>::const_iterator search = likelihoodMap.find(*it);
			outputfile1 << search->second << "\n";
			//outputfile1 << std::setprecision(15) << *it << "\n";
		}
    }//end for

	//outputfile.close();
	outputfile1.close();
	cout << "Output is written in Reads_with_highest_probable_location.sam\n";

    return 0;
}
