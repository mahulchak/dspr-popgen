//this program takes a GFF file and reconstructs the CDS. Then it takes the codon table and assigns amino acids to each
//codon in the CDS. if a VCF file is provided, this program will inform the user whether an SNP is synonymous or non-synonymous. finally, the program informs the 4 four degenrate sites.

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<algorithm>
#include<vector>
#include<map>

using namespace std;

vector<string> parentSplitter(string & str);//split the last column of the GFF
vector<string> fieldSplitter(string & str);//splits the row into individual columns

int main()
{
	ifstream fin;
	ofstream fout;
	fout.open("mycds.txt");
	string str;
	vector<string> fields,vs;
	vector<int> cord(2);//records start and end of a feature
	map<string,vector<int> > transSpan;//records transcript span
//vector<int> has 3 elements. 1st is start, 2nd is end, 3rd is 0 or 1. 0 is + strand and 1 is - strand
	fin.open("dmel-all-filtered-r6.09.gff");
	while(getline(fin,str))
	{
		if((str.find('#') == string::npos) && (str.find('>') == string::npos))//comments and sequences are ignored
		{
			//cout<<"Line found\t"<<str<<endl;
			fields = fieldSplitter(str);		
			cord[0] = stoi(fields[3]);
			cord[1] = stoi(fields[4]);
			//cout<<"start=\t"<<fields[3]<<"\t"<<"end=\t"<<fields[4]<<endl;
		}
		if(str.find("\tCDS") != string::npos)//when the line has CDS feature
		{
			fout<<str<<endl;
			vs = parentSplitter(str);
			for(unsigned int i=0;i<vs.size();i++)
			{
				cout<<vs[i]<<'\t'<<cord[0]<<'\t'<<cord[1]<<endl;
			}
		}
	}	
	fin.close();
	fout.close();	
	return 0;
}
vector<string> fieldSplitter(string & str)
{
	vector<string> vs;
	size_t pos1,pos2;
	string field;
	pos1 = 0;
	while(pos2 != string::npos)//until no more tabs are found
	{
		pos2 = str.find('\t',pos1);
	//cout<<pos2<<endl;
		if(pos2 == string::npos)//if last field
		{
			field = str.substr(pos1);
		}
		else
		{
			field = str.substr(pos1,pos2-pos1);
		}
		pos1 = pos2+1;
		cout<<field<<endl;
		vs.push_back(field);
	}
	return vs;
}
	
vector<string> parentSplitter(string & str)
{
	string col5, transName;
	size_t pos1 =0, pos2 =0;
	vector<string> vs;
	pos1 = str.find("Parent");//finds where the col starts
	pos2 = str.find('=',pos1+1);
	col5 = str.substr(pos2+1);
	pos1 = 0;
	//cout<<col5<<'\t'<<pos1<<'\t'<<pos2<<endl;
	if(col5.find(',') == string::npos)//no comma is present
	{
		transName = col5.substr(pos1);
		vs.push_back(transName);
		//cout<<transName<<endl;
	}
	else //if comma is present/>1 transcript
	{	
		while(pos2 != string::npos)
		{
			pos2 = col5.find(',',pos1+1);//look for comma
			if(pos2 == string::npos)//if no comma is found
			{
				//cout<<"came here 1\t"<<col5<<'\t'<<pos1<<endl;
				transName = col5.substr(pos1);
				vs.push_back(transName);
				cout<<transName<<endl;
			}
			if(pos2 != string::npos)//comma is found
			{
			//	cout<<"came here 2"<<endl;
				transName = col5.substr(pos1,pos2-pos1);	 	
				vs.push_back(transName);
				cout<<transName<<endl;
			}
			pos1 = pos2 +1;
		}
	}
	return vs;
}

