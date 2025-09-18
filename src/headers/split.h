//
//            Function to perform perl-style splitting of lines of text
//                   Copyright (C) 2010-2011 Matthew Penny
//                           mpenny@jb.man.ac.uk
//

#ifndef SPLIT_HEADER

#include<string>
#include<vector>
#include<sstream>

#include <algorithm> 
#include <cctype>
#include <locale>

using namespace std;

// Trim from the start (in place)
inline void ltrim(string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !isspace(ch);
    }));
}

// Trim from the end (in place)
inline void rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !isspace(ch);
    }).base(), s.end());
}

inline void trim(string &s) {
    rtrim(s);
    ltrim(s);
}

//use getline(ifstream, string) to grab line from file
//use split(string, vector) to split that line into parts stored in the vector

//do not use for T=string split by anything other than whitespace e.g. :
//ignores columns that do not fit format based on the istream conventions

template<class T> void split(const string& line, vector<T>& data, string splitat=" \t")
{
  data.clear();

  int start, end=-1;
  istringstream sub;
  T datum;

  do
    {
      start = end+1;
      end = line.find_first_of(splitat,start);

      if( (end-start>0) || (end==int(string::npos) && start<int(line.size())))
	{
	  //we have some content
	  sub.clear();
	  sub.str(line.substr(start,end-start));
	  sub >> datum;
	  if(!sub.fail())
	    data.push_back(datum);
	  else
	    sub.clear();
	}

    } while(end!=int(string::npos));
}

#define SPLIT_HEADER
#endif
