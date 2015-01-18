/*
  Converts the preposition lines in Joel format into a rectangular array
  that will then be processed to give word2vec coordinates.  Exports only those
  columns identified the tags file.

  Joel style (space delimited)
       into DY BGL#hurry_VBD BGR#the_ATI FFtag#, FFword#, FH#classroom_NN FHtag#NN FHword#classroom
       around DY BGL#look_VBD BGR#the_ATI FFtag#IN FFword#for PV#look_VBD

  Rectangular style (tab delimited, with uniform set of columns)
       Y        BGL          BGR       FFtag  FFword   FH            FHtag    FHword     PV
       into     hurry_VBD    the_ATI   ,      ,        classroom_NN  NN       classroom  NA
       around   look_VBD     the_ATI   IN     for      [blank]       NA        NA     look_VBD

*/

#include <getopt.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <map>
#include <set>

#include "string_trim.h"

const bool verbose = false;

void
parse_arguments(int argc, char** argv, bool &keepPOS, std::string& tagFileName);


int main(int argc, char** argv)
{
  using std::string;
  using std::set;
  using std::map;

  // defaults
  string tagFileName ("tags.txt");
  bool   keepPOS     (false);
  parse_arguments(argc, argv, keepPOS, tagFileName);

  // read set of tags
  std::ifstream tagStream (tagFileName.c_str(), std::ifstream::in);
  set<string> tags;
  if(tagStream)
  { while(!tagStream.eof())
    { string tag;
      tagStream >> tag;
      tag = trim(tag);
      if(tag.size() > 0) tags.insert(tag);
    }
  }
  else
  { std::cerr << "ERROR: Tag file `" << tagFileName << "' not found; terminating.\n";
    return 0;
  }
  std::clog << "--- Read " << tags.size() << " tags from file '" << tagFileName << "'.\n";

  // write tags as column headers
  std::cout << "Y";
  for(string tag:tags) std::cout << "\t" << tag;
  std::cout << std::endl;

  // process each line to std output
  while (!std::cin.eof())
  { string thePrep;
    string dlTag;
    std::cin >> thePrep >> dlTag;
    string theLine;
    std::getline(std::cin, theLine);
    trim(theLine);
    if(theLine.size() > 0)
    { if (verbose) std::clog << "\n\nRead prep=" << thePrep << " for line:\n" << theLine << std::endl;
      std::istringstream ss (theLine);
      map<string,string> tagValues;
      while(!ss.eof())
      { string token;
	ss >> token;
	token = trim(token);
	std::size_t pos = token.find('#');
	if(pos != string::npos)
	{ string key    = token.substr(0,pos);
	  if(tags.find(key) != tags.end())
	  { string value  = token.substr(pos+1);
	    if(!keepPOS)
	    { std::size_t p = value.find('_');
	      if(p != string::npos)
		value = value.substr(0,p);
	    }
	    if (verbose) std::clog << "    Key " << key << "=" << value << std::endl;
	    tagValues[key] = value;
	  }
	}
      }
      // write the line
      std::cout << thePrep;
      for (string tag : tags)
      { string value = tagValues[tag];
	if(value.size() == 0)
	  std::cout << "\t" << "NA";
	else
	  std::cout << "\t" << value;
      }
      std::cout << std::endl;
    }
  }
}

void
parse_arguments(int argc, char** argv, bool &keepPOS, std::string& tagFileName)
{
  static struct option long_options[] = {
    //    {"risk",               no_argument, 0, 'R'},
    {"keep_POS",     no_argument,       0, 'p'},
    {"tag_file",     required_argument, 0, 't'},
    {0, 0, 0, 0}                             // terminator
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "pt:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "PARSE: Key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'p' :  { keepPOS = true; 	break;      }
    case 't' :  { tagFileName = optarg; break;      }
      /*    case 'n' :
      {
	nRounds = read_utils::lexical_cast<int>(optarg);
	break;
      }
      */
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}

