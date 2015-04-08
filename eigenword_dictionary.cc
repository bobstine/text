#include "eigenword_dictionary.h"

#include <sstream>
#include <fstream>

static std::string messageTag ("EWDC: ");


int
EigenwordDictionary::init_from_file(string fileName)
{ 
  std::ifstream input(fileName);
  if (!input)
  { std::cerr << messageTag << " *** ERROR ***  File " << fileName << " not found.\n";
    return 0;
  }
  std::string line;
  int columnCount = 0;
  {
    getline(input,line);                  // check column count from headings
    std::istringstream is (line);
    std::string name;
    while(is >> name) ++columnCount;
  }
  if ((columnCount-1) != eigen_dim())     // first in header is name of words
  { std::cerr << messageTag << " *** ERROR *** Found " << columnCount << " columns, but expected " << eigen_dim() << std::endl
	      << "      Line read: <<<" << line << ">>>" << std::endl;
    return 0;
  }
  int lineNumber = 0;
  while(std::getline(input,line))
  { std::string::const_iterator it = line.begin();
    if (*it != '"')
      std::cerr << messageTag << "Input eigenword line " << line << " not properly formatted.\n";
    else
    { std::string token;
      int count = 2;
      ++it;
      while (*it != '"' && count < 75)
      { token.push_back(*it);
	++count; 
	++it;
      }
      if (count < 50)
      { mTypes.push_back(Text::Type(token));
	mIndexMap[Text::Type(token)] = lineNumber;
	std::istringstream is (line.substr(count));
	for(int j=0; j<mEigenwords.cols(); ++j)
	  is >> mEigenwords(lineNumber,j);
	++lineNumber;
      }
      else // did not find matching " for end of type within 75 chars
	std::cerr << messageTag << "*** ERROR *** Line around "
		  << lineNumber << " that begins -->" << line.substr(0,30) << "<-- appears incorrectly formatted.\n";
     }
  }
  std::clog << messageTag << "Read " << lineNumber << " lines from file " << fileName << std::endl;
  return lineNumber;
}

//     operator     operator     operator     operator     operator     operator     operator     operator

int
EigenwordDictionary::type_index(Text::Type const& t) const
{
  auto it = mIndexMap.find(t);
  if(it != mIndexMap.end())
    return it->second;
  else
    return -1;
}

EigenwordDictionary::Vector
EigenwordDictionary::operator[](Text::Type const& t) const
{
  int row = type_index(t);
  if(row < 0)
  { std::cerr << messageTag << "Word type " << t << " not in eigenword dictionary; return null vector.\n";
    return Vector::Zero(0);
  }
  else
    return mEigenwords.row(row);
}

EigenwordDictionary::Matrix
EigenwordDictionary::eigenwords(TypeVector const& types) const
{
  Matrix dict(types.size(), mEigenwords.cols());

  for(int row=0; row<(int)types.size(); ++row)
  { int index = type_index(types[row]);
    if (-1 < index)
      dict.row(row) = operator[](index);
    else
      dict.row(row) = Vector::Zero(dict.cols());
  }
  return dict;
}

//     print     print     print     print     print     print     print     print     print     print     print

void
EigenwordDictionary::print_to_stream(std::ostream &os) const
{
  os << "EigenwordDictionary with " << n_types() << " types and " << mEigenwords.cols() << " dimensions. Leading types are:";
  for (int i=0; i<10; ++i)
    os << " " << mTypes[i];
}

