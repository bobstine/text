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
  { getline(input,line);             // check column count from headings
    std::istringstream is (line);
    std::string name;
    int count = 0;
    while(is >> name) ++count;
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
      while (*it != '"' && count < 50)
      { token.push_back(*it);
	++count; 
	++it;
      }
      if (count < 50)
      { mTypes.push_back(Type(token));
	mIndexMap[Type(token)] = lineNumber;
	std::istringstream is (line.substr(count));
	for(int j=0; j<mEigenwords.cols(); ++j)
	  is >> mEigenwords(lineNumber,j);
	++lineNumber;
      }
      else // did not find matching " for end of type within 50 chars
	std::cerr << messageTag << "*** ERROR *** Line around " << lineNumber << " that begins -->" << line.substr(0,30) << "<-- appears incorrectly formatted.\n";
     }
  }
  std::clog << messageTag << "Read " << lineNumber << " lines from file " << fileName << std::endl;
  return lineNumber;
}

//     operator     operator     operator     operator     operator     operator     operator     operator

int
EigenwordDictionary::type_index(Type const& t) const
{
  auto it = mIndexMap.find(t);
  if(it != mIndexMap.end())
    return it->second;
  else
    return -1;
}

EigenwordDictionary::Vector
EigenwordDictionary::operator[](Type const& t) const
{
  int row = type_index(t);
  if(row < 0)
  { std::cerr << messageTag << "*** ERROR *** Type " << t << " not found; return 1 vector.\n";
    return Vector::Ones(mEigenwords.cols());
  }
  else
    return mEigenwords.row(row);
}

EigenwordDictionary::Matrix
EigenwordDictionary::eigenwords(TypeVector const& types) const
{
  Matrix dict(types.size(), mEigenwords.cols());

  for(int row=0; row<(int)types.size(); ++row)
    dict.row(row) = operator[](types[row]);
  return dict;
}

//     print     print     print     print     print     print     print     print     print     print     print

void
EigenwordDictionary::print_to_stream(std::ostream &os) const
{
  os << "EigenwordDictionary with " << n_types() << " types and " << mEigenwords.cols() << " dimensions. Leading types are:";
  for (int i=0; i<5; ++i)
    os << " " << mTypes[i];
}

