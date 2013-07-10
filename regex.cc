#include "regex.h"

// #include <iostream> // debug only
#include <sstream>

#include <boost/regex.hpp>

float
first_numeric_pattern (std::string s, boost::regex const& pattern)
{
  boost::smatch match;
  if (boost::regex_search(s, match, pattern))  // first occurence
  { std::string number;  // remove ,
    std::string str = match[1];
    for (std::string::const_iterator it=str.begin(); it!=str.end(); ++it)
      if (std::isdigit(*it) || (*it == '.'))
	number.push_back(*it);
    std::istringstream ss(number);
    float x;
    ss >> x;
    return x;
  }
  else return 0;
}
  
      
float
square_footage (std::string str)
{
  boost::regex sqftPattern ("(\\d+(,\\d{3})*)((sq?ft?)|( sq(uare)? ?f(ee)?t))");
  return first_numeric_pattern(str,sqftPattern);
}

float
number_bedrooms (std::string str)
{
  boost::regex  bdrmPattern ("(\\d(\\.\\d)?) ?((bd?rm?s?)|(bed(room)?s?))");
  return first_numeric_pattern(str,bdrmPattern);
}

float
number_bathrooms (std::string str)
{
  boost::regex  bathPattern ("(\\d(\\.\\d)?) ?((baths?)|bths?)");
  return first_numeric_pattern(str,bathPattern);
}


//     parse_numeric     parse_numeric     parse_numeric     parse_numeric     parse_numeric     parse_numeric

bool can_parse_as_numeric_string(std::string token)
{
  std::string::const_iterator it = token.begin();
  if (it == token.end())
    return false;
  if (*it == '$')
    ++it;
  else if (!std::isdigit(*it))                               // first must be $ or digit
    return false;
  while (it != token.end() && (std::isdigit(*it)
			   || (*it == '-')                   // used in phone number
			   || (*it == ',') || (*it == '.'))  // embedded in number
	 ) ++it;
  return !token.empty() && (it == token.end()
			    || ((it+1) == token.end() && *it == 'k')) ; // last char is k
}


std::string parse_numeric_string(std::string token)
{
  std::string prefix = "";
  std::string result = "";
  
  { std::string::const_iterator it = token.begin();
    if ( *it == '$' )
    { prefix = "USD";
      ++it;
    }
    while (it != token.end())
    { if (std::isdigit(*it))
      { result.push_back(*it);
	++it;
      }
      else
      { if (*it == '-')   // phone number
	{ prefix = "PH";
	  if (3 == result.size())
	    return prefix + result;
	  else
	    return prefix;
	}
	else if (*it == ',')  // skip over comma
	  ++it;
	else if (*it == 'k')  // add zeros
	{ result = result + "000";
	  break;
	}
	else if (*it == '.')  // decimals are truncated
	  break;
      }
    }
  }
  // encode the resulting value
  if (result.size() > 1)
  { std::string::iterator it = result.begin();
    *it = '1'; ++it;
    while (it != result.end())
    { *it = '0'; ++it;
    }
  }
  return prefix + result;
}
