#ifndef _REGEX_H_
#define _REGEX_H_

#include <string>

/*
  Temporary regex for converting

       numbers to log scale rounding
       dollars to log scale
       phone numbers to area codes

*/

class Parser: std::unary_function<std::string, float>
{
  const std::string   mName;                 // pure function with no state
  float             (*mFunc)(std::string); 

 public:
  
 Parser( std::string name, float (*f)(std::string) ): mName(name), mFunc(f) {} 

  std::string name()               const { return mName; }
  
  float operator()(std::string s)  const { return mFunc(s); }
};


float square_footage   (std::string str);     // all return 0 if not found; only find first in string
float number_bathrooms (std::string str);
float number_bedrooms  (std::string str);


bool         can_parse_as_numeric_string(std::string token);

std::string  parse_numeric_string(std::string token);      // assumes you have tested with prior function

#endif
