#include "regex.h"

#include <iostream>
#include <string>
#include <sstream>
#include <regex>


int main ()
{
  using std::string;
  using std::cout;
  using std::endl;


  // Parse for square footage, bedrooms, bathrooms
  string desc;

  desc = "3 bdrm   2bth  1800sf";
  cout << "TEST: Input " << desc << " gives "
       << number_bedrooms(desc)  << " bedrooms, " << number_bathrooms(desc) << " bathrooms, and " << square_footage(desc) << " square feet.\n";
  desc = "3 beds   2.5bths  1,800 sq ft";
  cout << "TEST: Input " << desc << " gives "
       << number_bedrooms(desc)  << " bedrooms, " << number_bathrooms(desc) << " bathrooms, and " << square_footage(desc) << " square feet.\n";

  cout << endl << endl;

  // wrapped in parser object
  Parser nBedrmFunc ("bedrooms", number_bedrooms);
  cout << "TEST: Wrapped as parser gives " << nBedrmFunc(desc) << " bedrooms " << endl;

  const std::vector<Parser> parsers { Parser("SqFt", square_footage), Parser("Bedrooms", number_bedrooms), Parser("Bathrooms", number_bathrooms) };
  std::cout << "TEST: vector of parsers...  using string    " << desc << endl << "         ";
  for (auto f : parsers)
    cout << f(desc) << " " << f.name() << "   ";
  cout << endl;


  desc = "3 bsdfdrm   2ssdbth  1800sdf";
  std::cout << "TEST: vector of parsers...  using string    " << desc << endl << "         ";
  for (auto f : parsers)
    cout << f(desc) << " " << f.name() << "   ";
  cout << endl;


    
  
  // remove OOV by encoding numeric fields, phone numbers, etc
  //  string     input         ("Axxxx asmdjb 3jdh 333-333-3333 234-087 46464 $100  $100k  $123,123.12  $98.62");

  string     input         ("4495000 the ultimate single family ! this magnificent large - scale residence overlooking the fountains of historic union park has had no expense spared . state of the art systems include central sound and video systems , a luxurious master suite with sitting / dressing room , amazing eat - in kitchen opening to a deck and garden , dining-room with butler's pantry , five large bedrooms , two studies , family / media room , home gym , sauna , roof terrace , wine-cellar , au - pair suite and garage parking for two . a must - see !");

  cout << "TEST: input is " << input << endl;

  //  std::regex dollarPattern (R"(\$(\d+[,\d+][.]?\d+)");
  /*  GCC implementation of regex not ready
      std::regex dollarPattern ("[:digit:]*");
  
      std::smatch sm;
      if(std::regex_match(token,dollarPattern))
  */
  
  string token;
  std::istringstream is (input);
  while(is >> token)
  { cout << "Token `" << token << "' -> ";
    if (can_parse_as_numeric_string(token))
      cout << " parses as " << parse_numeric_string(token) << std::endl;
    else
      cout << " does not match " << endl;
  }
}
  
