#include "eigenword_dictionary.h"
#include "file_utils.h"


int main(int, char**)
{
  std::string fileName ("text_src/eigenwords/test_dict.txt");            // head -n 11 google.txt in text_src/eigenwords/
  int cols = FileUtils::count_fields(fileName);
  std::cout << "TEST: Found " << cols << " columns in file " << fileName << std::endl;
  
  EigenwordDictionary dict(fileName, 10, cols-1);
  std::cout << "TEST: " << dict << std::endl;

  std::string s;
  EigenwordDictionary::TypeVector types;

  s = "<s>"; types.push_back(Type(s));
  std::cout << "TEST: index of " << s << " is "  << dict.type_index(Type(s))
	    << " with leading coordinates  " << dict[Type(s)].transpose().head(5) << std::endl;
  s = "the"; types.push_back(Type(s));
  std::cout << "TEST: index of " << s << " is "  << dict.type_index(Type(s))
	    << " with leading coordinates  " << dict[Type(s)].transpose().head(5) << std::endl;
  s = "sss"; types.push_back(Type(s));
  std::cout << "TEST: index of " << s << " is "  << dict.type_index(Type(s))
	    << " with leading coordinates  " << dict[Type(s)].transpose().head(5) << std::endl;

  std::cout << "TEST: Eigenwords for type list:\n" << dict.eigenwords(types) << std::endl;

  
  std::cout << "\n\nTEST: Done.\n";
}
