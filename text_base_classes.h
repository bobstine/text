#ifndef _TEXT_BASE_CLASSES_H_
#define _TEXT_BASE_CLASSES_H_

#include <string>
#include <iostream>

namespace Text {

  using std::string;

  typedef float Scalar;
  
  class POS: public string
  {  
  public:
    
    explicit POS()         : string("") { }
    explicit POS(string s) : string(s) { }
    explicit POS(char *s)  : string(s) { }
    
    void  print_to_stream(std::ostream &os) const   { os << "|" << static_cast<string>(*this) << "|"; }
  };
  
  
  class Type: public string
  {
  public:
    explicit Type()         : string("") { }
    explicit Type(string s) : string(s)       { }
    explicit Type(char *s)  : string(s)       { }
    
    void  print_to_stream(std::ostream &os) const   { os << "\"" << static_cast<string>(*this) << "\"" ; }
  };
}
  
inline
std::ostream&
operator<<(std::ostream &os, Text::POS const& p) { p.print_to_stream(os); return os; }


inline
std::ostream&
operator<<(std::ostream &os, Text::Type const& t) { t.print_to_stream(os); return os; }


#endif
