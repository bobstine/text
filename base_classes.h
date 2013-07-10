#ifndef _BASE_CLASSES_H_
#define _BASE_CLASSES_H_

#include <string>
#include <iostream>

using std::string;

class POS: public string
{  
 public:
  
  explicit POS()         : string("") { }
  explicit POS(string s) : string(s) { }
  explicit POS(char *s)  : string(s) { }
  
  void  print_to_stream(std::ostream &os) const   { os << "|" << static_cast<string>(*this) << "|"; }
};

inline
std::ostream&
operator<<(std::ostream &os, POS const& p) { p.print_to_stream(os); return os; }


class Type: public string
{
 public:
  explicit Type()         : string("") { }
  explicit Type(string s) : string(s)       { }
  explicit Type(char *s)  : string(s)       { }

  void  print_to_stream(std::ostream &os) const   { os << "`" << static_cast<string>(*this) << "'" ; }
};


inline
std::ostream&
operator<<(std::ostream &os, Type const& t) { t.print_to_stream(os); return os; }


#endif
