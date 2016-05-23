//#include <boost/algorithm/string/predicate.hpp>
//#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>

using namespace std;
int main(){
  bool header = true;
  for (std::string line; std::getline(std::cin, line);){
      if(!line.empty()){
		  if (line.substr(0, 2)=="##"){
		if (header)
        cout << line <<endl;
        continue;
      } else if (line.substr(0,1) == "#"){
		 if (header){
          cout << line << endl;
          header = false;
		 }
          continue;
	  }
      cout << line << endl;
	  }
  }
}
