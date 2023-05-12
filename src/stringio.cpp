#include <vector>
#include <string>
#include <sstream>
#include <iterator>

#include "string.h"

std::vector<std::string> split(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret(
        (std::istream_iterator<std::string>(buffer)), 
        std::istream_iterator<std::string>()
    );
    return ret;
}