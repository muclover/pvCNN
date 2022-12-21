#include "writer.hpp"

void Writer::keep(std::string key, std::string value)
{
    map[key] = value;
}

void Writer::write(int N)
{
    out << "\n"
        << __DATE__ << " " << __TIME__ << "************* N = "
        << N << "********************\n";
    out << "------------ Time ------------------------------\n";
    auto iter = map.begin();
    while (iter != map.end())
    {
        out << iter->first << "\t\t:\t" << iter->second << "  seconds"
            << "\n";
        iter++;
    }
    out << "\n";
}
