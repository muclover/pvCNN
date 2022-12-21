#pragma once
#include <string>
#include <map>
#include <iostream>
#include <fstream>

class Writer
{
public:
    Writer(std::string s)
    {
        out.open(s, std::ios::app);
    }
    ~Writer()
    {
        out.close();
    }
    void keep(std::string key, std::string value);
    void write(int N);

    void write_bit(std::string s, size_t bit)
    {
        out << s << ": " << bit << "\n";
    }

    void write_string(std::string s)
    {
        out << s << "\n";
    }

private:
    std::ofstream out;
    std::map<std::string, std::string> map;
};