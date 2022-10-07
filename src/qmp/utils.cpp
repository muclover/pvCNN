#include "utils.hpp"

std::string string_to_hex(const std::string &in)
{
    static const char *const ch = "0123456789ABCDEF";
    size_t size = in.size();
    std::string output;
    output.reserve(2 * size);
    for (size_t i = 0; i < size; ++i)
    {
        const unsigned char c = in[i];
        output.push_back(ch[c >> 4]);
        output.push_back(ch[c & 15]);
    }
    return output;
}

std::string SHA256(std::string data)
{
    typedef unsigned char byte;
    const byte *pb = (byte *)data.data();
    unsigned int len = data.length();
    byte digest[CryptoPP::SHA256::DIGESTSIZE];

    CryptoPP::SHA256().CalculateDigest(digest, pb, len);

    return std::string((char *)digest, CryptoPP::SHA256::DIGESTSIZE);
}

void SHA256_Bytes(std::string data, unsigned char *output)
{
    typedef unsigned char byte;
    byte const *pb = (byte *)data.data();
    unsigned int len = data.length();
    byte abDigest[CryptoPP::SHA256::DIGESTSIZE];

    CryptoPP::SHA256().CalculateDigest(abDigest, pb, len);

    output = new unsigned char[CryptoPP::SHA256::DIGESTSIZE];
}

void test_read_dir(std::string PATH)
{
    struct dirent *ptr;
    DIR *dir;
    dir = opendir(PATH.c_str());
    std::vector<std::string> files;

    while ((ptr = readdir(dir)) != nullptr)
    {
        if (ptr->d_name[0] == '.')
            continue;
        files.push_back(ptr->d_name);
    }
    for (int i = 0; i < files.size(); ++i)
    {
        std::cout << files[i] << std::endl;
    }
    closedir(dir);
}