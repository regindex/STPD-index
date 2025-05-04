#ifndef BITPACKED_TEXT_ORACLE_HPP
#define BITPACKED_TEXT_ORACLE_HPP

#include <common.hpp>

namespace stpd {

# define PACKED_INT_SIZE 2

class bitpacked_text_oracle
{
public:

    bitpacked_text_oracle(){};

    void build(std::string input_file_path)
    {
        std::ifstream file_text(input_file_path, std::ios::binary);
        file_text.seekg(0, std::ios::end);
        usafe_t N = file_text.tellg();
        file_text.seekg(0, std::ios::beg);

        this->T.resize(N);
        for(usafe_t i=0;i<N;++i)
        {
            char_t c;
            file_text.read(reinterpret_cast<char*>(&c), sizeof(char_t));
            if(dna_to_code_table[c] > 3)
                { std::cerr << "Non DNA character detected!" << std::endl; exit(1); }
            T[i] = dna_to_code_table[c];
        }

        //std::ofstream fout(input_file_path+".bitpacked", std::ios::binary);
        //store(fout);
        //fout.close();
    } 

    usafe_t total_length() const { return T.size(); }

    usafe_t serialize(std::ostream& out) const
    {
        //usafe_t w_bytes = (this->N * PACKED_INT_SIZE)/8;
        usafe_t w_bytes = T.serialize(out);

        return w_bytes;
    }

    void load(std::istream& in)
    {
        T.load(in);
    }

    void load(const std::string& input_file_path)
    {
        std::ifstream fin(input_file_path, std::ios::binary);
        load(fin);
        fin.close();
    }

    usafe_t size(){ return (T.size() * PACKED_INT_SIZE)/8; }

    unsigned char extract(usafe_t i)
    { 
        assert(i<T.size());
        return code_to_dna_table[this->T[i]]; 
    }
    
    usafe_t LCP(const std::string& pattern, usafe_t p, usafe_t t) const
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((pattern.size()-p),(this->T.size()-t));

        while(available_chars > 0)
        {
            if(pattern[p+matched_chars] != code_to_dna_table[this->T[t+matched_chars]])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }

    usafe_t LCS(const std::string& pattern, size_t p, size_t t) const
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != code_to_dna_table[this->T[t-matched_chars]])
                return matched_chars;

            matched_chars++;
            available_chars--;
        }

        return matched_chars;
    }

    std::pair<usafe_t,char_t> LCS_char(const std::string& pattern, usafe_t p, usafe_t t) const
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);

        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != code_to_dna_table[this->T[t-matched_chars]])
                return std::make_pair(matched_chars,code_to_dna_table[this->T[t-matched_chars]]);

            matched_chars++;
            available_chars--;
        }

        return std::make_pair(matched_chars,-1);
    }
    
private:

   sdsl::int_vector<PACKED_INT_SIZE> T;  
   // usafe_t N;
};

}  

#endif  // BITPACKED_TEXT_ORACLE_HPP