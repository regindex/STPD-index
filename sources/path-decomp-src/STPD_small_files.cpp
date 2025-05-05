#include <iostream>
#include <vector>
#include <string>
#include <sdsl/construct.hpp>
#include <cassert>
#include <set>
#include <stack>
#include <algorithm>

using namespace std;
using namespace sdsl;

size_t N;
int_vector< 8> T;
int_vector<32> SA;
int_vector<32> PA;
int_vector<32> Rank;
int_vector<32> ISA;
int_vector<32> LCP;
vector< pair<int,int> > ST;  // each node is a pair of parentID and depth
vector< int > STLeaf;

int_vector<32> suffix_lex( const int_vector< 8>& T, const int_vector<32>& ISA ) {
    return ISA;
}

int_vector<32> suffix_lex_r( const int_vector< 8>& T, const int_vector<32>& ISA ) {
    int_vector<32> ret;
    ret.resize( ISA.size() );
    for( size_t i = 0; i < ISA.size(); ++i ) {
        ret[i] = ISA.size()-1-ISA[i];
    }
    return ret;
}

int_vector<32> text_position( const int_vector< 8>& T ) {
    int_vector<32> ret;
    size_t n = T.size();
    ret.resize( n );
    for( size_t i = 0; i < n; ++i ) {
        ret[i] = i;
    }
    return ret;
}

int_vector<32> prefix_colex( const int_vector< 8>& T ) {
    int_vector<8> T_rev;
    size_t n = T.size();
    T_rev.resize( n );

    assert( T[n-1] == '\0' );
    T_rev[n-1] = '\0';
    for( size_t i = 0; i < n-1; ++i ) {
        T_rev[ n-2 - i ] = T[ i ];
    }
    
    int_vector<32> pa;
    algorithm::calculate_sa<>( (const unsigned char*) T_rev.data(), n, PA );

    int_vector<32> ret;
    ret.resize( n );
    for( size_t i = 1; i < n; ++i ) {
        ret[n-2 - PA[i]] = i;
    }
    ret[n-1] = 0;
    
    return ret;
}

int_vector<32> prefix_colex_r( const int_vector< 8>& T ) {
    int_vector<8> T_rev;
    size_t n = T.size();
    T_rev.resize( n );

    assert( T[n-1] == '\0' );
    T_rev[n-1] = '\0';
    for( size_t i = 0; i < n-1; ++i ) {
        T_rev[ n-2 - i ] = T[ i ];
    }
    
    algorithm::calculate_sa<>( (const unsigned char*) T_rev.data(), n, PA );

    int_vector<32> ret;
    ret.resize( n );
    for( size_t i = 1; i < n; ++i ) {
        ret[n-2 - PA[i]] = n-1-i;
    }
    ret[n-1] = n-1;
    
    return ret;
}

set<int> sampling( const int_vector<32>& rank ) {
    Rank = rank;
    vector< pair<int,int> > SortedV;

    for( int i = 0; i < N; ++i ) {
        SortedV.push_back( make_pair( Rank[i], i ) ); 
    }
    sort( SortedV.begin(), SortedV.end() );

    vector<bool> Marked( ST.size(), false );
    
    set<int> S;
    int Init = -1;
    for( int i = 0; i < N; ++i ) {      // rank order
        int j = ISA[SortedV[i].second]; // suffix order 
        int v = STLeaf[j]; 

        while( v >= 0 && !Marked[v] ) {
            Marked[v] = true;
            v = ST[v].first;
        }
        if( v >= 0 ) {
            S.insert( SA[j]+ST[v].second );
        } else {
            Init = SA[j];
            S.insert(SA[j]);
        }
    }

    return S;
}

void store_set_colex(const set<int>& sampling, const string output_file)
{
    // compute inverse prefix array
    int_vector<32> IPA; IPA.resize( N );
    for( size_t i = 0; i < N; ++i ) { IPA[ PA[i] ] = i; }

    vector<pair<int,int>> stpd_array; stpd_array.resize(sampling.size());
    int i=0;
    for(const uint64_t& m:sampling)
    { 
        if(m < N-1)
            stpd_array[i++] = make_pair(m,IPA[N-m-2]);
    }

    sort(stpd_array.begin(), stpd_array.end(), [](auto& left, auto& right)
        { return left.second < right.second; });

    ofstream output(output_file,ofstream::binary);

    for(size_t i=1;i<stpd_array.size();++i)
    { 
        uint64_t x = stpd_array[i].first;
        output.write((char*)&x,5);
    }

    output.close();
}

void store_set_lex(const set<int>& sampling, const string output_file)
{
    // compute prefix array and inverse prefix array
    {
        int_vector<8> T_rev;
        size_t n = T.size();
        T_rev.resize( n );

        assert( T[n-1] == '\0' );
        T_rev[n-1] = '\0';
        for( size_t i = 0; i < n-1; ++i ) {
            T_rev[ n-2 - i ] = T[ i ];
        }
        
        algorithm::calculate_sa<>( (const unsigned char*) T_rev.data(), n, PA );
    }

    // compute inverse prefix array
    int_vector<32> IPA; IPA.resize( N );
    for( size_t i = 0; i < N; ++i ) { IPA[ PA[i] ] = i; }

    vector<pair<int,int>> stpd_array; stpd_array.resize(sampling.size());
    int i=0;
    for(const uint64_t& m:sampling)
    { 
        if(m < N-1)
            stpd_array[i++] = make_pair(m,IPA[N-m-2]);
    }

    sort(stpd_array.begin(), stpd_array.end(), [](auto& left, auto& right)
        { return left.second < right.second; });

    ofstream output(output_file,ofstream::binary);

    for(size_t i=1;i<stpd_array.size();++i)
    { 
        uint64_t x = stpd_array[i].first;
        output.write((char*)&x,5);
    }

    output.close();
}

void output_Prefix_Array_BWT(const string output_file, const string output_file_BWT)
{
    // compute prefix array and inverse prefix array
    if(PA.size() == 0)
    {
        int_vector<8> T_rev;
        size_t n = T.size();
        T_rev.resize( n );

        assert( T[n-1] == '\0' );
        T_rev[n-1] = '\0';
        for( size_t i = 0; i < n-1; ++i ) {
            T_rev[ n-2 - i ] = T[ i ];
        }
        
        algorithm::calculate_sa<>( (const unsigned char*) T_rev.data(), n, PA );
    }

    ofstream output_pa(output_file,ofstream::binary);
    ofstream output_bwt(output_file_BWT,ofstream::binary);

    for(size_t i=0;i<PA.size();++i)
    { 
        uint64_t x = T.size() - PA[i] - 1;
        char c = T[x];
        output_pa.write((char*)&x,5);
        output_bwt.write((char*)&c,1);
    }

    output_pa.close();
    output_bwt.close();
}

void help()
{
    cout << "stpd [options]" << endl <<
    "Input: One text. Output: ST path decomposition of the text (stored using 5 bytes per element)." << endl <<
    "Options:" << endl <<
    "-h          Print this help" << endl <<
    "-i <arg>    Input Text (REQUIRED)" << endl <<
    "-o <arg>    Output file name (REQUIRED)" << endl <<
    "-c          Compute ST colex- sampling (DEFAULT)" << endl <<
    "-C          Compute ST colex+- sampling" << endl <<
    "-l          Compute ST lex- sampling" << endl <<
    "-L          Compute ST lex+- sampling" << endl <<
    "-P          Output the Prefix Array and the BWT of the reversed text" << endl;
    exit(0);
}

int main( int argc, char **argv ) {

    if(argc < 3) help();

    string input_filename, output_file;
    bool colexM = false, colexP = false, lexM = false, lexP = false;
    bool outPA_BWT = false;

    int opt;
    while ((opt = getopt(argc, argv, "hi:o:cClLP")) != -1){
        switch (opt){
            case 'h':
                help();
            break;
            case 'i':
                input_filename = string(optarg);
            break;
            case 'o':
                output_file = string(optarg);
            break;
            case 'c':
                colexM = true;
            break;
            case 'C':
                colexP = true;
            break;
            case 'l':
                lexM = true;
            break;
            case 'L':
                lexP = true;
            break;
            case 'P':
                outPA_BWT = true;
            break;
            default:
                help();
            return -1;
        }
    }

    { // load T
        load_vector_from_file( T, input_filename, 1 );
        append_zero_symbol( T );
        N = T.size();
    }

    // compute SA
    algorithm::calculate_sa<>( (const unsigned char*) T.data(), N, SA );

    { // compute LCP
        ISA.resize( N );
        for( size_t i = 0; i < N; ++i ) {
            ISA[ SA[i] ] = i;
        }

        LCP.resize( N );
        size_t m = 0;
        for( size_t i = 0; i < N-1; ++i ) {
            size_t j = SA[ISA[i]-1];
            while( m < N ) {
                if( T[i+m] != T[j+m] ) break;
                ++m;
            }
            LCP[ISA[i]] = m;
            if( m > 0 ) --m;
        }

        /** check LCP
        for( size_t i = 0; i < N-1; ++i ) {
            size_t m = 0;
            while( m < N ) {
                if( T[SA[i]+m] != T[SA[i+1]+m] ) break;
                ++m;
            }
            cout << LCP[i+1] << ' ' << m << endl;
        }
        **/
    }

    { // compute ST
        ST.push_back( make_pair( -1, 0 ) );
        ST.push_back( make_pair(  0, 1 ) );
        STLeaf.push_back( 1 );

        stack<int> stk;
        stk.push(0);
        stk.push(1);
        for( int i = 1; i < N; ++i ) {
            int last_u = stk.top();
            while( !stk.empty() && ST[stk.top()].second > LCP[i] ) {
                last_u = stk.top();
                stk.pop();
            }
            if( ST[stk.top()].second == LCP[i] ) {
                int new_v = ST.size();
                ST.push_back( make_pair( stk.top(), N-SA[i] ) );
                stk.push(new_v);
                STLeaf.push_back(new_v);
            } else {
                int new_u = ST.size();
                ST.push_back( make_pair( stk.top(), LCP[i] ) );
                ST[last_u] = make_pair( new_u, ST[last_u].second );
                stk.push(new_u);
                int new_v = ST.size();
                ST.push_back( make_pair( stk.top(), N-SA[i] ) );
                stk.push(new_v);
                STLeaf.push_back(new_v);
            }
        }
        // free LCP array
        LCP.resize(0);
    }

    { // ST colex sampling
        if(colexM or colexP)
        {
            set<int> col_set_m, col_set_p;

            col_set_m = sampling( prefix_colex( T ) );
            if(colexP)
            { 
                col_set_p = sampling( prefix_colex_r( T ) );
                col_set_m.insert( col_set_p.begin(), col_set_p.end() );
            }
            // free suffix tree
            STLeaf.clear();
            ST.clear();
            // free rank vector
            Rank.resize(0);
            // free sa/isa vectors
            SA.resize(0);
            ISA.resize(0);

            store_set_colex(col_set_m,output_file);
        }
    }

    { // ST lex sampling
        if(lexM or lexP)
        {
            set<int> lex_set_m, lex_set_p;

            lex_set_m = sampling( suffix_lex( T, ISA ) );
            if(lexP)
            { 
                lex_set_p = sampling( suffix_lex_r( T, ISA ) );
                lex_set_m.insert( lex_set_p.begin(), lex_set_p.end() );
            }
            // free suffix tree
            STLeaf.clear();
            ST.clear();
            // free rank vector
            Rank.resize(0);
            // free sa/isa vectors
            SA.resize(0);
            ISA.resize(0);

            store_set_lex(lex_set_m,output_file);
        }
    }

    if(outPA_BWT){ output_Prefix_Array_BWT(input_filename+".pa",input_filename+".rbwt"); }

    return 0;
}
