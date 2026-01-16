#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "avl.h"

using namespace std;

typedef std::pair<int64_t,int64_t> pair_t;

template <class POS_T, class LEN_T>
void break_long_blocks( AVLTree<POS_T>& T_P, std::map<POS_T,LEN_T>& B, LEN_T TH ) { 
    // break blocks longer than TH
    std::vector< std::pair<POS_T,LEN_T> > breakpoints;
    int64_t n = T_P.size( T_P.root() )-1;
    POS_T u = T_P.key( T_P.select(n-1) );
    for( int64_t i = 0; i < n; ++i ) {
        POS_T p = T_P.key(T_P.select(i));
        LEN_T l = T_P.key(T_P.select(i+1))-p;
        LEN_T b = B.at(p);

        LEN_T cut = TH;
        for( int64_t k = cut; k < l; k += cut ) {
            breakpoints.push_back( std::make_pair( p+k, b+k ) );
        }
    }
    for( int64_t i = 0; i < breakpoints.size(); ++i ) {
        POS_T p = breakpoints[i].first;
        LEN_T q = breakpoints[i].second;
        T_P.insert( p );
        B[p]=q;
    }
}

template <class POS_T, class LEN_T>
void break_every_d_blocks( AVLTree<POS_T>& T_P, std::map<POS_T,LEN_T>& B, LEN_T d ) {
    // if an image contains too many blocks (>d), we break it
    while( 1 ) {
        std::vector< std::pair<POS_T, LEN_T> > breakpoints;
        long long int n = T_P.size( T_P.root() )-1;
        for( long long int i = 0; i < n; ++i ) {
            POS_T p = T_P.key(T_P.select(i));
            LEN_T l = T_P.key(T_P.select(i+1))-p;
            POS_T b = B.at(p);

            int64_t blk_i = T_P.rank( b+1 );
            int64_t blk_j = T_P.rank( b+l );

            // if we want to have a threshold for breaking, uncomment the following line.
            //if( blk_i + 2*d < blk_j ) continue;

            // break the blocks
            for( int64_t k = blk_i + d - 1; k < blk_j; k += d ) {
                POS_T qk = T_P.key( T_P.select(k) );
                breakpoints.push_back( make_pair( qk-b+p, qk ) );
            }
        }

        if( breakpoints.empty() ) break;

        for( int64_t i = 0; i < breakpoints.size(); ++i ) {
            POS_T p = breakpoints[i].first;
            LEN_T q = breakpoints[i].second;
            T_P.insert( p );
            B[p]=q;
        }
    }
}

void break_phi_blocks(std::vector<pair_t>& blocks, int64_t M = 8)
{
    AVLTree<int64_t> T_P;
    map<int64_t,int64_t> B;

    {
        int64_t i = 0;
        vector<int64_t> P;
        for(const auto& b : blocks)
        {
            P.push_back( i );
            B[i] = b.first;
            i += b.second;
        }

        P.push_back( i );
        B[i]=i;

        // T_P can be initialized with a vector containing the starting positions of blocks, plus the size of the universe.
        T_P.init( P );
    }

    // calling the following two functions will break the blocks accordingly.
    //break_long_blocks( T_P, B, 4 ); // 4096, 8192, 16384, ... (instead of 4) for bigger data:
    break_every_d_blocks( T_P, B, M ); // 8 (instead of 2) for bigger data
   
    {
        int64_t n = T_P.size( T_P.root() )-1;
        blocks.resize(n);
        for( int64_t i = 0; i < n; ++i ) {
            int64_t p = T_P.key(T_P.select(i));
            int64_t l = T_P.key(T_P.select(i+1))-p;
            int64_t b = B.at(p);

            blocks[i] = make_pair(b,l);
        }
    }
}