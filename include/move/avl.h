// modified from https://github.com/sunghwank/spmindex/blob/main/avl.h
#ifndef __MY_AVL_H__
#define __MY_AVL_H__

#include <vector>
#include <cassert>

template <class K>
class AVLTree {
    public:
    typedef K   key_t;
    typedef int64_t node_t;

    AVLTree() {

    }

    void init( int64_t init_buf_size ) {
        m_root = -1;
        m_key       .resize( init_buf_size );
        m_left      .resize( init_buf_size );
        m_right     .resize( init_buf_size );
        m_height    .resize( init_buf_size );
        m_size      .resize( init_buf_size );
        m_free_nodes.resize( init_buf_size );
        for( int64_t i = 0; i < init_buf_size; ++i ) {
            m_free_nodes[i] = init_buf_size-i-1;
        }
    }

    void init( const std::vector<K>& keys, float buf_factor = 1.2 ) {
        init( (1+keys.size())*buf_factor );
        m_root = m_init( keys, 0, keys.size() );
    }

    node_t m_init( const std::vector<K>& keys, size_t i, size_t j ) {
        if( i == j ) return -1;
        size_t m = (i+j)/2;
        node_t v_m = m_new_node( keys[m] );
        node_t l = m_init( keys, i  , m );
        node_t r = m_init( keys, m+1, j );
        size_t h_l = height(l);
        size_t h_r = height(r);
        m_height[v_m] = (h_l>h_r?h_l:h_r)+1;
        m_size  [v_m] = size(l)+size(r)+1;
        m_left  [v_m] = l;
        m_right [v_m] = r;
        return v_m;
    }
 
    inline node_t root  ( void     ) const { return m_root;    }
    inline size_t size  ( node_t v ) const { if( v == -1 ) return 0; else return m_size  [v]; }
    inline size_t height( node_t v ) const { if( v == -1 ) return 0; else return m_height[v]; }
    inline node_t left  ( node_t v ) const { return m_left[v]; }
    inline node_t right ( node_t v ) const { return m_right[v];}
    inline key_t  key   ( node_t v ) const { if( v == -1 ) return K(); return m_key[v];}
    
    void insert( key_t k ) {
        if( m_root == -1 ) m_root = m_new_node( k );
        else m_root = m_insert( m_root, k );
    }

    void remove( key_t k ) {
        m_root = m_remove( m_root, k );
    }

    node_t lower_bound( K k ) const {
        node_t last_left = -1;
        node_t v = m_root;
        while( v != -1 ) {
            if( k < m_key[v] ) {
                v = left(v);
            } else if( k > m_key[v] ) {
                last_left = v;
                v = right(v);
            } else {
                return v;
            }
        }
        return last_left;
    }

    size_t rank( K k ) {
        size_t count = 0;
        node_t v = m_root;
        while( v != -1 ) {
            if( k < m_key[v] ) {
                v = left(v);
            } else if( k > m_key[v] ) {
                count += size(left(v))+1;
                v = right(v);
            } else {
                return count+size(left(v));
            }
        }
        return count;
    }

    node_t select( size_t s ) {
        size_t count = 0;
        node_t v = m_root;
        while( v != -1 ) {
            size_t sz = size(left(v));
            if( s < sz ) {
                v = left(v);
            } else if( s > sz ) {
                s -= sz+1;
                v = right(v);
            } else {
                return v;
            }
        }
        assert(0);
    }

    private:
    int64_t bf( node_t v ) {
        if( v == -1 ) return 0;
        int64_t h_l = height(left (v));
        int64_t h_r = height(right(v));
        return h_l - h_r;
    }

    void m_update( node_t v ) {
        m_size[v] = size(left(v))+size(right(v))+1;
        int64_t h_l = height(left(v));
        int64_t h_r = height(right(v));
        int64_t h = h_l>h_r?h_l:h_r;
        m_height[v] = h + 1;
    }

    node_t m_balance( node_t v ) {
        if( bf(v) > 1 ) { // left heavy
            if( bf(left(v)) < 0 ) {
                m_left[v] = rotate_l(left(v));
            }
            return rotate_r(v);
        } else if( bf(v) < -1 ) { // right heavy
            if( bf(right(v)) > 0 ) {
                m_right[v] = rotate_r(right(v));
            }
            return rotate_l(v);
        }
        return v;
    }

    node_t m_insert( node_t v, key_t k ) {
        if( v == -1 ) return m_new_node( k );
        if( k < m_key[v] ) {
            m_left[v] = m_insert( left(v), k );
        } else if( k > m_key[v] ) {
            m_right[v] = m_insert( right(v), k );
        } else {
            return v;
        }
        m_update(v);
        return m_balance(v);
    }

    node_t m_remove( node_t v, key_t k ) {
        if( v == -1 ) return v;

        if( k < m_key[v] ) {
            m_left[v] = m_remove( left(v), k );
        } else if( k > m_key[v] ) {
            m_right[v] = m_remove( right(v), k );
        } else {
            if( left(v) != -1 && right(v) != -1 ) {
                node_t u = right(v);
                while( left(u) != -1 ) {
                    u = left(u);
                }
                int64_t k_u = m_key[u];
                m_right[v] = m_remove( right(v), k_u );
                m_key  [v] = k_u;
            } else {
                m_free_node(v);
                if( left(v) == -1 ) {
                    return right(v);
                } else if( right(v) == -1 ) {
                    return left(v);
                }
            }
        }
        m_update(v);
        return m_balance(v);
    }

    void m_new_free_nodes() {
        size_t new_node_id = m_key.size();
        m_key.push_back(-1);
        size_t ubound_node_id = m_key.capacity();

        for( size_t v = new_node_id; v < ubound_node_id; ++v ) {
            m_free_nodes.push_back( ubound_node_id-1 - v + new_node_id );
        }
        m_key       .resize( ubound_node_id );
        m_left      .resize( ubound_node_id );
        m_right     .resize( ubound_node_id );
        m_height    .resize( ubound_node_id );
        m_size      .resize( ubound_node_id );
    }

    node_t m_new_node( key_t k ) { 
        size_t n = m_free_nodes.size();
        if( n == 0 ) {
            m_new_free_nodes();
            n = m_free_nodes.size();
        }
        node_t v = m_free_nodes[n-1];
        m_free_nodes.resize(n-1);
        m_key   [v] = k;
        m_left  [v] = m_right[v] = -1;
        m_height[v] = 1;
        m_size  [v] = 1;
        return v;
    }

    void m_free_node( node_t v ) {
        m_free_nodes.push_back( v ); 
    }

    node_t rotate_r( node_t v ) {
        node_t l  = left (v);
        node_t lr = right(l);
        m_right[l] = v;
        m_left [v] = lr;
        m_update(v);
        m_update(l);
        return l;
    }

    node_t rotate_l( node_t v ) {
        node_t r  = right(v);
        node_t rl = left (r);
        m_left [r] = v;
        m_right[v] = rl;
        m_update(v);
        m_update(r);
        return r;
    }

    node_t m_root;

    std::vector<node_t> m_key;

    std::vector<node_t> m_left;
    std::vector<node_t> m_right;
    std::vector<size_t> m_height;
    std::vector<size_t> m_size;

    std::vector<node_t> m_free_nodes;
};

#endif
