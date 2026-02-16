// Copyright (c) 2025, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#pragma once

#include <string>
#include <vector>
#include <stack>

namespace stpd{

struct Node;
Node* prev_leaf = nullptr;

struct Node{

	Node(uint_t b, uint_t e) : text_begin(b), text_end(e) {}
	Node(uint_t b, uint_t e, uint_t d) : text_begin(b), text_end(e), depth(d) {}

	usafe_t serialize(std::ostream& out)
	{
		usafe_t tot_bytes = 0;

		out.write((char*)&text_begin, sizeof(text_begin));
		out.write((char*)&text_end, sizeof(text_end));
		out.write((char*)&depth, sizeof(depth));

		for(auto& child : children)
		{
			bool empty_child = child == nullptr ? true : false;
			out.write((char*)&empty_child, sizeof(empty_child));
			if(not empty_child)
				tot_bytes += child->serialize(out);
		}

		tot_bytes += sizeof(text_begin) + sizeof(text_end) + sizeof(depth) + sizeof(next_leaf)*4;
		return tot_bytes;
	}

	void load(std::istream& in)
	{
		in.read((char*)&text_begin, sizeof(text_begin));
		in.read((char*)&text_end, sizeof(text_end));
		in.read((char*)&depth, sizeof(depth));

		bool is_leaf = (text_begin == 0) and (text_end > 0) ? true : false;
		if(is_leaf)
		{
			if(prev_leaf == nullptr){ prev_leaf = this; }
			else
			{
				prev_leaf->next_leaf = this;
				prev_leaf = this;
			}
		}

		bool empty_child;
		for(auto& child : children)
		{
			in.read((char*)&empty_child, sizeof(empty_child));
			if(not empty_child)
			{
				child = new Node(0,0,0);
				child->load(in);
			}
		}
	}

	Node*  children[4] = {nullptr, nullptr, nullptr, nullptr};
	Node*  next_leaf = nullptr;
	uint_t text_begin = 0;
	uint_t text_end = 0;
	uint_t depth = 0;
};

template<class text_oracle = RLZ_DNA_sux<>>
class PrefixTree{

private:
	text_oracle* text;
	Node* root;

public:
	PrefixTree(){};

	void build(text_oracle* text_,
		       const std::string& pa_file,const std::string& lcs_file,
		       bool verbose = true)
	{
		// set pointer to the text oracle
		text = text_;
		// open files
		std::ifstream PA (pa_file,  std::ios::binary);
		if (!PA.is_open()) { std::cerr << "Error1: Could not open file '" << pa_file << std::endl; }
		std::ifstream LCS(lcs_file, std::ios::binary);
		if (!LCS.is_open()) { std::cerr << "Error2: Could not open file '" << lcs_file << std::endl; }
		// skip dollar character
	    PA.seekg(5, std::ios::beg);
	    LCS.seekg(5, std::ios::beg);
		// set $ leaf
		root = new Node(0,0,0);
		root->next_leaf = new Node(0,0,1);
		Node* prev_leaf = root->next_leaf;
		char c;
		uint64_t pa = 0, lcs = 0;

		// set all other leaves
		std::stack<Node*> st; // stack for pt construction
		st.push(root);

		size_t no_nodes = 1;

		for (usafe_t i = 1; i < text->text_length()+1; ++i)
		{
			PA.read(reinterpret_cast<char*>(&pa), 5);
			LCS.read(reinterpret_cast<char*>(&lcs), 5);

			// go until you get to the right node
	        while (st.top()->depth > lcs){ st.pop(); }

	        auto parent = st.top();

	        // create a new node if necesarry
	        if (parent->depth < lcs)
	        {
	        	Node* internal = new Node(pa-lcs,pa-parent->depth,lcs);
	        	no_nodes++;
	        	c = text->extract(pa-parent->depth-1);

	        	Node* tmp = parent->children[dna_to_code_table[c]];
	        	tmp->text_end -= (lcs - parent->depth);

	        	parent->children[dna_to_code_table[c]] = internal;

	        	c = text->extract(tmp->text_end-1);
	        	internal->children[dna_to_code_table[c]] = tmp;

	        	st.push(internal);
	        	parent = internal;
	        }

	        // create new leaf
	        Node* leaf = new Node(0, pa-lcs, pa);
	        c = text->extract(pa-lcs-1);
	        parent->children[dna_to_code_table[c]] = leaf;
	        st.push(leaf);

	        prev_leaf->next_leaf = leaf;
	        prev_leaf = leaf;
	        no_nodes++;
		}

		PA.close();
		LCS.close();

		return;
	}

	Node* locate_locus(const std::string& pattern) const
	{
		usafe_t i = pattern.size();
		Node* curr = root;
		while(i > 0)
		{
			curr = curr->children[dna_to_code_table[pattern[i-1]]];
			if(curr == nullptr){ return curr; }
			i--;

			usafe_t edge_len = curr->text_end - curr->text_begin;
			if(edge_len > 1)
			{
				usafe_t to_match = std::min(i,edge_len-1);
				auto j = text->LCS_char(pattern,i-1,curr->text_end-2); 
				usafe_t mlen = std::min(static_cast<usafe_t>(j.first),edge_len-1);

				if(mlen < to_match){ return nullptr; }

				i -= mlen;
			}
		}

		return curr;
	}

	Node* locate_smallest_leaf(Node* locus) const
	{
		Node* curr = locus;
		while(curr->text_begin != 0)
		{
			for(usafe_t i=0;i<4;++i)
			{
				Node* tmp = curr->children[i];
				if(tmp != nullptr){ curr = tmp; break; }
			}
		}

		return curr;
	}

	Node* locate_largest_leaf(Node* locus) const
	{
		Node* curr = locus;
		bool is_leaf = std::all_of(std::begin(curr->children), std::end(curr->children), [](Node* ptr) {
    					return ptr == nullptr; });
		while(not is_leaf)
		{
			for(usafe_t i=4;i>0;--i)
			{
				Node* tmp = curr->children[i-1];
				if(tmp != nullptr){ curr = tmp; break; }
			}
			is_leaf = std::all_of(std::begin(curr->children), std::end(curr->children), [](Node* ptr) {
			    					return ptr == nullptr; });
		}

		return curr;
	}

	std::vector<uint_t> get_SA_range(Node* b, Node* e) const
	{
		std::vector<uint_t> res;

		Node* curr = b;
		uint_t p = curr->depth - 1;
		uint_t p_e = e->depth - 1;
		res.push_back(p);

		while(p != p_e)
		{
			curr = curr->next_leaf;
			p = curr->depth - 1;

			res.push_back(p);
		}

		return res;
	}

	std::vector<uint_t>
		get_SA_range_thr(Node* b, Node* e, usafe_t thr) const
	{
		std::vector<uint_t> res { };

		Node* curr = b;
		res.push_back(curr->depth - 1);

		while(--thr > 0 and curr != e)
		{
			curr = curr->next_leaf;
			res.push_back(curr->depth - 1);
		}

		return res;
	}

	uint_t operator[](Node* i) const
	{
		return i->depth - 1;
	}

	void print_tree(Node* node, int level = 0) const
	{
	    if (node != nullptr)
	    {
	    	std::cout << std::string(level * 2, ' ') << node->text_begin << "-" << node->text_end << std::endl;
	    	std::string t(node->text_end - node->text_begin,0);
	    	for(usafe_t i=0;i<(node->text_end - node->text_begin);++i)
	    	{
	    		t.push_back(text->extract(node->text_end-i-1));
	    	}
	        std::cout << std::string(level * 2, ' ') << t << "\n";

		    for (Node* child : node->children)
		    {
		        print_tree(child, level + 1);
		    }
	    }
	}

	usafe_t serialize(std::ostream& out)
	{
		usafe_t w_bytes = 0;

		w_bytes += root->serialize(out);

		return w_bytes;
	}

	void load(std::istream& in, text_oracle* text_)
	{
		text = text_;

		root = new Node(0,0,0);
		root->load(in);
	}
};
}