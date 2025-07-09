#ifndef BITPACKED_VECTOR_HPP_
#define BITPACKED_VECTOR_HPP_

#include <cstdint>
#include <string>

std::string uint128_to_string(__uint128_t value)
{
    if (value == 0) return "0";

    std::string result;
    while (value > 0) {
        result = "0123456789"[value % 10] + result;
        value /= 10;
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, __uint128_t value) {
    return os << uint128_to_string(value);
}

class bitpacked_vector {

public:
    bitpacked_vector() : width_(0), data(nullptr) {}

    void width(uint8_t width__) { width_ = width__; }

    void resize(size_t size)
    {
    	size_ = size;
        size_t no_bits = size_ * width_;
        size_t total_bytes = (no_bits + 7) / 8; 
        delete[] data; 
        data = new char[total_bytes](); 
    }

	void resize_aligned(size_t size)
	{
	    size_ = size;
	    byte_width = (width_ + 7) / 8;       
	    size_t total_bytes = size_ * byte_width;    

	    delete[] data;
	    data = new char[total_bytes]();             
	}

    ~bitpacked_vector()
    {
        delete[] data;
    }

	void set(size_t i, __uint128_t value) {
	    assert(i < size_);
	    assert(width_ > 0 && width_ <= 128);
	    assert(value < (__uint128_t(1) << width_));

	    size_t bit_pos = i * width_;
	    size_t byte_pos = bit_pos / 8;
	    size_t bit_offset = bit_pos % 8;

	    // Read enough bytes to cover the 128-bit value + possible offset
	    const size_t bits_needed = bit_offset + width_;
	    const size_t bytes_needed = (bits_needed + 7) / 8;

	    // Load current bytes into 128-bit block
	    __uint128_t block = 0;
	    for (size_t j = 0; j < bytes_needed; ++j) {
	        block |= (__uint128_t(static_cast<unsigned char>(data[byte_pos + j])) << (8 * j));
	    }

	    // Clear the target bit region
	    __uint128_t mask = ((__uint128_t(1) << width_) - 1);
	    block &= ~(mask << bit_offset);

	    // Set the new value
	    block |= (value & mask) << bit_offset;

	    // Store back
	    for (size_t j = 0; j < bytes_needed; ++j) {
	        data[byte_pos + j] = static_cast<char>((block >> (8 * j)) & 0xFF);
	    }
	}

	__uint128_t get(size_t i) const {
	    assert(i < size_);
	    assert(width_ > 0 && width_ <= 128);

	    size_t bit_pos = i * width_;
	    size_t byte_pos = bit_pos / 8;
	    size_t bit_offset = bit_pos % 8;

	    const size_t bits_needed = bit_offset + width_;
	    const size_t bytes_needed = (bits_needed + 7) / 8;

	    __uint128_t block = 0;
	    for (size_t j = 0; j < bytes_needed; ++j) {
	        block |= (__uint128_t(static_cast<unsigned char>(data[byte_pos + j])) << (8 * j));
	    }

	    block >>= bit_offset;
	    block &= (__uint128_t(1) << width_) - 1;
	    return block;
	}

	void set_aligned(size_t i, __uint128_t value)
	{
	    //size_t byte_width = (width_ + 7) / 8;
	    assert(i < size_);
	    assert(value < (__uint128_t(1) << width_));

	    size_t byte_pos = i * byte_width;
	    for (size_t j = 0; j < byte_width; ++j) {
	        data[byte_pos + j] = static_cast<char>((value >> (8 * j)) & 0xFF);
	    }
	}

	__uint128_t get_aligned(size_t i) const
	{
	    //size_t byte_width = (width_ + 7) / 8;
	    assert(i < size_);

	    size_t byte_pos = i * byte_width;
	    __uint128_t value = 0;

	    for (size_t j = 0; j < byte_width; ++j) {
	        value |= (static_cast<__uint128_t>(
	                     static_cast<unsigned char>(data[byte_pos + j])
	                 ) << (8 * j));
	    }

	    return value;
	}
	/*
    __uint128_t operator[](size_t i) const
    {
        return get(i);
    }
    */
	size_t serialize(std::ostream& out) const
	{
	    size_t total_bytes = size_ * byte_width;

	    out.write(reinterpret_cast<const char*>(&width_), sizeof(width_));
	    out.write(reinterpret_cast<const char*>(&size_), sizeof(size_));
	    out.write(reinterpret_cast<const char*>(&byte_width), sizeof(byte_width));
	    out.write(data, total_bytes);

	    return sizeof(width_) + sizeof(size_) + sizeof(byte_width) + total_bytes;
	}

	void load(std::istream& in)
	{
	    in.read(reinterpret_cast<char*>(&width_), sizeof(width_));
	    in.read(reinterpret_cast<char*>(&size_), sizeof(size_));
	    in.read(reinterpret_cast<char*>(&byte_width), sizeof(byte_width));

	    size_t total_bytes = size_ * byte_width;

	    delete[] data;
	    data = new char[total_bytes]();
	    in.read(data, total_bytes);
	}

private:
    char* data;
    uint8_t width_;
    uint8_t byte_width;
    size_t size_;
};

#endif 