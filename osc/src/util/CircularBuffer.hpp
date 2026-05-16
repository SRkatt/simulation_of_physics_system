#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>

template<typename T>
class CircularBuffer {
    std::vector<T> data_;
    std::size_t cap_;
    std::size_t head_ = 0;
    std::size_t count_ = 0;
public:
    explicit CircularBuffer(std::size_t cap) : data_(cap), cap_(cap) {}
    void push(T v) { data_[head_] = v; head_ = (head_ + 1) % cap_; if (count_ < cap_) ++count_; }
    std::size_t size() const { return count_; }
    bool empty() const { return count_ == 0; }
    T operator[](std::size_t i) const {
        if (i >= count_) throw std::out_of_range("CircularBuffer index");
        std::size_t idx = (head_ + cap_ - count_ + i) % cap_;
        return data_[idx];
    }
    T back() const { return (*this)[count_ - 1]; }
    void clear() { count_ = 0; head_ = 0; }
};
