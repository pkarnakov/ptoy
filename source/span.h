#pragma once

#include <cstddef>
#include <type_traits>

template <class T>
class span {
 public:
  // Types.
  using value_type = std::remove_cv_t<T>;
  using iterator = T*;
  using const_iterator = const T*;

  // Constructors, copy, and assignment.
  span() noexcept = default;
  template <class It>
  span(It first, size_t size) : data_(&(*first)), size_(size) {}
  template <class It, class End>
  span(It first, End last) : data_(&(*first)), size_(last - first) {}
  template <class Container>
  explicit span(Container& other) : span(other.begin(), other.end()) {}
  template <class OtherT>
  explicit span(const span<OtherT>& other) : span(other.data_, other.size_) {}
  ~span() noexcept = default;
  span& operator=(const span&) noexcept = default;

  // Observers.
  size_t size() const noexcept {
    return size_;
  }
  size_t size_bytes() const noexcept {
    return size_ * sizeof(value_type);
  }
  bool empty() const noexcept {
    return size_;
  }

  // Element access.
  T& operator[](size_t i) const {
    return data_[i];
  }
  T& front() const {
    return data_[0];
  }
  T& back() const {
    return data_[size_ - size_t(1)];
  }
  T* data() const noexcept {
    return data_;
  }

  // Iterator support.
  iterator begin() const noexcept {
    return data_;
  }
  iterator end() const noexcept {
    return data_ + size_;
  }
  const_iterator cbegin() const noexcept {
    return data_;
  }
  const_iterator cend() const noexcept {
    return data_ + size_;
  }

 private:
  T* data_ = nullptr;
  size_t size_ = 0;
};
