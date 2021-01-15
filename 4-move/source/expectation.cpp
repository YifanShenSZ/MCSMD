#include "../include/expectation.hpp"

ExpectationSet::ExpectationSet() {}
ExpectationSet::ExpectationSet(const std::vector<at::Tensor> & _expectations) : expectations_(_expectations) {}
ExpectationSet::~ExpectationSet() {}

at::Tensor ExpectationSet::operator[](const size_t & index) const {
    return expectations_[index];
}
ExpectationSet ExpectationSet::operator+(const ExpectationSet & ES2) const {
    assert(("Expectation sets 1 and 2 must have same size", expectations_.size() == ES2.expectations_.size()));
    std::vector<at::Tensor> news(expectations_.size());
    for (size_t i = 0; i < expectations_.size(); i++)
    if (expectations_[i].defined() && ES2.expectations_[i].defined()) {
        news[i] = expectations_[i].clone();
        news[i] += ES2.expectations_[i];
    }
    return ExpectationSet(news);
}
ExpectationSet ExpectationSet::operator*(const double & coeff) const {
    std::vector<at::Tensor> news(expectations_.size());
    for (size_t i = 0; i < expectations_.size(); i++)
    if (expectations_[i].defined()) {
        news[i] = expectations_[i].clone();
        news[i] *= coeff;
    }
    return ExpectationSet(news);
}

void ExpectationSet::add_(const ExpectationSet & ES2) {
    assert(("Expectation sets 1 and 2 must have same size", expectations_.size() == ES2.expectations_.size()));
    for (size_t i = 0; i < expectations_.size(); i++)
    if (expectations_[i].defined() && ES2.expectations_[i].defined())
    expectations_[i] += ES2.expectations_[i];
}

WindowedSet::WindowedSet() {}
WindowedSet::WindowedSet(const std::vector<at::Tensor> & _expectations, const tchem::gaussian::Gaussian & _window)
: ExpectationSet(_expectations), window_(_window) {}
WindowedSet::~WindowedSet() {}

WindowedSet WindowedSet::operator+(const ExpectationSet & ES2) const {
    assert(("Expectation sets 1 and 2 must have same size", size() == ES2.size()));
    std::vector<at::Tensor> news(size());
    for (size_t i = 0; i < size(); i++)
    if ((*this)[i].defined() && ES2[i].defined()) {
        news[i] = (*this)[i].clone();
        news[i] += ES2[i];
    }
    return WindowedSet(news, window_);
}
WindowedSet WindowedSet::operator*(const double & coeff) const {
    std::vector<at::Tensor> news(size());
    for (size_t i = 0; i < size(); i++)
    if ((*this)[i].defined()) {
        news[i] = (*this)[i].clone();
        news[i] *= coeff;
    }
    return WindowedSet(news, window_);
}
