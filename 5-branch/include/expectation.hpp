#include <torch/torch.h>

#include <tchem/gaussian.hpp>

class ExpectationSet {
    private:
        std::vector<at::Tensor> expectations_;
    public:
        ExpectationSet();
        ExpectationSet(const std::vector<at::Tensor> & _expectations);
        ~ExpectationSet();

        inline std::vector<at::Tensor> expectations() const {return expectations_;}
        inline size_t size() const {return expectations_.size();}

        at::Tensor operator[](const size_t & index) const;
        ExpectationSet operator+(const ExpectationSet & ES2) const;
        ExpectationSet operator*(const double & coeff) const;

        void add_(const ExpectationSet & ES2);
};

class WindowedSet : public ExpectationSet {
    private:
        tchem::gaussian::Gaussian window_;
    public:
        WindowedSet();
        WindowedSet(const std::vector<at::Tensor> & _expectations, const tchem::gaussian::Gaussian & _window);
        ~WindowedSet();

        inline tchem::gaussian::Gaussian window() const {return window_;}

        WindowedSet operator+(const ExpectationSet & ES2) const;
        WindowedSet operator*(const double & coeff) const;
};
