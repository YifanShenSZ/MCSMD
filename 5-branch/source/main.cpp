#include <torch/torch.h>

#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>
#include <tchem/linalg.hpp>
#include <tchem/gaussian.hpp>

#include "../include/basic.hpp"

void update_expectations(
std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window,
std::vector<at::Tensor> & windowed_expectations, const double & dt);

std::vector<at::Tensor>  update_windows(
const std::vector<at::Tensor> & expectations,
const tchem::gaussian::Gaussian & window);

argparse::ArgumentParser parse_args(const int & argc, const char ** & argv) {
    CL::utility::EchoCommand(argc, argv); std::cout << '\n';
    argparse::ArgumentParser parser("2nd order semiclassical Moyal dynamics");

    parser.add_argument("--miu_q" , 1, false);
    parser.add_argument("--miu_p" , 1, false);
    parser.add_argument("--sigma_q", 1, false);
    parser.add_argument("--sigma_p", 1, false);

    parser.add_argument("-t","--total_time", 1, true, "default = 10000");
    parser.add_argument("-d","--time_step", 1, true, "default = 1");
    parser.add_argument("-o","--output_interval", 1, true, "default = 10");

    parser.parse_args(argc, argv);
    return parser;
}

int main(int argc, const char** argv) {
    argparse::ArgumentParser args = parse_args(argc, argv);
    initialize_random();
    auto top = at::TensorOptions().dtype(torch::kFloat64);

    at::Tensor miu = at::empty(2, top),
               var = at::empty({2, 2}, top);
    miu[0] = args.retrieve<double>("miu_q");
    miu[1] = args.retrieve<double>("miu_p");
    var[0][0] = pow(args.retrieve<double>("sigma_q"), 2.0);
    var[0][1] = 0.0;
    var[1][1] = pow(args.retrieve<double>("sigma_p"), 2.0);

    std::vector<at::Tensor> expectations(3);
    expectations[0] = at::ones(1, top)[0];
    expectations[1] = miu.new_empty(miu.sizes());
    expectations[1].copy_(miu);
    expectations[2] = var.new_empty(var.sizes());
    expectations[2].copy_(var + miu.outer(miu));

    tchem::gaussian::Gaussian distribution(miu, var);
    tchem::gaussian::Gaussian window      (miu, var);
    at::Tensor coeff;
    tchem::gaussian::Gaussian product;
    std::tie(coeff, product) = distribution * window;
    tchem::polynomial::PolynomialSet variable_set(2 * dimension, 2);
    at::Tensor integrals = coeff * product.integral(variable_set);
    std::vector<at::Tensor> windowed_expectations = variable_set.views(integrals);
    windowed_expectations[0] = windowed_expectations[0][0];
    // 1 is a vector, so no need to transform to scalor nor tensor
    windowed_expectations[2] = tchem::LA::vec2sytensor(windowed_expectations[2], 2, 2 * dimension);

    double total_time = 10000.0;
    if (args.gotArgument("total_time")) total_time = args.retrieve<double>("total_time");
    double time_step = 1.0;
    if (args.gotArgument("time_step")) time_step = args.retrieve<double>("time_step");
    double output_interval = 10.0;
    if (args.gotArgument("output_interval")) output_interval = args.retrieve<double>("output_interval");

    std::ofstream ofs0; ofs0.open("SMD.txt");
    std::ofstream ofs1; ofs1.open("window.txt");
    std::ofstream ofs2; ofs2.open("prediction.txt");
    size_t NSnapshots = total_time / time_step;
    size_t OutputStep = output_interval / time_step;
    double time = 0.0;
    // Output initial expectations
    double miu_q = expectations[1][0].item<double>(),
           miu_p = expectations[1][1].item<double>();
    double sigmaq = sqrt(expectations[2][0][0].item<double>() - miu_q * miu_q),
           sigmap = sqrt(expectations[2][1][1].item<double>() - miu_p * miu_p);
    double rho = (expectations[2][0][1].item<double>() - miu_q * miu_p) / sigmaq / sigmap;
    ofs0 << "time    <q>    <p>    sigma_x    rho    sigma_p\n"
         << time << "    " << miu_q << "    " << miu_p << "    "
         << sigmaq << "    " << rho << "    " << sigmap << std::endl;
    // Output initial windowed expectations
    double population = windowed_expectations[0].item<double>();
    miu_q = windowed_expectations[1][0].item<double>() / population;
    miu_p = windowed_expectations[1][1].item<double>() / population;
    sigmaq = sqrt(windowed_expectations[2][0][0].item<double>() / population - miu_q * miu_q);
    sigmap = sqrt(windowed_expectations[2][1][1].item<double>() / population - miu_p * miu_p);
    rho = (windowed_expectations[2][0][1].item<double>() / population  - miu_q * miu_p) / sigmaq / sigmap;
    ofs1 << "time    population    <q>    <p>    sigma_x    rho    sigma_p\n"
         << time << "    " << population << "    " << miu_q << "    " << miu_p << "    "
         << sigmaq << "    " << rho << "    " << sigmap << std::endl;
    // Output initial predictions
    std::vector<at::Tensor> predictions = update_windows(expectations, window);
    population = predictions[0].item<double>();
    miu_q = predictions[1][0].item<double>() / population;
    miu_p = predictions[1][1].item<double>() / population;
    sigmaq = sqrt(predictions[2][0][0].item<double>() / population - miu_q * miu_q);
    sigmap = sqrt(predictions[2][1][1].item<double>() / population - miu_p * miu_p);
    rho = (predictions[2][0][1].item<double>() / population  - miu_q * miu_p) / sigmaq / sigmap;
    ofs2 << "time    population    <q>    <p>    sigma_x    rho    sigma_p\n"
         << time << "    " << population << "    " << miu_q << "    " << miu_p << "    "
         << sigmaq << "    " << rho << "    " << sigmap << std::endl;
    for (size_t i = 1; i <= NSnapshots; i++) {
        update_expectations(expectations, window, windowed_expectations, time_step);
        time += time_step;
        if (i % OutputStep == 0) {
            // Output expectations
            double miu_q = expectations[1][0].item<double>(),
                   miu_p = expectations[1][1].item<double>();
            double sigmaq = sqrt(expectations[2][0][0].item<double>() - miu_q * miu_q),
                   sigmap = sqrt(expectations[2][1][1].item<double>() - miu_p * miu_p);
            double rho = (expectations[2][0][1].item<double>() - miu_q * miu_p) / sigmaq / sigmap;
            ofs0 << time << "    " << miu_q << "    " << miu_p << "    "
                 << sigmaq << "    " << rho << "    " << sigmap << std::endl;
            // Output windowed expectations
            double population = windowed_expectations[0].item<double>();
            miu_q = windowed_expectations[1][0].item<double>() / population;
            miu_p = windowed_expectations[1][1].item<double>() / population;
            sigmaq = sqrt(windowed_expectations[2][0][0].item<double>() / population - miu_q * miu_q);
            sigmap = sqrt(windowed_expectations[2][1][1].item<double>() / population - miu_p * miu_p);
            rho = (windowed_expectations[2][0][1].item<double>() / population  - miu_q * miu_p) / sigmaq / sigmap;
            ofs1 << time << "    " << population << "    " << miu_q << "    " << miu_p << "    "
                 << sigmaq << "    " << rho << "    " << sigmap << std::endl;
            // Output predictions
            std::vector<at::Tensor> predictions = update_windows(expectations, window);
            population = predictions[0].item<double>();
            miu_q = predictions[1][0].item<double>() / population;
            miu_p = predictions[1][1].item<double>() / population;
            sigmaq = sqrt(predictions[2][0][0].item<double>() / population - miu_q * miu_q);
            sigmap = sqrt(predictions[2][1][1].item<double>() / population - miu_p * miu_p);
            rho = (predictions[2][0][1].item<double>() / population  - miu_q * miu_p) / sigmaq / sigmap;
            ofs2 << time << "    " << population << "    " << miu_q << "    " << miu_p << "    "
                 << sigmaq << "    " << rho << "    " << sigmap << std::endl;
        }
    }
    ofs0.close(); ofs1.close(); ofs2.close();
}