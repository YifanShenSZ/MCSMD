#include <torch/torch.h>

#include <CppLibrary/argparse.hpp>
#include <CppLibrary/utility.hpp>
#include <tchem/polynomial.hpp>

#include "../include/basic.hpp"

void update_expectations(std::vector<at::Tensor> & expectations, const double & dt);

argparse::ArgumentParser parse_args(const int & argc, const char ** & argv) {
    CL::utility::EchoCommand(argc, argv); std::cout << '\n';
    argparse::ArgumentParser parser("2nd order semiclassical Moyal dynamics");

    parser.add_argument("-q", "--average_q" , 1, false);
    parser.add_argument("-p", "--average_p" , 1, false);
    parser.add_argument("--qq", 1, false);
    parser.add_argument("--pp", 1, false);

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

    std::vector<at::Tensor> expectations(3);
    expectations[1] = at::empty(2, top);
    expectations[1][0] = args.retrieve<double>("average_q");
    expectations[1][1] = args.retrieve<double>("average_p");
    expectations[2] = at::empty({2, 2}, top);
    expectations[2][0][0] = args.retrieve<double>("qq");
    expectations[2][0][1] = 0.0;
    expectations[2][1][1] = args.retrieve<double>("pp");

    double total_time = 10000.0;
    if (args.gotArgument("total_time")) total_time = args.retrieve<double>("total_time");
    double time_step = 1.0;
    if (args.gotArgument("time_step")) time_step = args.retrieve<double>("time_step");
    double output_interval = 10.0;
    if (args.gotArgument("output_interval")) output_interval = args.retrieve<double>("output_interval");

    std::ofstream ofs; ofs.open("SMD.txt");
    size_t NSnapshots = total_time / time_step;
    size_t OutputStep = output_interval / time_step;
    double time = 0.0;
    ofs << time << "    "
        << expectations[1][0].item<double>() << "    "
        << expectations[1][1].item<double>() << "    "
        << expectations[2][0][0].item<double>() << "    "
        << expectations[2][0][1].item<double>() << "    "
        << expectations[2][1][1].item<double>() << '\n';
    for (size_t i = 1; i <= NSnapshots; i++) {
        update_expectations(expectations, time_step);
        time += time_step;
        if (i % OutputStep == 0)
        ofs << time << "    "
            << expectations[1][0].item<double>() << "    "
            << expectations[1][1].item<double>() << "    "
            << expectations[2][0][0].item<double>() << "    "
            << expectations[2][0][1].item<double>() << "    "
            << expectations[2][1][1].item<double>() << '\n';
    }
    ofs.close();
}