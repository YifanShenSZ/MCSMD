#include <torch/torch.h>

#include <CppLibrary/argparse.hpp>

#include <CppLibrary/utility.hpp>

argparse::ArgumentParser parse_args(const int & argc, const char ** & argv) {
    CL::utility::EchoCommand(argc, argv); std::cout << '\n';
    argparse::ArgumentParser parser("Multi-centred semiclassical Moyal dynamics");

    parser.add_argument("-q" ,"--<q>" , 1, false);
    parser.add_argument("-p" ,"--<p>" , 1, false);
    parser.add_argument("-qq","--<qq>", 1, false);
    parser.add_argument("-pp","--<pp>", 1, false);

    parser.parse_args(argc, argv);
    return parser;
}

int main(int argc, const char** argv) {
    argparse::ArgumentParser args = parse_args(argc, argv);

    double mass = 2000.0;

    at::Tensor expectations = at::empty(5, at::TensorOptions().dtype(torch::kFloat64));
    expectations[0] = args.retrieve<double>("<q>");
    expectations[1] = args.retrieve<double>("<p>");
    expectations[2] = args.retrieve<double>("<qq>");
    expectations[3] = 0.0;
    expectations[4] = args.retrieve<double>("<pp>");

    
}