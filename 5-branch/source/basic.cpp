#include <chrono>
#include <random>

std::default_random_engine generator;

void initialize_random() {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator = std::default_random_engine(seed);
}