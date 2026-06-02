#include <iostream>
#include <matplot/matplot.h>
#include "Labs/Special/Lab8/Tasks/ChannelStudy.h"

int main() {
    try {
        ChannelStudy::runAll();
        matplot::show();
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
}