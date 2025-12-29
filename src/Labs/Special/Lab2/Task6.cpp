#include <iostream>
#include "Labs/Special/Additional/NoiseGenerator/DichotomicNoise.h"
#ifdef _WIN32
#include <windows.h>
#endif

int main() {

#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    std::cout << "=== СИММЕТРИЧНЫЙ ШУМ ===" << std::endl;
    {
        DichotomicNoise gen(1.0, 0.05, 0.001);
        auto prof = gen.generate(100'000'000, 100'000, true);

        std::cout << "a=" << gen.a() << ", b=" << gen.b()
                  << ", gamma_a=" << gen.gamma_a()
                  << ", gamma_b=" << gen.gamma_b() << "\n";
        std::cout << "Теория: <σ> = " << gen.mean_theoretical() << "\n";
        std::cout << "Числено: <σ> = " << prof.m1
                  << ", <σ²> = " << prof.m2 << "\n\n";
    }

    std::cout << "=== АСИММЕТРИЧНЫЙ ШУМ С НУЛЕВЫМ СРЕДНИМ ===" << std::endl;
    {
        auto gen = DichotomicNoise::ZeroMeanAsymmetric(2.0, -1.0, 1.0, 0.001);
        auto prof = gen.generate(100'000'000, 100'000, true);

        std::cout << "a=" << gen.a() << ", b=" << gen.b()
                  << ", gamma_a=" << gen.gamma_a()
                  << ", gamma_b=" << gen.gamma_b() << "\n";
        std::cout << "Теория: <σ> = " << gen.mean_theoretical() << "\n";
        std::cout << "Числено: <σ> = " << prof.m1
                  << ", <σ²> = " << prof.m2 << "\n\n";
    }

    std::cout << "=== ЕЩЁ ОДИН АСИММЕТРИЧНЫЙ С НУЛЕВЫМ СРЕДНИМ ===" << std::endl;
    {
        auto gen = DichotomicNoise::ZeroMeanAsymmetric(6.0, -3.0, 0.5, 0.001);
        auto prof = gen.generate(1'000'000'000, 1'000'000, true);

        std::cout << "a=" << gen.a() << ", b=" << gen.b()
                  << ", gamma_a=" << gen.gamma_a()
                  << ", gamma_b=" << gen.gamma_b() << "\n";
        std::cout << "Теория: <σ> = " << gen.mean_theoretical() << "\n";
        std::cout << "Числено: <σ> = " << prof.m1
                  << ", <σ²> = " << prof.m2 << "\n\n";
    }
}