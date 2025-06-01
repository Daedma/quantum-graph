# Quantum Graph

Библиотека для численного поиска собственных значений и собственных функций квантовых графов с тремя рёбрами и произвольными потенциалами.

## Возможности

- Поиск собственных значений (спектра) для заданного квантового графа.
- Вычисление собственных функций.
- Поддержка различных граничных условий (Дирихле, Неймана).
- Гибкая работа с потенциалами на рёбрах.

## Установка

1. Клонируйте репозиторий с подмодулями:
   ```
   git clone --recurse-submodules https://github.com/Daedma/quantum-graph.git
   ```
2. Перейдите в директорию проекта:
   ```
   cd quantum-graph
   ```
3. Соберите проект с помощью CMake:
   ```
   mkdir build
   cd build
   cmake ..
   cmake --build .
   ```

**Требования:**
- CMake >= 3.26
- Компилятор с поддержкой C++17
- Boost (odeint)
- Подмодуль Mathter (автоматически подтягивается)

## Пример использования

```cpp
#include "QuantumGraph.hpp"
#include <cmath>
#include <iostream>

int main() {
    // Определяем потенциалы на рёбрах
    auto q1 = [](double x) { return 0.0; };
    auto q2 = [](double x) { return x; };
    auto q3 = [](double x) { return std::sin(x); };

    // Создаём граф с граничными условиями по умолчанию
    QuantumGraph graph(q1, q2, q3);

    // Находим собственные значения на интервале [0, 20]
    auto eigenvalues = graph.calcEigenvalues(0.0, 20.0);

    std::cout << "Eigenvalues:\n";
    for (double val : eigenvalues) {
        std::cout << val << std::endl;
    }

    // Вычисляем собственную функцию для первого собственного значения
    if (!eigenvalues.empty()) {
        auto eigenfunctions = graph.calcEigenfunction(eigenvalues[0], 100, 1e-6);
        // eigenfunctions - вектор функций, которые можно использовать для анализа
    }
    return 0;
}
```

## Лицензия

MIT License

