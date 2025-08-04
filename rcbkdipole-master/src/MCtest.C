#include <TF1.h>
#include <Math/GSLIntegrator.h>
#include <iostream>

// 定义被积函数
double f(double x) {
    return x;
}

// 包装器类，用于与ROOT::Math::GSLIntegrator接口
class FuncWrapper {
public:
    double operator()(double x) const {
        return f(x);
    }
};

int main() {
    // 创建函数包装器
    FuncWrapper func;

    // 创建GSL积分器对象
    ROOT::Math::GSLIntegrator integrator;

    // 设置积分参数
    integrator.SetRelTolerance(1e-6); // 设置相对误差容限

    // 执行积分
    double result, error;
    integrator.Integral(func, 0, 1, result, error);

    // 输出结果
    std::cout << "Integral of y = x from 0 to 1 is " << result << " +/- " << error << std::endl;

    return 0;
}