# PyLinearCode

PyLinearCode 是一个纯 Python 的有限域线性码计算库雏形，目标是提供接近 Magma
线性码模块的易用 API，同时比 SageMath 更轻、更容易嵌入 Python 项目，并为后续
Numba/Cython/Rust 加速留下清晰边界。

> 现实一点说：当前版本不是完整 SageMath 替代品，而是一个可运行、可测试、可扩展的
> 开源内核。要“整体优于 SageMath/Magma”需要持续补齐大规模最小距离、同构判定、
> 自同构群、BCH/AG/Goppa 等高级算法。

## 已实现

- 有限域 `GF(p)` 与 `GF(p^m)`，扩域用多项式基和整数编码表示。
- 域元素加减乘除、逆元、幂、格式化、默认不可约多项式搜索。
- 有限域矩阵行化简、秩、零空间、正交补、矩阵乘法。
- 线性码核心：
  - `LinearCode`
  - `generator_matrix`
  - `parity_check_matrix`
  - `encode`
  - `contains`
  - `syndrome`
  - `dual`
  - `hull`
  - `puncture`
  - `shorten`
  - `extend_with_parity`
  - `minimum_distance`
  - `weight_distribution`
  - `syndrome_decode`
  - `systematic_form`
- 典型构造：
  - 零码、全集码、重复码、单校验码
  - Hamming 码
  - Reed-Solomon 码
  - 随机线性码
  - 循环码生成多项式构造
- 常见界：
  - Singleton
  - Hamming 球体积/上界检查
  - Griesmer 下界

## 快速开始

```python
from pylinearcode import GF, hamming_code, reed_solomon_code

F2 = GF(2)
C = hamming_code(F2, r=3)

print(C.parameters())          # [7, 4, 3] over GF(2)
print(C.minimum_distance())    # 3

word = C.encode([1, 0, 1, 1])
received = word[:]
received[2] ^= 1
decoded, error = C.syndrome_decode(received, max_weight=1)
assert decoded == word

F5 = GF(5)
RS = reed_solomon_code(F5, length=5, dimension=3)
print(RS.parameters())         # [5, 3, 3] over GF(5)
```

## 命令行

```bash
python -m pylinearcode hamming --q 2 --r 3
python -m pylinearcode rs --q 5 --n 5 --k 3
```

## 开发

```bash
python -m unittest discover -s tests
python -m pip install -e ".[dev]"
pytest
```

## 设计方向

PyLinearCode 的长期目标不是简单复刻 SageMath，而是做一个更适合现代 Python 生态的
编码理论计算库：

1. 核心算法保持小而清晰，默认无重依赖。
2. 热点路径保留整数矩阵接口，便于切换到 Numba、Cython、Rust 或 GPU 后端。
3. API 对齐 Magma/Sage 常用术语，但返回 Python 原生数据结构。
4. 大规模困难问题使用显式计算预算，避免无提示卡死。
5. 后续优先补齐信息集译码、Leon/Brouwer-Zimmermann 最小距离、码等价与自同构群。

