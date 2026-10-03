# Morley-LSG-3D

三维圆环柱 pull-out 算例，使用 Morley 元与 JS 投影格式。
本上传包来自 `Morley_pullout_stage5f_fresh_N8_exactH.zip`（2026-10-03 保存的版本）。
全部 36 个 MATLAB 源文件、6 份原有说明和 1 个基准 CSV 均原样保留。

当前主入口是 **Stage 5F：分两次 MATLAB 会话运行的 N=8、exact-H 能量坐标 GMRES 求解**。
网格由 `mesh/build_annular_structured_tet_mesh.m` 在 MATLAB 内生成，当前版本不需要 Gmsh。
这是 JS pull-out 项目；包内没有 SP 格式或之前 Gmsh 制造解项目的入口。

## 上传到 GitHub

1. 解压本压缩包，打开 `Morley-LSG-3D` 文件夹。
2. 创建普通 GitHub 仓库，例如 `Morley-LSG-3D`；私人项目可选择 Private。
3. 将此文件夹**里面的文件和子目录**上传到仓库根目录。
   根目录应直接出现本 `README.md`、`setup_stage5.m` 和各个 `RUN_*.m`。
4. 一起上传 `.gitignore`、`.gitattributes` 和 `AGENTS.md`。不要仅将 ZIP 作为一个文件上传。
5. 提交代码后，在 Codex 的仓库选择界面选择该仓库。

`.gitignore` 会排除运行输出、缓存和生成的 `.mat` 文件，
并特意保留 `results/stage5a_verified_N4_baseline.csv`。
这个 CSV 是主程序读取的输入，不能作为普通结果文件删掉。
使用网页上传时，也请只上传本包内容；以后不要把新生成的大型输出手动选入上传列表。

## 目录与文件

| 路径 | 作用 |
| --- | --- |
| `setup_stage5.m` | 设置项目路径，创建输出目录并返回配置 |
| `config/` | 几何、材料、网格和求解器参数 |
| `mesh/`、`fem/` | MATLAB 网格生成、Morley 自由度和基函数 |
| `assembly/` | 体项、JS 投影面项和非齐次边界处理 |
| `reference/` | 一维径向 Aifantis 参考解 |
| `solver/` | 直接求解、块预条件与 exact-H GMRES |
| `post/`、`utils/` | 误差、径向剖面和辅助函数 |
| `results/stage5a_verified_N4_baseline.csv` | 随包保留的 N=4 历史基准输入 |
| `README_STAGE5*.txt` | 原始阶段说明；当前运行顺序以本页和 Stage5F 说明为准 |
| `AGENTS.md` | 给 Codex 的项目工作约定 |

## 当前默认参数

以 `config/stage5_pullout_config.m` 为准：

| 参数 | 数值 |
| --- | --- |
| 内半径、外半径、高度 | `RI=0.01`、`RE=1`、`H=0.05` |
| 内侧位移参数 | `ub=0.005` |
| 材料与长度参数 | `mu=1`、`iota=0.1` |
| 比较的 lambda | `1`、`1e5` |
| 配置内已有网格层级 | `N=4`、`N=8` |

这些参数属于当前 Stage5F 源码，不沿用早期工程算例的其他材料参数。
修改几何、材料、iota 或网格后，应重新生成对应基准和初值，
不能直接将随包 CSV 或旧 `.mat` 当成新参数下的结果。

## 运行前检查

进入项目根目录，在 MATLAB 命令窗口运行：

```matlab
cfg = setup_stage5();
disp(cfg.codeVersion)
which bvp4c -all
which gmres -all
which ichol -all
which symbfact -all
```

需要实际可运行且已授权的 MATLAB。用户此前环境为 R2024b；
本次打包环境没有 MATLAB，因此没有重新验证 MATLAB 版本兼容性或求解结果。
包内不包含 MATLAB、授权文件或 Agentic Toolkit 安装文件。

## 推荐运行顺序

### 1. 新计算环境先做 N=4 检查

```matlab
RUN_STAGE5C_N4_CALIBRATION
```

此入口比较 N=4 的直接 LU 与块预条件 GMRES，供检查求解环境和误差差异。
它是诊断程序；执行完成本身不代表所有精度指标已经满足要求。
若只是阅读或修改代码，不必运行数值求解。

### 2. 准备 N=8、lambda=1 的初值

```matlab
RUN_STAGE5F_PREPARE_X0
```

成功后会生成：

```text
results/stage5f_N8_lambda1_x0.mat
```

这是下一步必需的本机运行产物，初始上传包没有包含它。
请保留文件，然后**完全关闭 MATLAB**。

### 3. 用新的 MATLAB 会话计算 N=8、lambda=1e5

重新打开 MATLAB，进入同一个项目目录，只运行：

```matlab
RUN_STAGE5F_N8_HIGHLAMBDA
```

该程序读取上一步初值以及随包 N=4 基准 CSV，完成求解后输出：

```text
results/stage5f_N8_exactH_fresh.csv
results/stage5f_N8_exactH_fresh.mat
```

如果触发内存停止或精度不满足，应按程序报告记录为停止或诊断结果，
不能据此宣称 N=8 求解成功或 lambda 鲁棒性已经得到验证。

## 云端运行与内存

GitHub 仓库保存源码。实际执行上述程序的环境还需具备 MATLAB 和有效授权。
本地 Windows 上安装的工具及其绝对路径，不会随本压缩包成为云端的 MATLAB 安装。

上传包没有启动时自动运行 N=8 的脚本，也没有配置自动计算的工作流。
请先确认环境，再按上面的顺序执行。

当前源码用 `ispc` 分支读取 Windows 可用物理内存；
Linux 上不会执行这部分检测，但仍有 Cholesky 填充估计和配置中的固定上限检查。
这些估计不等于完整峰值内存保证。云端运行前还需核实实际内存及容器限制，
不要仅因使用云端就提高网格层级或放宽停止阈值。

## 历史入口说明

- `RUN_STAGE5_LAMBDA_N8.m` 与 `RUN_STAGE5_PLOT.m` 属于早期 Stage5B 路线。
  其中绘图程序读取 `stage5b_lambda_N8.mat`，并不读取 Stage5F 输出。
- `RUN_STAGE5C_SOLVER_DIAGNOSTIC.m` 在当前源码中会转到 `RUN_STAGE5C_N4_CALIBRATION.m`。
- `RUN_STAGE5E_ENERGY_TRANSFORM.m` 会先做 N=4 校准，再尝试 N=8；它不是仅运行 N=4 的入口。
- 历史 `README_STAGE5D.txt` 提到的 `RUN_STAGE5D_N4_FULLH_CALIBRATION.m`
  原始压缩包没有提供，本包也没有补造该程序；不要按那份历史说明启动。

## 本次打包的核验范围

交付前核对源码与原始压缩包的字节一致性、必需输入、Git 忽略规则和 ZIP 完整性。
本次没有重新执行 MATLAB 数值计算；原有文档与基准中的数值是继承数据。
`docs/SOURCE_SHA256SUMS.txt` 列出原始 43 个文件的 SHA-256，可用于核对。
