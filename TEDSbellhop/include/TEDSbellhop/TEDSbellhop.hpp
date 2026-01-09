#pragma once

/**
 * @file TEDSbellhop.hpp
 * @brief TEDSbellhop：对 bellhopcxx / bellhopcuda 的“纯内存输入、C++17 友好”封装。
 *
 * ## 设计目标
 * - **只暴露 tedsbellhop 命名空间内的类型**：对外不直接暴露底层 bellhop 的 bhc:: 结构体/参数，避免 ABI/版本耦合。
 * - **不读取任何输入文件**：强制走 no-file 模式（所有环境、边界、SSP、射线、网格参数均由内存结构体提供）。
 * - **同步 + 异步两套调用方式**：既能阻塞式计算，也能启动后台任务并轮询进度。
 *
 * ## 使用方式概览
 * 1. 构造并填充 @ref tedsbellhop::TL2DInput
 * 2. （可选）调用 @ref tedsbellhop::validate 做参数校验（也可以让 compute/start 在内部校验）
 * 3. 调用：
 *    - 同步：@ref tedsbellhop::compute_tl_2d
 *    - 异步：@ref tedsbellhop::start_tl_2d -> @ref tedsbellhop::TL2DJob
 *
 * ## 单位与网格约定（非常重要）
 * - 深度统一使用 **米 (m)**。
 * - 距离（水平 range）既可用 m 也可用 km：
 *   - TL 输出网格使用 @ref tedsbellhop::Positions2D::receiver_ranges + @ref tedsbellhop::Positions2D::receiver_ranges_in_km
 *   - Q-SSP 的 range 网格使用 @ref tedsbellhop::SSPQuad::ranges + @ref tedsbellhop::SSPQuad::ranges_in_km
 * - Q-SSP 声速矩阵布局：@ref tedsbellhop::SSPQuad::c_mps 采用 **行优先（depth 行）**：
 *   `c_mps[iz * Nr + ir]`
 *
 * ## 与 bellhop 选项字符串（OPTIONS / RunType）关系
 * - bellhop 传统通过 7 字符 RunType + OPTIONS 控制算法分支。
 * - 本库提供：
 *   - “安全易用”的枚举组合（@ref tedsbellhop::RunType）
 *   - 以及用于完全复刻原行为的 raw 字符串（@ref tedsbellhop::RunType::raw + use_raw）。
 * - @ref tedsbellhop::TL2DInput::options 对应 env 的 OPTIONS1（如 "QVMT"），但库保证不读取任何文件。
 *
 * ## 线程安全与并行
 * - 多个 @ref tedsbellhop::TL2DJob 可并行启动（典型用于多剖面/多源深度批处理）。
 * - 具体并行度与底层实现（CPU/CUDA）及其内部资源管理有关。
 *
 * ## 说明
 * 本库内部仍会链接到底层 bellhop 静态库实现（bellhopcxxstatic / bellhopcudastatic）。
 */



#include "TEDSbellhop_Export.h"
#include "TEDSbellhop_Data.h"
#include "TEDSbellhop_API.h"


