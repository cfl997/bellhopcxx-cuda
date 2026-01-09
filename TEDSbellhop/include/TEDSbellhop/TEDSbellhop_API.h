#pragma once

#include "TEDSbellhop_Export.h"
#include "TEDSbellhop_Data.h"

// -----------------------------
// API
// -----------------------------
namespace tedsbellhop {

// -----------------------------
// 异步任务（Job）接口
// -----------------------------

struct TL2DJobImpl;

class TEDSBELLHOP_API TL2DJob {
public:
    TL2DJob() = default;

    // 可移动，不可拷贝
    TL2DJob(TL2DJob &&) noexcept;
    TL2DJob &operator=(TL2DJob &&) noexcept;
    TL2DJob(const TL2DJob &)            = delete;
    TL2DJob &operator=(const TL2DJob &) = delete;

    ~TL2DJob();

    // 当前进度：0~100；未开始/不可用返回 -1
    int progress() const;

    // 是否完成（成功或失败都算完成）
    bool ready() const;

    // 若失败，返回错误信息；未失败返回空串
    std::string error() const;

    // 等待完成并取结果；如果任务失败会抛 tedsbellhop::Error
    TL2DResult get();

    // 主动等待完成（不取结果）
    void wait();

    // 是否持有有效任务
    explicit operator bool() const noexcept { return static_cast<bool>(impl_); }

private:
    friend TEDSBELLHOP_API TL2DJob start_tl_2d(const TL2DInput &in);
    explicit TL2DJob(std::shared_ptr<TL2DJobImpl> impl) : impl_(std::move(impl)) {}

    std::shared_ptr<TL2DJobImpl> impl_;
};

TEDSBELLHOP_API void validate(const TL2DInput &in);

// 同步（阻塞）接口
TEDSBELLHOP_API TL2DResult compute_tl_2d(const TL2DInput &in);

// 异步接口：启动计算并返回 job
TEDSBELLHOP_API TL2DJob start_tl_2d(const TL2DInput &in);

} // namespace tedsbellhop