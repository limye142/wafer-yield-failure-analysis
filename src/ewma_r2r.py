"""
EWMA Run-to-Run(R2R) 보정 시뮬레이션

공정 모델 (예: 식각/노광 CD 제어)
    y_t = a + b * u_t + d_t + e_t
      y_t : 측정 CD (nm)
      u_t : 레시피 조절 변수 (예: dose / etch time)
      b   : 공정 gain (컨트롤러는 추정치 b_hat 사용)
      d_t : 외란 = 장비 drift + PM 후 step shift
      e_t : 측정/공정 노이즈

EWMA 컨트롤러
    a_hat_t = λ (y_t - b_hat * u_t) + (1 - λ) a_hat_{t-1}
    u_{t+1} = (T - a_hat_t) / b_hat

λ 값에 따라 보정 속도와 노이즈 민감도가 trade-off 관계임을 비교한다.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from font_setup import set_korean_font

set_korean_font()

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"

TARGET = 45.0
A_TRUE = 30.0
B_TRUE = 1.5
B_HAT = 1.3          # 모델 gain 오차(실제보다 ~13% 작게 추정)
NOISE = 0.35
DRIFT_PER_RUN = 0.02  # run 당 장비 drift (nm)
STEP_RUN = 120        # PM 이후 step shift 발생 시점
STEP_SIZE = 1.5
N_RUNS = 250
USL, LSL = TARGET + 1.5, TARGET - 1.5

LAMBDAS = [0.1, 0.3, 0.6]


def disturbance(n: int, rng: np.random.Generator):
    t = np.arange(n)
    d = DRIFT_PER_RUN * t + np.where(t >= STEP_RUN, STEP_SIZE, 0.0)
    e = rng.normal(0, NOISE, n)
    return d, e


def simulate(lam: float | None, d: np.ndarray, e: np.ndarray) -> np.ndarray:
    """lam=None 이면 보정 없이 초기 레시피 고정(open-loop)."""
    u = (TARGET - A_TRUE) / B_TRUE  # 초기 레시피는 qual 시점에 target으로 맞춰져 있다고 가정
    a_hat = TARGET - B_HAT * u  # 컨트롤러 모델과 일관된 초기 절편 추정
    y = np.empty(len(d))
    for t in range(len(d)):
        y[t] = A_TRUE + B_TRUE * u + d[t] + e[t]
        if lam is not None:
            a_hat = lam * (y[t] - B_HAT * u) + (1 - lam) * a_hat
            u = (TARGET - a_hat) / B_HAT
    return y


def metrics(y: np.ndarray) -> dict:
    mu, s = y.mean(), y.std(ddof=1)
    return {
        "mean": mu,
        "std": s,
        "RMSE_vs_target": np.sqrt(np.mean((y - TARGET) ** 2)),
        "Cpk": min(USL - mu, mu - LSL) / (3 * s),
        "spec_out_pct": 100 * np.mean((y > USL) | (y < LSL)),
        # step 이후 20 run 내 평균 편차: 회복 속도 지표
        "post_step_bias": np.mean(y[STEP_RUN:STEP_RUN + 20] - TARGET),
    }


def lambda_sweep(n_rep: int = 200, seed: int = 7) -> pd.DataFrame:
    """λ 별 RMSE를 Monte-Carlo 반복으로 평가해 최적 λ 범위를 확인."""
    rng = np.random.default_rng(seed)
    lams = np.round(np.arange(0.05, 1.0001, 0.05), 2)
    out = np.zeros((n_rep, len(lams)))
    for r in range(n_rep):
        d, e = disturbance(N_RUNS, rng)
        for j, lam in enumerate(lams):
            y = simulate(lam, d, e)
            out[r, j] = np.sqrt(np.mean((y - TARGET) ** 2))
    return pd.DataFrame({"lambda": lams, "RMSE_mean": out.mean(0), "RMSE_std": out.std(0)})


def main() -> pd.DataFrame:
    FIG_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(2024)
    d, e = disturbance(N_RUNS, rng)

    runs = {"No control": simulate(None, d, e)}
    for lam in LAMBDAS:
        runs[f"EWMA λ={lam}"] = simulate(lam, d, e)

    summary = pd.DataFrame({k: metrics(v) for k, v in runs.items()}).T
    summary.index.name = "scenario"
    summary.round(4).to_csv(DATA_DIR / "ewma_r2r_summary.csv")
    pd.DataFrame(runs).round(4).rename_axis("run").to_csv(DATA_DIR / "ewma_r2r_runs.csv")

    # --- 시계열 비교 ---
    t = np.arange(N_RUNS)
    colors = ["lightgray", "tab:orange", "tab:blue", "tab:green"]
    fig, ax = plt.subplots(figsize=(11, 5))
    for (name, y), c in zip(runs.items(), colors):
        ax.plot(t, y, lw=1.6 if name == "No control" else 1.1, color=c, label=name)
    ax.axhline(TARGET, color="k", lw=1, label="Target")
    ax.axhspan(LSL, USL, color="green", alpha=0.06)
    ax.axhline(USL, color="r", ls="--", lw=1, label="USL / LSL")
    ax.axhline(LSL, color="r", ls="--", lw=1)
    ax.axvline(STEP_RUN, color="gray", ls=":", lw=1)
    ax.text(STEP_RUN + 2, ax.get_ylim()[1] * 0.995, f"PM step +{STEP_SIZE}nm", fontsize=8, va="top")
    ax.set_title(f"EWMA R2R 보정 시뮬레이션 (drift {DRIFT_PER_RUN}nm/run, gain 모델비 b_hat/b={B_HAT / B_TRUE:.2f})")
    ax.set_xlabel("Run")
    ax.set_ylabel("CD (nm)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ewma_r2r_simulation.png", dpi=130)
    plt.close(fig)

    # --- λ sweep ---
    sweep = lambda_sweep()
    sweep.round(4).to_csv(DATA_DIR / "ewma_lambda_sweep.csv", index=False)
    best = sweep.loc[sweep["RMSE_mean"].idxmin()]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.errorbar(sweep["lambda"], sweep["RMSE_mean"], yerr=sweep["RMSE_std"], fmt="o-", ms=4, capsize=3)
    ax.axvline(best["lambda"], color="orange", ls="--", lw=1, label=f"최적 λ≈{best['lambda']:.2f}")
    ax.set_xlabel("λ (EWMA weight)")
    ax.set_ylabel("RMSE vs target (nm)")
    ax.set_title("λ 별 보정 성능 (Monte-Carlo 200회)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ewma_lambda_sweep.png", dpi=130)
    plt.close(fig)

    print("[r2r] summary\n" + summary.round(3).to_string())
    print(f"[r2r] best lambda ≈ {best['lambda']:.2f} (RMSE {best['RMSE_mean']:.3f} nm)")
    return summary


if __name__ == "__main__":
    main()
