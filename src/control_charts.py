"""
웨이퍼 데이터 기반 SPC 관리도

1. I-MR 관리도   : lot 평균 CD (lot 간 변동을 반영한 관리한계)
2. 기존 X-bar/R  : 비교용 - lot 내 범위만으로 한계를 잡으면 오알람이 생기는 사례
3. p 관리도      : wafer 별 불량 die 비율 (n=400)
4. EWMA 관리도   : lot 평균 CD의 작은 shift 검출

관리한계는 Phase I 구간(정상 lot 1~20)으로 산출하고 이후 lot에 적용한다.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from font_setup import set_korean_font

set_korean_font()

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "wafer_data.csv"
FIG_DIR = ROOT / "figures"

PHASE1_LOTS = 20

# 관리도 상수
D2_MR, D4_MR = 1.128, 3.267            # moving range (n=2)
A2_25, D3_25, D4_25 = 0.153, 0.459, 1.541  # subgroup n=25


def western_electric(x: np.ndarray, cl: float, sigma: float) -> np.ndarray:
    """Western Electric 규칙 1~4 위반 index를 반환."""
    z = (x - cl) / sigma
    flag = np.zeros(len(x), dtype=bool)
    for i in range(len(x)):
        # Rule 1: 1점이 3σ 밖
        if abs(z[i]) > 3:
            flag[i] = True
        # Rule 2: 연속 3점 중 2점이 같은 쪽 2σ 밖
        if i >= 2:
            w = z[i - 2:i + 1]
            if (w > 2).sum() >= 2 or (w < -2).sum() >= 2:
                flag[i] = True
        # Rule 3: 연속 5점 중 4점이 같은 쪽 1σ 밖
        if i >= 4:
            w = z[i - 4:i + 1]
            if (w > 1).sum() >= 4 or (w < -1).sum() >= 4:
                flag[i] = True
        # Rule 4: 연속 8점이 중심선 한쪽
        if i >= 7:
            w = z[i - 7:i + 1]
            if (w > 0).all() or (w < 0).all():
                flag[i] = True
    return np.where(flag)[0]


def _first(idx: np.ndarray, labels: np.ndarray):
    return int(labels[idx[0]]) if len(idx) else None


def _plot_limits(ax, x, cl, ucl, lcl, label, phase1_x=PHASE1_LOTS + 0.5):
    ax.plot(x, np.broadcast_to(ucl, x.shape), "r--", lw=1, label="UCL / LCL")
    ax.plot(x, np.broadcast_to(lcl, x.shape), "r--", lw=1)
    ax.axhline(cl, color="g", lw=1, label="CL")
    ax.axvline(phase1_x, color="gray", ls=":", lw=1)
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)


def lot_means(df: pd.DataFrame):
    s = df.groupby("lot_seq")["cd_mean"].mean()
    return s.index.to_numpy(), s.to_numpy()


def imr_chart(df: pd.DataFrame) -> dict:
    lots, x = lot_means(df)
    mr = np.abs(np.diff(x))

    cl = x[:PHASE1_LOTS].mean()
    mrbar = mr[:PHASE1_LOTS - 1].mean()
    sigma = mrbar / D2_MR
    ucl, lcl = cl + 3 * sigma, cl - 3 * sigma

    viol = western_electric(x, cl, sigma)
    ooc = np.where((x > ucl) | (x < lcl))[0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax1.plot(lots, x, "o-", ms=4, color="tab:blue")
    ax1.plot(lots[viol], x[viol], "o", ms=9, mfc="none", mec="orange", label="WE rule 위반")
    ax1.plot(lots[ooc], x[ooc], "rx", ms=10, label="관리한계 이탈")
    _plot_limits(ax1, lots, cl, ucl, lcl, "Lot 평균 CD (nm)")
    ax1.set_title("I-MR 관리도 (lot 평균 CD)")
    ax1.legend(loc="upper left", fontsize=8)

    ax2.plot(lots[1:], mr, "o-", ms=4, color="tab:purple")
    _plot_limits(ax2, lots[1:], mrbar, D4_MR * mrbar, 0.0, "Moving Range (nm)")
    ax2.set_xlabel("Lot sequence  (점선 왼쪽: Phase I 한계 산출 구간)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "imr_chart.png", dpi=130)
    plt.close(fig)

    return {"CL": cl, "UCL": ucl, "LCL": lcl, "sigma_lot": sigma,
            "first_OOC_lot": _first(ooc, lots), "first_WE_lot": _first(viol, lots)}


def xbar_r_chart_naive(df: pd.DataFrame) -> dict:
    """lot 내 wafer 범위(R)로만 한계를 잡은 X-bar/R. lot 간 변동이 빠져 한계가 과도하게 좁다."""
    g = df.groupby("lot_seq")["cd_mean"]
    lots, xbar = lot_means(df)
    rng = (g.max() - g.min()).to_numpy()

    xbarbar = xbar[:PHASE1_LOTS].mean()
    rbar = rng[:PHASE1_LOTS].mean()
    ucl, lcl = xbarbar + A2_25 * rbar, xbarbar - A2_25 * rbar
    ooc = np.where((xbar > ucl) | (xbar < lcl))[0]
    false_alarms = int((ooc < PHASE1_LOTS).sum())

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax1.plot(lots, xbar, "o-", ms=4, color="tab:blue")
    ax1.plot(lots[ooc], xbar[ooc], "rx", ms=10, label="관리한계 이탈")
    _plot_limits(ax1, lots, xbarbar, ucl, lcl, "CD X-bar (nm)")
    ax1.set_title(f"(비교) X-bar/R 관리도 - lot 내 범위 기반 한계 → Phase I 오알람 {false_alarms}건")
    ax1.legend(loc="upper left", fontsize=8)

    ax2.plot(lots, rng, "o-", ms=4, color="tab:purple")
    _plot_limits(ax2, lots, rbar, D4_25 * rbar, D3_25 * rbar, "CD Range (nm)")
    ax2.set_xlabel("Lot sequence")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "xbar_r_chart_naive.png", dpi=130)
    plt.close(fig)

    return {"UCL": ucl, "LCL": lcl, "phase1_false_alarms": false_alarms}


def p_chart(df: pd.DataFrame) -> dict:
    p = (df["fail_die"] / df["total_die"]).to_numpy()
    n = df["total_die"].to_numpy()
    seq = np.arange(1, len(df) + 1)
    n1 = PHASE1_LOTS * df["wafer_id"].max()

    pbar = df["fail_die"].iloc[:n1].sum() / df["total_die"].iloc[:n1].sum()
    s = np.sqrt(pbar * (1 - pbar) / n)
    ucl, lcl = pbar + 3 * s, np.clip(pbar - 3 * s, 0, None)
    ooc = np.where((p > ucl) | (p < lcl))[0]
    lot_p = df.groupby("lot_seq")[["fail_die", "total_die"]].sum()
    lot_p = (lot_p["fail_die"] / lot_p["total_die"]).to_numpy()

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(seq, p * 100, ".", ms=3, color="tab:brown", alpha=0.6, label="wafer 불량률")
    wpl = df["wafer_id"].max()
    edges = np.arange(len(lot_p) + 1) * wpl + 1
    ax.stairs(lot_p * 100, edges, color="k", lw=1, label="lot 평균")
    ax.plot(seq[ooc], p[ooc] * 100, "rx", ms=6, label="관리한계 이탈")
    _plot_limits(ax, seq, pbar * 100, ucl * 100, lcl * 100, "불량 die 비율 (%)", phase1_x=n1 + 0.5)
    ax.set_title("p 관리도 (wafer 별 불량률, n=400 die)")
    ax.set_xlabel("Wafer sequence")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p_chart.png", dpi=130)
    plt.close(fig)

    lot_of = df["lot_seq"].to_numpy()
    ooc_p1 = int((ooc < n1).sum())
    return {"p_CL_pct": pbar * 100, "phase1_OOC_wafers": ooc_p1,
            "phase2_OOC_wafers": int(len(ooc) - ooc_p1),
            "first_phase2_OOC_lot": _first(ooc[ooc >= n1], lot_of)}


def ewma_chart(df: pd.DataFrame, lam: float = 0.2, L: float = 3.0) -> dict:
    """lot 평균 CD에 대한 EWMA 관리도 (σ는 MR 기반)."""
    lots, x = lot_means(df)
    mu0 = x[:PHASE1_LOTS].mean()
    sigma = np.abs(np.diff(x[:PHASE1_LOTS])).mean() / D2_MR

    z = np.empty_like(x)
    prev = mu0
    for k, xi in enumerate(x):
        prev = lam * xi + (1 - lam) * prev
        z[k] = prev
    i = np.arange(1, len(x) + 1)
    width = L * sigma * np.sqrt(lam / (2 - lam) * (1 - (1 - lam) ** (2 * i)))
    ucl, lcl = mu0 + width, mu0 - width
    ooc = np.where((z > ucl) | (z < lcl))[0]

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(lots, x, "o", ms=4, color="lightgray", label="lot 평균 CD")
    ax.plot(lots, z, "o-", ms=3, lw=1.2, color="tab:blue", label=f"EWMA (λ={lam})")
    ax.plot(lots, ucl, "r--", lw=1, label="UCL / LCL")
    ax.plot(lots, lcl, "r--", lw=1)
    ax.axhline(mu0, color="g", lw=1, label="CL")
    ax.axvline(PHASE1_LOTS + 0.5, color="gray", ls=":", lw=1)
    if len(ooc):
        ax.axvline(lots[ooc[0]], color="orange", lw=1, label=f"첫 알람 (lot {lots[ooc[0]]})")
    ax.set_title("EWMA 관리도 (lot 평균 CD)")
    ax.set_xlabel("Lot sequence")
    ax.set_ylabel("CD (nm)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ewma_chart.png", dpi=130)
    plt.close(fig)

    return {"lambda": lam, "first_alarm_lot": _first(ooc, lots)}


def main() -> dict:
    FIG_DIR.mkdir(exist_ok=True)
    df = pd.read_csv(DATA)
    res = {
        "imr": imr_chart(df),
        "xbar_r_naive": xbar_r_chart_naive(df),
        "p": p_chart(df),
        "ewma": ewma_chart(df),
    }
    for k, v in res.items():
        print(f"[chart] {k}: " + ", ".join(
            f"{a}={b:.3f}" if isinstance(b, float) else f"{a}={b}" for a, b in v.items()))
    return res


if __name__ == "__main__":
    main()
