"""
가상 웨이퍼 공정 데이터 생성기

- 40 lot x 25 wafer, wafer 당 9 site CD(Critical Dimension, nm) 측정
- wafer 당 die 수(400)와 불량 die 수 -> 수율(yield)
- 공정 이상 시나리오를 의도적으로 삽입
    * lot 25~ : 장비 PM 이후 CD 평균 +1.5nm step shift
    * lot 32~ : 국소 파티클 이슈로 불량률 증가
"""
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
N_LOTS = 40
WAFERS_PER_LOT = 25
SITES = 9
DIES_PER_WAFER = 400

CD_TARGET = 45.0      # nm
CD_SIGMA_SITE = 0.6   # wafer 내 site 간 산포
CD_SIGMA_WAFER = 0.4  # wafer 간 산포
CD_SIGMA_LOT = 0.3    # lot 간 산포

SHIFT_LOT = 25
SHIFT_SIZE = 1.5
DEFECT_LOT = 32

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def generate(seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for lot in range(1, N_LOTS + 1):
        lot_offset = rng.normal(0, CD_SIGMA_LOT)
        if lot >= SHIFT_LOT:
            lot_offset += SHIFT_SIZE
        base_p = 0.06 if lot < DEFECT_LOT else 0.10
        for wafer in range(1, WAFERS_PER_LOT + 1):
            wafer_offset = rng.normal(0, CD_SIGMA_WAFER)
            sites = CD_TARGET + lot_offset + wafer_offset + rng.normal(0, CD_SIGMA_SITE, SITES)
            # CD가 target에서 벗어날수록 불량률 증가
            p_fail = base_p + 0.01 * (sites.mean() - CD_TARGET) ** 2
            n_fail = rng.binomial(DIES_PER_WAFER, min(p_fail, 1.0))
            row = {
                "lot_id": f"LOT{lot:03d}",
                "lot_seq": lot,
                "wafer_id": wafer,
                **{f"cd_site{i + 1}": round(v, 3) for i, v in enumerate(sites)},
                "cd_mean": round(sites.mean(), 3),
                "cd_range": round(sites.max() - sites.min(), 3),
                "total_die": DIES_PER_WAFER,
                "fail_die": n_fail,
                "yield_pct": round(100 * (1 - n_fail / DIES_PER_WAFER), 2),
            }
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> Path:
    DATA_DIR.mkdir(exist_ok=True)
    df = generate()
    out = DATA_DIR / "wafer_data.csv"
    df.to_csv(out, index=False)
    print(f"[data] {len(df)} wafers -> {out}")
    return out


if __name__ == "__main__":
    main()
