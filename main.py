"""데이터 생성 -> 관리도 -> EWMA R2R 시뮬레이션을 순서대로 실행."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import control_charts  # noqa: E402
import ewma_r2r  # noqa: E402
import generate_wafer_data  # noqa: E402

if __name__ == "__main__":
    generate_wafer_data.main()
    control_charts.main()
    ewma_r2r.main()
