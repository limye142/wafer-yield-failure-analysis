"""그래프 한글 표시용 폰트 설정 (macOS / Windows / Linux)."""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib import font_manager

CANDIDATES = ["AppleGothic", "Malgun Gothic", "NanumGothic", "Noto Sans CJK KR", "Noto Sans KR"]


def set_korean_font() -> str | None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in CANDIDATES:
        if name in available:
            plt.rcParams["font.family"] = name
            plt.rcParams["axes.unicode_minus"] = False
            return name
    return None
