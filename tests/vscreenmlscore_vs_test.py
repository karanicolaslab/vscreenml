import pandas as pd
from vscreenmlscore_vs import get_vs_score


def test_get_vs_score():
    dir = "/Users/yusufadeshina/vscreenml/data"
    result = get_vs_score(dir)
    assert isinstance(result, pd.DataFrame)
