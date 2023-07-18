from dataclasses import dataclass

import pandas as pd
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    plotter: stk.ProgressPlotter
    plot_data: pd.DataFrame
