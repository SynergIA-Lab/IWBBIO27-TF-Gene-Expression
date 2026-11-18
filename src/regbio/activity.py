from __future__ import annotations

import pandas as pd


def _run_decoupler_method(expression: pd.DataFrame, regulon: pd.DataFrame, method_name: str, tmin: int = 5) -> pd.DataFrame:
    import decoupler as dc

    method = getattr(dc.mt, method_name)
    estimates, _ = method(expression.T, regulon, tmin=tmin, verbose=False)
    return estimates.T


def ulm_activity(expression: pd.DataFrame, regulon: pd.DataFrame, tmin: int = 5) -> pd.DataFrame:
    return _run_decoupler_method(expression, regulon, "ulm", tmin=tmin)


def viper_activity(expression: pd.DataFrame, regulon: pd.DataFrame, tmin: int = 5) -> pd.DataFrame:
    return _run_decoupler_method(expression, regulon, "viper", tmin=tmin)


def tf_expression_control(expression: pd.DataFrame, tf_set: list[str]) -> pd.DataFrame:
    available = [tf for tf in tf_set if tf in expression.index]
    return expression.loc[available].copy()
