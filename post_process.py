#!/usr/bin/python3

import sys
from typing import Callable
import pandas as pd
import matplotlib.pyplot as plt

# =================  Parser  =================
class RuntimeVariable:
    def __init__(self, name: str, prefix: str, type: Callable = float) -> None:
        self.name = name
        self.prefix = prefix
        self.type = type


# ragulatedPatchTargetValue_prefix = "Regulator: targetValue = "
# regulatedPatchValue_prefix = "Regulator: currentRegulatedPatchValue = "
# error_prefix = "Regulator: error = "
# outputSignal_prefix = "Regulator: outputSignal = "

RUNTIME_VARIABLES = [
    RuntimeVariable(name="time", prefix="Time = "),
    RuntimeVariable(
        name="regulatedPatchValue", prefix="Regulator: currentRegulatedPatchValue = "
    ),
    RuntimeVariable(name="targetValue", prefix="Regulator: targetValue = "),
    RuntimeVariable(name="inletValue", prefix="Regulator: value at inlet = "),
    RuntimeVariable(name="error", prefix="Regulator: error = "),
    RuntimeVariable(name="outputSignal", prefix="Regulator: outputSignal = "),
]


def dataframe_from_logs(file: str) -> pd.DataFrame:
    data = {var.name: [] for var in RUNTIME_VARIABLES}
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            for var in RUNTIME_VARIABLES:
                if line.startswith(var.prefix):
                    value = var.type(line.removeprefix(var.prefix).strip())
                    data[var.name].append(value)

    df = pd.DataFrame(data, columns=[var.name for var in RUNTIME_VARIABLES])
    return df

def plot_results(df: pd.DataFrame) -> None:
    t = df["time"]

    ax1 = plt.subplot(311)
    plt.plot(t, df["inletValue"])
    plt.setp(ax1.get_xticklabels())
    plt.ylabel("inlet T [C]")

    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(t, df["regulatedPatchValue"])
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.ylabel("target T [C]")

    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(t, df["error"])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.ylabel("error")

    plt.xlabel("time [s]")
    plt.show()

# ======================================================
log_file = sys.argv[1]

df = dataframe_from_logs(log_file)

try:
    # write csv if out file is specified
    out_file = sys.argv[2]
    df.to_csv(out_file, float_format="%.2f", index=False, sep=",")
except IndexError:
    # else plot
    plot_results(df)