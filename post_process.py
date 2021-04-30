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

RUNTIME_VARIABLES = [
    RuntimeVariable(name="time", prefix="Time = "),
    RuntimeVariable(
        name="sensorValue", prefix="Regulator: sensorValue = "
    ),
    RuntimeVariable(name="targetValue", prefix="Regulator: targetValue = "),
    RuntimeVariable(name="error", prefix="Regulator: error = "),
    RuntimeVariable(name="outputSignal", prefix="Regulator: outputSignal = "),
]

INLET_VAR = RuntimeVariable(name="inletValue", prefix="Regulator: value at inlet = ")


def dataframe_from_logs(file: str) -> pd.DataFrame:

    with open(file, "r", encoding="utf-8") as f:
        content = f.read()

    # Add inlet variable if exist in logs
    if INLET_VAR.prefix in content:
        RUNTIME_VARIABLES.append(INLET_VAR)

    data = {var.name: [] for var in RUNTIME_VARIABLES}
    for line in content.splitlines():
        for var in RUNTIME_VARIABLES:
            is_updated = var.name != "time" and len(data[var.name]) == len(data["time"])
            if line.startswith(var.prefix) and not is_updated:
                    value = var.type(line.removeprefix(var.prefix).strip())
                    data[var.name].append(value)

    df = pd.DataFrame(data, columns=[var.name for var in RUNTIME_VARIABLES])
    return df

def plot_results(df: pd.DataFrame) -> None:
    t = df["time"]

    has_inlet_data = INLET_VAR.name in df.columns
    nrows = 3 if has_inlet_data else 2

    fig, axs = plt.subplots(nrows=nrows, ncols=1, sharex=True)

    if has_inlet_data:
        axs[0].plot(t, df[INLET_VAR.name])
        axs[0].set_ylabel("Inlet Value")
        axs[0].grid(True)

    l1, = axs[-2].plot(t, df["targetValue"], "k--")
    axs[-2].plot(t, df["sensorValue"])
    axs[-2].set_ylabel("Process Variable")
    axs[-2].grid(True)
    axs[-2].legend((l1, ), ("Set Point", ))

    axs[-1].plot(t, df["outputSignal"])
    axs[-1].set_ylabel("Controller Output")
    axs[-1].grid(True)
    axs[-1].set_ylim(-0.1, 1.1)

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