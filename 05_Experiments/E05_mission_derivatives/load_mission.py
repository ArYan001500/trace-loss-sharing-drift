"""Load the mission-profile electrical log (1 Hz) and thermocouple log, average the
thermal rows within each recorded second, and join on absolute clock time with an
optional integer lag (seconds) applied to the electrical log."""
import os
import numpy as np, pandas as pd

def load(root, lag=0):
    e = pd.read_csv(os.path.join(root, "Mission profile", "fr06final.csv"), sep=";")
    e.columns = [c.replace("0x301:", "").strip() for c in e.columns]
    tcol = e.columns[0]
    e["time"] = pd.to_datetime("2022-05-17 " + e[tcol].astype(str)) + pd.to_timedelta(lag, unit="s")
    th = pd.read_excel(os.path.join(root, "Mission profile", "fr06final_temperature.xlsx"))
    th = th.rename(columns={th.columns[0]: "time"})
    th.columns = ["time"] + [str(c).split(" - ")[-1].strip() for c in th.columns[1:]]
    th["time"] = pd.to_datetime(th["time"]).dt.floor("s")
    tm = th.groupby("time").mean(numeric_only=True)
    keep = ["Ia_ef", "Ib_ef", "Ic_ef", "u_dc", "torque", "rpm_mech", "fs_set"]
    keep = [k for k in keep if k in e.columns]
    df = e.set_index("time")[keep].join(tm, how="inner")
    return df

if __name__ == "__main__":
    import sys
    df = load(sys.argv[1])
    print(df.shape, df.index.min(), df.index.max())
    print(df.columns.tolist())
    print(df[["Ia_ef", "Ib_ef", "Ic_ef", "u_dc", "fs_set"]].describe().round(1).to_string())
