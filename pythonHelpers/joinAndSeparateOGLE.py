import pandas as pd
import os
import numpy as np

local = os.path.split(os.path.abspath(__file__))[0]
files = [f for f in os.listdir(os.path.join(local, "OGLE-III")) if f.endswith(".dat")]

for file in files:
    a = pd.read_csv(os.path.join(local, "OGLE-III", file), names = ["HJD","MAG","ERR"], sep = r"\s+")
    a["SRV"] = "OGLE-III"
    a["MAG"] -= np.median(a["MAG"].values)
    b = pd.read_csv(os.path.join(local, "OGLE-IV", file), names = ["HJD","MAG","ERR"], sep = r"\s+")
    b["SRV"] = "OGLE-IV"
    b["MAG"] -= np.median(b["MAG"].values)
    c = pd.concat([a, b])
    c.to_csv(os.path.join(local, "testFiles", file), index = False)
