#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Configuration
# -------------------------------
mpicmd = ["mpirun", "-n", "10", "../qc_project"]
input_file = "input.dat"
data_file = "avcond.dat"

# -------------------------------
# Run the MPI program
# -------------------------------
print("Running MPI job...")
with open(input_file, "r") as fin:
    result = subprocess.run(
        mpicmd,
        stdin=fin,
        stdout=sys.stdout,
        stderr=sys.stderr,
        text=True
    )

if result.returncode != 0:
    print(f"MPI run failed with return code {result.returncode}.")
    sys.exit(1)

print("MPI job finished successfully.")

# -------------------------------
# Load data
# -------------------------------
data_path = Path(data_file)
if not data_path.exists():
    raise FileNotFoundError(f"Data file not found: {data_file}")

data = np.loadtxt(data_path)

if data.ndim == 1:
    data = data.reshape(1, -1)

if data.shape[1] < 4:
    raise ValueError("Expected at least 4 columns: x y1 y2 y3")

x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]

# -------------------------------
# Plot
# -------------------------------
plt.figure(figsize=(8, 5))

# y1 and y2 solid lines
line1, = plt.plot(x, y1, label=r"$\mathrm{Re}[\sigma(\omega)]$")
plt.plot(x, y2, label=r"$\omega\,\mathrm{Im}[\sigma(\omega)]$")

# y3 dashed, same color as y1
plt.plot(
    x, y3,
    linestyle="--",
    color=line1.get_color(),
    label=r"$\mathrm{Re}[\sigma_N(\omega)]$"
)

plt.xlabel(r"$\omega / 2\Delta_0$")
plt.ylabel("")

plt.legend()
plt.grid(True)

# Axis limits
plt.xlim(0.0, 6.0)
ymax = 1.2 * np.max(y1) if y1.size > 0 else 1.0
plt.ylim(0.0, ymax)

plt.tight_layout()
plt.show()
