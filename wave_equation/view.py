import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation


def natural_sort_key(s):
    import re

    return [int(text) if text.isdigit() else text.lower() for text in re.split("(\d+)", s)]


csv_files = sorted(glob.glob("./output/*.csv"), key=natural_sort_key)

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection="3d")


def update(frame):
    ax.clear()
    df = pd.read_csv(csv_files[frame], header=None)
    df.columns = ["x", "y", "u"]
    ax.scatter(df["x"], df["y"], df["u"], s=2, alpha=0.7)

    ax.set_zlim(-12, 12)


ani = FuncAnimation(fig, update, frames=len(csv_files), interval=900)

ani.save("wave.gif", writer="pillow", fps=10)
