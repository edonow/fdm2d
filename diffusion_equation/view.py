import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from PIL import Image


def natural_sort_key(s):
    import re

    return [int(text) if text.isdigit() else text.lower() for text in re.split("(\d+)", s)]


csv_files = sorted(glob.glob("./output/*.csv"), key=natural_sort_key)

frames = []

for csv_file in csv_files:
    df = pd.read_csv(csv_file, header=None)
    df.columns = ["x", "y", "u"]

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(df["x"], df["y"], 0, c=df["u"], cmap=plt.get_cmap("jet"), s=2, alpha=0.7)
    ax.set_zlim(-12, 12)

    plt.tight_layout()
    image_path = csv_file.replace(".csv", ".png")
    plt.savefig(image_path)
    plt.close()

    frames.append(Image.open(image_path))

frames[0].save("animation.gif", save_all=True, append_images=frames[1:], loop=0, duration=100)

for file_path in glob.glob("./output/*.png"):
    os.remove(file_path)
