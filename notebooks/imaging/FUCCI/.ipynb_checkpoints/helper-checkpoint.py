
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np

import btrack
# from btrack import datasets

from skimage.io import imread
from skimage.util import montage
from skimage.transform import resize

from matplotlib.animation import FuncAnimation
from IPython.display import HTML


def plot_track_montage(
    image_stack,
    track_df,
    track_id,
    crop_size=200,
    ncols=8,
    id_col="ID",
    time_col="t",
    x_col="x",
    y_col="y",
    percentile_range=(1, 99.5),
    show=True,
):

    if image_stack.ndim != 3:
        raise ValueError(
            f"Expected image_stack shape (T, Y, X), "
            f"received {image_stack.shape}"
        )

    trajectory = (
        track_df.loc[track_df[id_col] == track_id]
        .sort_values(time_col)
        .reset_index(drop=True)
    )

    if trajectory.empty:
        raise ValueError(f"Track ID {track_id} was not found.")

    crops = []
    times = []

    half = crop_size // 2
    n_frames, height, width = image_stack.shape

    for _, row in trajectory.iterrows():
        t = int(round(row[time_col]))
        x = int(round(row[x_col]))
        y = int(round(row[y_col]))

        if not 0 <= t < n_frames:
            print(f"Skipping invalid time point t={t}")
            continue

        # Desired crop boundaries
        x0 = x - half
        x1 = x0 + crop_size
        y0 = y - half
        y1 = y0 + crop_size

        # Clip to image boundaries
        cx0 = max(0, x0)
        cx1 = min(width, x1)
        cy0 = max(0, y0)
        cy1 = min(height, y1)

        crop = image_stack[t, cy0:cy1, cx0:cx1]

        # Pad cells near the image boundary
        pad_left = cx0 - x0
        pad_right = x1 - cx1
        pad_top = cy0 - y0
        pad_bottom = y1 - cy1

        crop = np.pad(
            crop,
            (
                (pad_top, pad_bottom),
                (pad_left, pad_right),
            ),
            mode="constant",
            constant_values=0,
        )

        crops.append(crop)
        times.append(t)

    if not crops:
        raise ValueError(f"No valid frames found for track {track_id}.")

    crops = np.stack(crops)

    # Use common contrast limits for every panel
    vmin, vmax = np.percentile(crops, percentile_range)

    nrows = math.ceil(len(crops) / ncols)

    
    if show:
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(2.5 * ncols, 2.7 * nrows),
            squeeze=False,
        )
    
        

        for i, ax in enumerate(axes.flat):
            if i < len(crops):
                ax.imshow(
                    crops[i],
                    cmap="gray",
                    vmin=vmin,
                    vmax=vmax,
                )
    
                ax.scatter(
                    crop_size / 2,
                    crop_size / 2,
                    marker="+",
                    s=60,
                    linewidth=1.5,
                )
    
                ax.set_title(f"t = {times[i]}", fontsize=10)
            else:
                ax.set_visible(False)
    
            ax.axis("off")
    
        fig.suptitle(
            f"Track {track_id}: {len(crops)} time points",
            fontsize=16,
            y=1.01,
        )
    
        plt.tight_layout()

        plt.show()

    return crops, trajectory
