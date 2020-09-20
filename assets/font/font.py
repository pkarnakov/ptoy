#!/usr/bin/env python3

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import io
import numpy as np

def HideAxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)


def DrawToArray(fig, **kwargs):
    buf = io.BytesIO()
    fig.savefig(buf, format='raw', **kwargs)
    buf.seek(0)
    img = np.reshape(np.frombuffer(buf.getvalue(), dtype=np.uint8),
                     newshape=(int(fig.bbox.height), int(fig.bbox.width), -1))
    buf.close()
    return img

dpi = 200
fontsize_px = 32
fontsize_in = fontsize_px / dpi
fontsize_pt = fontsize_in * 72.

figsize = [20, 1]
fig = plt.figure(dpi=dpi, figsize=figsize)

ax = plt.Axes(fig, [0., 0.] + figsize)
fig.add_axes(ax)
ax.set_xlim(0, figsize[0])
ax.set_ylim(0, figsize[1])
HideAxis(ax)

chars = ''.join(chr(i) for i in range(0x20, 0x7f))
print("{:} characters".format(len(chars)))

text = ax.text(0, 0, chars, fontfamily="Dejavu Sans Mono", fontsize=fontsize_pt)

bbox = text.get_tightbbox(fig.canvas.get_renderer())
bboxi = mpl.transforms.Bbox.from_bounds(bbox.x0 / dpi, bbox.y0 / dpi,
                                        bbox.width / dpi, bbox.height / dpi)
print(bbox.width, bbox.height)
print("character size: ", bbox.width / len(chars), bbox.height)

fig.bbox = bbox
fig.savefig("font.png", bbox_inches=bboxi)
img = DrawToArray(fig, bbox_inches=bboxi)
print(img.shape)

img.tofile("font.data")

with open("font.size", 'w') as f:
    # size_x size_y stride_x
    f.write("{:} {:} {:}".format(img.shape[1], img.shape[0],
                                 img.shape[1] / len(chars)))
