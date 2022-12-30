#!/usr/bin/env python3

import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT = 'DejaVuSans.ttf'
#FONT = 'Arial.TTF'

font = ImageFont.truetype(FONT, 30)

# Printable characters.
chars = ''.join(chr(ic) for ic in range(0x20, 0x7f))

sizes = dict()

nx = 0
ny = 0
for c in chars:
    bbox = font.getbbox(c)
    cx, cy = int(np.ceil(bbox[2])), int(np.ceil(bbox[3]))
    sizes[c] = (cx, cy)
    nx += cx
    ny = max(ny, cy)

img = Image.new('L', (nx, ny))
draw = ImageDraw.Draw(img)

dx = 0
for c in chars:
    draw.text((dx, 0), c, font=font, fill='white', spacing=5)
    dx += sizes[c][0]

img.save('font.png', "PNG")
with open('font.bin', 'wb') as f:
    f.write(img.tobytes())

dx = 0
with open("font.geom", 'w') as f:
    f.write("{:} {:}\n".format(nx, ny))
    for ic in range(256):
        c = chr(ic)
        if c in sizes:
            cx, _ = sizes[c]
            f.write("{:} {:}\n".format(dx, cx))
            dx += cx
        else:
            f.write("{:} {:}\n".format(0, 0))
