from __future__ import division
import matplotlib

import trackpy.tracking as pt
import trackpy.identification as tid
import matplotlib.pyplot as plt
import random
import numpy as np
import string, Tkinter, tkFileDialog
from decimal import *
import Image
import numpy as np

root = Tkinter.Tk()
root.withdraw()
f_path,f_name = os.path.split(tkFileDialog.askopenfilename())
im= plt.imread(file_path)
'''X = range(10, 200, 20)
Y = range(10, 200, 20)
Z = range(100)
random.shuffle(X)
random.shuffle(Y)

img = tid.gen_fake_data(np.vstack([X, Y]), 5, 2.5, (210, 210))'''
bp_img = tid.band_pass(im, 2, 2.5)


res_lm = tid.find_local_max(bp_img, 3, .5)

locs, mass, r2 = tid.subpixel_centroid(bp_img, res_lm, 3)

# make figure
fig = plt.figure()
ax = fig.gca()
# display image
ax.imshow(bp_img, cmap='gray', interpolation='nearest')

# add pixval like output
ax.format_coord = lambda x, y: 'r=%d,c=%d,v=%0.2f' % (int(x + .5),
                                                      int(y + .5),
                                                      bp_img[int(x + .5), int(y + .5)] if
                                                      int(x + .5) < bp_img.shape[0] and int(y + .5) < bp_img.shape[1] else 0)

ax.plot(*locs, linestyle='none', marker='o')

plt.show()
