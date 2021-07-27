"""This function segments nuclei with otsu segmentation and cleans the results.

"""

import numpy as np
from skimage.filters import threshold_otsu, gaussian
from scipy.ndimage.morphology import binary_fill_holes
from skimage.morphology import remove_small_objects, label


class NucsSegmenter2:
    def __init__(self, nucs, g_kern, thr_fctr):

        zlen    =  nucs.shape[0]
        nucs_f  =  gaussian(nucs, g_kern)
        val     =  threshold_otsu(nucs_f)
        nucs_t  =  (nucs_f > thr_fctr * val).astype(np.uint8)

        nucs_sgm  =  np.zeros(nucs_t.shape, np.uint32)

        for z in range(zlen):
            nucs_sgm[z]  =  binary_fill_holes(nucs_t[z])

        for z in range(zlen):
            nucs_sgm[z]  =  label(nucs_sgm[z], connectivity=1)

        for z in range(zlen):
            nucs_sgm[z]  =  remove_small_objects(nucs_sgm[z], 160)

        self.nucs_lbls  =  nucs_sgm
