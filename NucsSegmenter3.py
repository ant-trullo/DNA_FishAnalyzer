"""This function segments nuclei with otsu segmentation and cleans the results.

"""

import multiprocessing
import numpy as np
from scipy.ndimage.morphology import binary_fill_holes
from scipy import ndimage as ndi
from skimage.morphology import remove_small_objects, label, disk, remove_small_holes
from skimage.filters import median, threshold_local   # threshold_otsu
from skimage.measure import regionprops_table
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from PyQt5 import QtWidgets


def mmddnn(nnn):
    """Function defined just to feed the multiprocessing"""
    return median(nnn[0], nnn[1])


class NucsSegmenter3:
    """Class for nuclei detection"""
    def __init__(self, nucs_m, z_st, z_end):

        zlen        =  nucs_m.shape[0]
        block_size  =  int(nucs_m.shape[1] / 5)  # define a block_size, 5 block is generally good
        nucs_sgm    =  np.zeros(nucs_m.shape, dtype=np.uint8)

        if block_size // 2 - block_size / 2 == 0:  # check block_size is odd
            block_size += 1

        val  =  np.zeros(nucs_sgm.shape)             # initialization of threshold matrix (adaptive)
        for zz in range(z_st, z_end + 1):
            val[zz]  =  threshold_local(nucs_m[zz], block_size=block_size)  # threshod for each frame

        pbar  =  ProgressBar(total1=zlen)
        pbar.show()
        pbar.update_progressbar(0)

        nucs_sgm[z_st:z_end]  =  nucs_m[z_st:z_end] > val[z_st:z_end]     # thresholding

        for z in range(zlen):
            pbar.update_progressbar(z)
            nucs_sgm[z]  =  binary_fill_holes(nucs_sgm[z])

        self.nucs_sgm  =  nucs_sgm
        self.vals      =  val
        self.val_step  =  (nucs_m.max() - val.max()) / 50   # step thyreshols output


class NucsSegmenter3Pre:
    """Class for the nuclei raw data pre-processing"""
    def __init__(self, nucs):

        n_jobs    =  int((nucs.shape[0] - 4) / 3)
        job_args  =  list()
        job_args.append([nucs[:5], np.ones((3, 15, 15))])  # organize data for multiprocessing
        for jj in range(1, n_jobs):
            job_args.append([nucs[3 * jj:3 * (jj + 1) + 2], np.ones((3, 15, 15))])
        job_args.append([nucs[n_jobs * 3:], np.ones((3, 15, 15))])

        pool     =  multiprocessing.Pool()
        results  =  pool.map(mmddnn, job_args)
        pool.close()

        nucs_m      =  np.zeros_like(nucs)   # concatenate results
        nucs_m[:4]  =  results[0][:-1]
        for jj in range(1, n_jobs):
            nucs_m[3 * jj + 1:3 * (jj + 1) + 1]  =  results[jj][1:-1]
        nucs_m[n_jobs * 3 + 1:]  =  results[-1][1:]

        nucs_m    =  np.log(nucs_m)

        self.nucs_m  =  nucs_m


class NucsSegmenter3Split:
    """Class for nuclei segmentation (frame by frame in 2D)"""
    def __init__(self, nucs_sgm, min_dist):

        zlen  =  nucs_sgm.shape[0]
        nucs  =  np.zeros(nucs_sgm.shape, dtype=np.uint32)
        for z in range(zlen):        # frame by frame
            nucs[z]  =  label(nucs_sgm[z], connectivity=1)   # label image
            nucs[z]  =  remove_small_objects(nucs[z], 500)   # remove small objects

        pbar  =  ProgressBar(total1=zlen)
        pbar.show()
        pbar.update_progressbar(0)

        for z in range(zlen):    # frame by frame
            pbar.update_progressbar(z)
            rgp_bff  =  regionprops_table(nucs[z], properties=["label", "area"])   # regionprops of the labeled objects in each frame
            ws_idxs  =  np.where(rgp_bff["area"] > 15000)[0]                       # find the objects to treat with watershed
            for jj in ws_idxs:
                img                    =  (nucs[z] == rgp_bff["label"][jj])         # relabel
                distance               =  ndi.distance_transform_edt(img)           # matrix distance
                coords                 =  peak_local_max(distance, min_distance=min_dist, labels=img)    # determine peacks of the distance matrix
                mask                   =  np.zeros(distance.shape, dtype=bool)
                mask[tuple(coords.T)]  =  True
                markers, _             =  ndi.label(mask)
                lbls                   =  watershed(-distance, markers, mask=img)                   # watershed
                nucs[z]                =  nucs[z] * (1 - np.sign(lbls)) + (nucs[z].max() + 1) * np.sign(lbls) + lbls # organize labels

        pbar.close()

        self.nucs  =  nucs


class NucsSegmenter3Smooth:
    """Smooth nuclei detection before segmentation"""
    def __init__(self, nucs_2d_sgm, z_st, z_end):

        pbar  =  ProgressBar(total1=z_end-z_st)
        pbar.show()
        pbar.update_progressbar(0)

        nucs_smooth  =  np.zeros_like(nucs_2d_sgm)
        for zz in range(z_st, z_end):   # frame by frame
            # print(zz)
            pbar.update_progressbar(zz - z_st)
            nucs_smooth[zz]  =  median(nucs_2d_sgm[zz], disk(15))    # work with median filter
            nucs_smooth[zz]  =  remove_small_holes(nucs_smooth[zz], area_threshold=200)   # remove small objects

        pbar.close()
        self.nucs_smooth  =  nucs_smooth


class ProgressBar(QtWidgets.QWidget):
    """Simple progress bar widget"""
    def __init__(self, parent=None, total1=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar1  =  QtWidgets.QProgressBar()
        self.progressbar1.setMinimum(1)
        self.progressbar1.setMaximum(total1)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar1, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        """progress bar updater"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()
