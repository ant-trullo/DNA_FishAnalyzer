"""This function piles up nuclei segmented in a z-stack.

Given as input the z-stack of nuclei segmented frame by frame,
this function gives as inpuut the same z-stack in which each
nucleus keeps the same tag z by z.
"""



import numpy as np
from skimage.measure import regionprops
from skimage.morphology import label
from PyQt5 import QtWidgets

import NucsPileUpUtility


class NucsPileUp:
    """Main function, does all the job"""
    def __init__(self, nucs_lbls_input):

        nucs_lbls        =  np.copy(nucs_lbls_input).astype(np.int)
        zlen             =  nucs_lbls.shape[0]                                      # number of z frames
        # nucs_lbls_piled  =  np.zeros(nucs_lbls.shape, dtype=np.int)                 # output matrix to fill
        nucs_lbls2pile   =  np.zeros(nucs_lbls.shape, dtype=np.int32)                 # matrix for correcction

        frames2work  =  list(np.where(np.sum(nucs_lbls, axis=(1, 2)) != 0)[0])       # z frames with nuclei to work on

        pbar  =  ProgressBar(total=len(frames2work))
        pbar.show()
        # pbar_idx  =  1

        [nucs_lbls, nucs_lbls_piled]  =  NucsPileUpUtility.NucsPileUpUtility(frames2work, nucs_lbls)

        rgp  =  regionprops(nucs_lbls_piled)                                            # regionprops of the reconstructed nuclei to check if the are some not assigned (single in tag in only 1 or 2 z frames)
        pbar.progressbar.setMaximum(len(rgp))
        pbar_idx2  =  0

        for ll in range(len(rgp)):
            pbar.update_progressbar(pbar_idx2 + 1)
            pbar_idx2     +=  1
#            print(np.diff(rgp[ll]['coords'][:, 0]).sum())
            if np.diff(rgp[ll]['coords'][:, 0]).sum() <= 1:                             # recognition of the bad reconstruction based on the z-coordinate of the 3D connected components
                nucs_bff          =  nucs_lbls_piled == rgp[ll]['label']
                nucs_lbls2pile   +=  nucs_bff                                           # if only 1 or 2 z frames are involved, the nucleus is added to the matrix for corrections
                nucs_lbls_piled  *=  1 - nucs_bff                                       # if only 1 or 2 z frames are involved, the nucleus is removed from the output 

        lbl_update  =  1                                                                # labeling of the matrix for correction: first label
        for k in range(zlen):
            nucs_lbls2pile[k]  =  label(nucs_lbls2pile[k])  + lbl_update * np.sign(nucs_lbls2pile[k])    # label is managed in order to not work in 3D, but frame by frame having new tags in each frames
            lbl_update         =  np.max((nucs_lbls2pile[k].max(), lbl_update))  + 1

        rgp2pile  =  regionprops(nucs_lbls2pile)                                                         # regionprops for the labeled matrix for corrections

        for j in range(len(rgp2pile)):
            z_coord     =  rgp2pile[j]['coords'][:, 0][0]                                                # for each components extract the z frame in which is present
            mask_bff_z  =  (nucs_lbls2pile[z_coord] == rgp2pile[j]['label'])                             # calc of the object
            if np.sum(mask_bff_z * nucs_lbls_piled[z_coord - 1]) > 0:                                    # check if there is overlapping in the previous frame of the output matrix
                idx                        =  np.unique(mask_bff_z * nucs_lbls_piled[z_coord - 1])[1]    # if there is overlapping, extract the proper tag from the output matrix
                nucs_lbls_piled[z_coord]  +=  idx * mask_bff_z                                           # add in the output matrix
                # print(np.min([z_coord + 1, nucs_lbls_piled.shape[0]]))
            elif np.sum(mask_bff_z * nucs_lbls_piled[np.min([z_coord + 1, nucs_lbls_piled.shape[0] - 1])]) > 0:      # if no overlapp in the previous frame, check in the following with the same protocol
#                 if z_coord + 1 <= nucs_lbls_piled.shape[0]:
                idx                        =  np.unique(mask_bff_z * nucs_lbls_piled[z_coord + 1])[1]
                nucs_lbls_piled[z_coord]  +=  idx * mask_bff_z

            pbar.close()

        rgp_final  =  regionprops(nucs_lbls_piled)                                      # small objects, coming from noise or small real object of the image, are removed
        volumes    =  np.zeros((2, len(rgp_final)))

        for rr in range(len(rgp_final)):
            volumes[:, rr]  =  rgp_final[rr]['label'], rgp_final[rr]['area']            # array with volume and label of each piled nucleus

        vol_thr  =  np.median(volumes[1, :]) / 3                                        # volume threshold (aribitrary)
        ii       =  np.where(volumes[1, :] < vol_thr)[0]                                # index of small objects

        for nn in volumes[0, ii]:
            nucs_lbls_piled  *=  (1 - (nucs_lbls_piled == nn))                          # removing small objects


        self.nucs_lbls_piled  =  nucs_lbls_piled.astype(np.int32)



class ProgressBar(QtWidgets.QWidget):
    """Simple progress bar widget"""
    def __init__(self, parent=None, total=20):
        super(ProgressBar, self).__init__(parent)
        self.name_line1  =  QtWidgets.QLineEdit()

        self.progressbar  =  QtWidgets.QProgressBar()
        self.progressbar.setMinimum(1)
        self.progressbar.setMaximum(total)

        main_layout  =  QtWidgets.QGridLayout()
        main_layout.addWidget(self.progressbar, 0, 0)

        self.setLayout(main_layout)
        self.setWindowTitle("Progress")
        self.setGeometry(500, 300, 300, 50)

    def update_progressbar(self, val1):
        """Progress bar updater"""
        self.progressbar.setValue(val1)
        QtWidgets.qApp.processEvents()







#class ProgressBar(QtGui.QWidget):
#    """Double progress bar widget"""
#    def __init__(self, parent=None, total1=20, total2=20):
#        super(ProgressBar, self).__init__(parent)
#        self.name_line1  =  QtGui.QLineEdit()
#
#        self.progressbar1  =  QtWidgets.QProgressBar()
#        self.progressbar1.setMinimum(1)
#        self.progressbar1.setMaximum(total1)
#
#        self.progressbar2  =  QtWidgets.QProgressBar()
#        self.progressbar2.setMinimum(1)
#        self.progressbar2.setMaximum(total2)
#
#        main_layout  =  QtGui.QGridLayout()
#        main_layout.addWidget(self.progressbar1, 0, 0)
#        main_layout.addWidget(self.progressbar2, 1, 0)
#
#        self.setLayout(main_layout)
#        self.setWindowTitle("Progress")
#        self.setGeometry(500, 300, 300, 50)
#
#    def update_progressbar1(self, val1):
#        self.progressbar1.setValue(val1)
#        QtWidgets.qApp.processEvents()
#
#    def update_progressbar2(self, val2):
#        self.progressbar2.setValue(val2)
#        QtWidgets.qApp.processEvents()
