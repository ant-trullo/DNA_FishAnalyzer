"""This function is used in the manual modification of the labels.

It starts from a matrix-image with labeled connected objects and cuts or merges
labels depending on the array_region and end_pts. Array_region is the vector
with the values of the pixels above the segment in the matrix-image. end_pts are
the coordinate of the extremal point of the segment. The out put is the updated
version of the labels_img.
"""


import numpy as np
from skimage.morphology import label
from skimage.morphology import remove_small_objects

import BresenhamLine


class LabelsModify:
    """Only class, does all the job"""
    def __init__(self, labels_img, end_pts):
        self.labels_img    =  labels_img
        self.end_pts       =  end_pts

        line_mask   =  np.zeros(self.labels_img.shape)
        pts_coords  =  BresenhamLine.BresenhamLine(self.end_pts[0, 0], self.end_pts[0, 1], self.end_pts[1, 0], self.end_pts[1, 1]).coords     # coordinate of the points overlapped by the ROI segment
        for k in range(len(pts_coords)):
            line_mask[pts_coords[k][0], pts_coords[k][1]] = 1                                                           # the segment is going to be 1 in the mask

        idxs  =  np.unique((line_mask * labels_img))[1:]                                                                # labels of the speckles overlapped by the segment

        """ If after this operations there is only one value left (derivatives is zero always), means that the segment overlaps
            only one speckle: it is going to cut. If there is more than one value, it has to join two or more speckles. """

        if idxs.size == 1:
            speckle     =  (labels_img == idxs).astype(np.int)                                             # mask of the speckle to work on
            labels_fin  =  self.labels_img * (1 - speckle)                                                 # cut

            speckle       *=  (1 - line_mask.astype(np.int))
            speckle_lbls   =  label(speckle, connectivity=1)                                              # assigns labels again
            speckle_lbls   =  remove_small_objects(speckle_lbls, 10)                                      # in case3 of strange shapes, you can have single pixels (or very small groups) that takes lower labels.
            speckle_lbls   =  label(speckle_lbls, connectivity=1)                                                          # With these two lines we get read of them
            labels_fin    +=  (speckle_lbls == 1).astype(np.int) * idxs[0].astype(np.int)
            labels_fin    +=  (speckle_lbls == 2).astype(np.int) * (self.labels_img.max() + 1)

        else:
            labels_mask  =  np.zeros(self.labels_img.shape)

            for l in idxs:
                labels_mask  +=  (self.labels_img == l).astype(np.int)                                                  # finds the mask with all the involved speckles

            labels_fin  =  self.labels_img * (1 - labels_mask) + idxs[0] * labels_mask                                  # assigns one single label to them

        self.labels_fin     =  labels_fin
