"""This function removes the nuclei regions on the border."""


import numpy as np


class RemoveBorderNuclei:
    def __init__(self, nucs_dil):

        mask                 =  np.zeros_like(nucs_dil)    # define mask array
        mask[:, 1:-1, 1:-1]  =  1                          # put a border of 1 on the external surface

        idxs2rem  =  np.unique((1 - mask) * nucs_dil)[1:]  # find the tag of the nuclei touching the border that must be removed
        for k in idxs2rem:
            nucs_dil  *=  (1 - (nucs_dil == k)).astype(np.uint16)    # remove all these nuclei one by one

        self.nucs_dil  =  nucs_dil
