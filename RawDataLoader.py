"""This function loads raw data.

The channel are associated following the standard of the acquisition of Maelle.
"""


import numpy as np
import czifile


class RawDataLoader:
    def __init__(self, raw_data_fname):

        raw_data  =  np.squeeze(czifile.imread(raw_data_fname))

        a  =  czifile.CziFile(raw_data_fname)  # read info about pixel size
        b  =  a.metadata()

        start      =  b.find("ScalingZ")
        end        =  b[start + 9:].find("ScalingZ")
        pix_sizeZ  =  float(b[start + 9:start + 7 + end]) * 1000000

        start      =  b.find("ScalingX")
        end        =  b[start + 9:].find("ScalingX")
        pix_sizeX  =  float(b[start + 9:start + 7 + end]) * 1000000

        self.ch1        =  raw_data[0]
        self.ch2        =  raw_data[1]
        self.ch3        =  raw_data[2]
        self.dapi       =  raw_data[3]
        self.pix_sizeX  =  pix_sizeX
        self.pix_sizeZ  =  pix_sizeZ
