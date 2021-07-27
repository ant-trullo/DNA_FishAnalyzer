"""This function crosses information of the several channels.

After loading the segmented spots and nuclei, it calculates nucleus per nucleus
the matrix of the distances between all the spots of the three chennels; we can have
a maximum of 2 dots per channel, and all the six possible distances are calculated
in order the minimu possible distance for un association, not for the sum of the distances.
Ex: c1, c2, e1, e2 are goi,ng to be associated like c1-e1, c2-e2 if d(c1-e1) (or d(c2-e2)) is
the shortest possbile distance, we don't care if at the end
d(c1-e1) + d(c2-e2) > d(c1-e2) + d(c2-e1)
"""

import numpy as np
from skimage.measure import label, regionprops_table
import xlsxwriter
# import pyqtgraph as pg
from PyQt5 import QtWidgets


def dist_choose(x1, x2, pix_sizeZ, pix_sizeX):
    """Calculates the distances between to 3D points: if one is (0,0,0), gives 1000000000"""
    if x1.sum() == 0 or x2.sum() == 0:
        dist  =  1000000000
    else:
        dist  =  np.sqrt(((x1[0] - x2[0]) * pix_sizeZ) ** 2 + ((x1[1] - x2[1]) * pix_sizeX) ** 2 + ((x1[2] - x2[2]) * pix_sizeX) ** 2)
    return dist


class CrossInfo:
    """Only class, does all the job"""
    def __init__(self, analysis_folder, spts_chs_modif):

        spts_ch1              =  spts_chs_modif[:, :, :, 0]                                     # load analysis results
        spts_ch2              =  spts_chs_modif[:, :, :, 1]
        spts_ch3              =  spts_chs_modif[:, :, :, 2]
        nucs_dapi             =  np.load(analysis_folder + '/nucs_dapi.npy')
        pix_sizeZ, pix_sizeX  =  np.load(analysis_folder + '/pix_sizes.npy')
        try:
            spts_clusters_flag  =  np.load(analysis_folder + '/spts_clstrs_flag.npy')
        except FileNotFoundError:
            spts_clusters_flag  =  1

        spts_ch1_lbls  =  label(spts_ch1, connectivity=1)                                           # spots labeling
        spts_ch2_lbls  =  label(spts_ch2, connectivity=1)
        spts_ch3_lbls  =  label(spts_ch3, connectivity=1)

        rgp_ch1  =  regionprops_table(spts_ch1_lbls, properties=["label", "coords", "area"])        # dictionary with useful info of spots
        rgp_ch2  =  regionprops_table(spts_ch2_lbls, properties=["label", "coords", "area"])
        rgp_ch3  =  regionprops_table(spts_ch3_lbls, properties=["label", "coords", "area"])

        idxs2rm_ch1  =  np.where(rgp_ch1["area"] < 4)[0]                                            # check id of spots smaller than 4 pixels to delete them
        idxs2rm_ch2  =  np.where(rgp_ch2["area"] < 4)[0]
        idxs2rm_ch3  =  np.where(rgp_ch3["area"] < 4)[0]

        for kk in idxs2rm_ch1:
            spts_ch1[rgp_ch1["coords"][kk][:, 0], rgp_ch1["coords"][kk][:, 1], rgp_ch1["coords"][kk][:, 2]]  =  0           # remove small spots using the coordinates of the pixels of the spots itself

        for kk in idxs2rm_ch2:
            spts_ch2[rgp_ch2["coords"][kk][:, 0], rgp_ch2["coords"][kk][:, 1], rgp_ch2["coords"][kk][:, 2]]  =  0

        for kk in idxs2rm_ch3:
            spts_ch3[rgp_ch3["coords"][kk][:, 0], rgp_ch3["coords"][kk][:, 1], rgp_ch3["coords"][kk][:, 2]]  =  0

        spts_ch1_lbls  =  label(spts_ch1, connectivity=1)                                           # relabel spots after removal of small objects
        spts_ch2_lbls  =  label(spts_ch2, connectivity=1)
        spts_ch3_lbls  =  label(spts_ch3, connectivity=1)

        rgp_ch1  =  regionprops_table(spts_ch1_lbls, properties=["label", "centroid", "area", "coords"])      # new dictionary of the selected spots
        rgp_ch2  =  regionprops_table(spts_ch2_lbls, properties=["label", "centroid", "area", "coords"])
        rgp_ch3  =  regionprops_table(spts_ch3_lbls, properties=["label", "centroid", "area", "coords"])

        idxs_nucs  =  np.unique(nucs_dapi)[1:]                                                              # tags of the nuclei
        dists_mtx  =  np.zeros((6, idxs_nucs.size))                                                         # matrix of the distances: for each nucleus you have at maximum 6 distances (c1-e1, c2-e2, c1-s1, c2-s2, e1-, s1, e2-s2)

        spts_ch1_fin  =  np.zeros(spts_ch1.shape, np.uint8)
        spts_ch2_fin  =  np.zeros(spts_ch2.shape, np.uint8)
        spts_ch3_fin  =  np.zeros(spts_ch3.shape, np.uint8)

        pbar  =  ProgressBar(total1=idxs_nucs.size)
        pbar.show()
        pbar.update_progressbar1(0)

            # print(counts)
        for counts, jj in enumerate(idxs_nucs):
            pbar.update_progressbar1(counts)
            sing_nuc           =  (nucs_dapi == jj)                                                         # for each nucleus
            # CH1
            spts_ch1_sing_tags  =  spts_ch1_lbls * sing_nuc                                              # find the tag of the ch1 spots inside the single nucleus
            spts_ch1_sing_tags  =  spts_ch1_sing_tags[spts_ch1_sing_tags != 0]
            spts_ch1_sing_tags  =  np.unique(spts_ch1_sing_tags)
            spts_ch1_sing_idxs  =  list()
            for gg in spts_ch1_sing_tags:
                spts_ch1_sing_idxs.append(np.where(rgp_ch1["label"] == gg)[0][0])

            ctrs_ch1  =  np.zeros((2, 3))

            if len(spts_ch1_sing_idxs) >= 2:
                area_idxs           =  np.zeros((2, len(spts_ch1_sing_idxs)))
                area_idxs[0, :]     =  np.take(rgp_ch1["area"], spts_ch1_sing_idxs)
                area_idxs[1, :]     =  np.asarray(spts_ch1_sing_idxs)
                area_idxs           =  (area_idxs[:, area_idxs[0].argsort()][1, -2:]).astype(np.int64)
                ctrs_ch1[0, :]      =  np.array([rgp_ch1["centroid-0"][area_idxs[0]], rgp_ch1["centroid-1"][area_idxs[0]], rgp_ch1["centroid-2"][area_idxs[0]]])
                ctrs_ch1[1, :]      =  np.array([rgp_ch1["centroid-0"][area_idxs[1]], rgp_ch1["centroid-1"][area_idxs[1]], rgp_ch1["centroid-2"][area_idxs[1]]])
                spts_ch1_fin[rgp_ch1["coords"][area_idxs[0]][:, 0], rgp_ch1["coords"][area_idxs[0]][:, 1], rgp_ch1["coords"][area_idxs[0]][:, 2]]  =  1
                spts_ch1_fin[rgp_ch1["coords"][area_idxs[1]][:, 0], rgp_ch1["coords"][area_idxs[1]][:, 1], rgp_ch1["coords"][area_idxs[1]][:, 2]]  =  1

            if len(spts_ch1_sing_idxs) == 1:
                ctrs_ch1[0, :]  =  np.array([rgp_ch1["centroid-0"][spts_ch1_sing_idxs[0]], rgp_ch1["centroid-1"][spts_ch1_sing_idxs[0]], rgp_ch1["centroid-2"][spts_ch1_sing_idxs[0]]])
                spts_ch1_fin[rgp_ch1["coords"][spts_ch1_sing_idxs[0]][:, 0], rgp_ch1["coords"][spts_ch1_sing_idxs[0]][:, 1], rgp_ch1["coords"][spts_ch1_sing_idxs[0]][:, 2]]  =  1

            # CH2
            spts_ch2_sing_tags  =  spts_ch2_lbls * sing_nuc                                              # find the tag of the ch1 spots inside the single nucleus
            spts_ch2_sing_tags  =  spts_ch2_sing_tags[spts_ch2_sing_tags != 0]
            spts_ch2_sing_tags  =  np.unique(spts_ch2_sing_tags)
            spts_ch2_sing_idxs  =  list()
            for gg in spts_ch2_sing_tags:
                spts_ch2_sing_idxs.append(np.where(rgp_ch2["label"] == gg)[0][0])

            ctrs_ch2  =  np.zeros((2, 3))

            if len(spts_ch2_sing_idxs) >= 2:
                area_idxs        =  np.zeros((2, len(spts_ch2_sing_idxs)))
                area_idxs[0, :]  =  np.take(rgp_ch2["area"], spts_ch2_sing_idxs)
                area_idxs[1, :]  =  np.asarray(spts_ch2_sing_idxs)
                area_idxs        =  (area_idxs[:, area_idxs[0].argsort()][1, -2:]).astype(np.int64)
                ctrs_ch2[0, :]   =  np.array([rgp_ch2["centroid-0"][area_idxs[0]], rgp_ch2["centroid-1"][area_idxs[0]], rgp_ch2["centroid-2"][area_idxs[0]]])
                ctrs_ch2[1, :]   =  np.array([rgp_ch2["centroid-0"][area_idxs[1]], rgp_ch2["centroid-1"][area_idxs[1]], rgp_ch2["centroid-2"][area_idxs[1]]])
                spts_ch2_fin[rgp_ch2["coords"][area_idxs[0]][:, 0], rgp_ch2["coords"][area_idxs[0]][:, 1], rgp_ch2["coords"][area_idxs[0]][:, 2]]  =  1
                spts_ch2_fin[rgp_ch2["coords"][area_idxs[1]][:, 0], rgp_ch2["coords"][area_idxs[1]][:, 1], rgp_ch2["coords"][area_idxs[1]][:, 2]]  =  1

            if len(spts_ch2_sing_idxs) == 1:
                ctrs_ch2[0, :]  =  np.array([rgp_ch2["centroid-0"][spts_ch2_sing_idxs[0]], rgp_ch2["centroid-1"][spts_ch2_sing_idxs[0]], rgp_ch2["centroid-2"][spts_ch2_sing_idxs[0]]])
                spts_ch2_fin[rgp_ch2["coords"][spts_ch2_sing_idxs[0]][:, 0], rgp_ch2["coords"][spts_ch2_sing_idxs[0]][:, 1], rgp_ch2["coords"][spts_ch2_sing_idxs[0]][:, 2]]  =  1

            # CH3
            spts_ch3_sing_tags  =  spts_ch3_lbls * sing_nuc                                              # find the tag of the ch1 spots inside the single nucleus: mask the nucleus on the labeled spots
            spts_ch3_sing_tags  =  spts_ch3_sing_tags[spts_ch3_sing_tags != 0]
            spts_ch3_sing_tags  =  np.unique(spts_ch3_sing_tags)                                         # find the surviving tags
            spts_ch3_sing_idxs  =  list()
            for gg in spts_ch3_sing_tags:
                spts_ch3_sing_idxs.append(np.where(rgp_ch3["label"] == gg)[0][0])                        # find the indexes in the dictionary of the surviving tags

            ctrs_ch3  =  np.zeros((2, 3))                                                                  # initialize the matrix of centroids

            if len(spts_ch3_sing_idxs) >= 2:                                                               # if we have 2 spots or more in the mask, search the biggest 2 (it's redundant in case of only 2 spots, but works anyway reducing code complexity)
                area_idxs        =  np.zeros((2, len(spts_ch3_sing_idxs)))                              # initialize with volume and relative indexes
                area_idxs[0, :]  =  np.take(rgp_ch3["area"], spts_ch3_sing_idxs)
                area_idxs[1, :]  =  np.asarray(spts_ch3_sing_idxs)
                area_idxs        =  (area_idxs[:, area_idxs[0].argsort()][1, -2:]).astype(np.int64)        # sort the matrix with respect to the volume
                ctrs_ch3[0, :]   =  np.array([rgp_ch3["centroid-0"][area_idxs[0]], rgp_ch3["centroid-1"][area_idxs[0]], rgp_ch3["centroid-2"][area_idxs[0]]])   # record the coordinates of the centroids of the biggest 2
                ctrs_ch3[1, :]   =  np.array([rgp_ch3["centroid-0"][area_idxs[1]], rgp_ch3["centroid-1"][area_idxs[1]], rgp_ch3["centroid-2"][area_idxs[1]]])
                spts_ch3_fin[rgp_ch3["coords"][area_idxs[0]][:, 0], rgp_ch3["coords"][area_idxs[0]][:, 1], rgp_ch3["coords"][area_idxs[0]][:, 2]]  =  1
                spts_ch3_fin[rgp_ch3["coords"][area_idxs[1]][:, 0], rgp_ch3["coords"][area_idxs[1]][:, 1], rgp_ch3["coords"][area_idxs[1]][:, 2]]  =  1

            if len(spts_ch3_sing_idxs) == 1:                                                              # if there is only 1 spot, no sorting is needed
                ctrs_ch3[0, :]  =  np.array([rgp_ch3["centroid-0"][spts_ch3_sing_idxs[0]], rgp_ch3["centroid-1"][spts_ch3_sing_idxs[0]], rgp_ch3["centroid-2"][spts_ch3_sing_idxs[0]]])
                spts_ch3_fin[rgp_ch3["coords"][spts_ch3_sing_idxs[0]][:, 0], rgp_ch3["coords"][spts_ch3_sing_idxs[0]][:, 1], rgp_ch3["coords"][spts_ch3_sing_idxs[0]][:, 2]]  =  1

            # FILLING DISTANCE MATRIX

            all_d_ce  =  [dist_choose(ctrs_ch1[0, :], ctrs_ch2[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[0, :], ctrs_ch2[1, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[1, :], ctrs_ch2[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[1, :], ctrs_ch2[1, :], pix_sizeZ, pix_sizeX)]

            if np.argmin(all_d_ce) == 1 or np.argmin(all_d_ce) == 2:
                dists_mtx[:2, counts]  =  all_d_ce[1], all_d_ce[2]
            else:
                dists_mtx[:2, counts]  =  all_d_ce[0], all_d_ce[3]

            all_d_cs  =  [dist_choose(ctrs_ch1[0, :], ctrs_ch3[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[0, :], ctrs_ch3[1, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[1, :], ctrs_ch3[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch1[1, :], ctrs_ch3[1, :], pix_sizeZ, pix_sizeX)]

            if np.argmin(all_d_cs) == 1 or np.argmin(all_d_cs) == 2:
                dists_mtx[2:4, counts]  =  all_d_cs[1], all_d_cs[2]
            else:
                dists_mtx[2:4, counts]  =  all_d_cs[0], all_d_cs[3]

            all_d_es  =  [dist_choose(ctrs_ch2[0, :], ctrs_ch3[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch2[0, :], ctrs_ch3[1, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch2[1, :], ctrs_ch3[0, :], pix_sizeZ, pix_sizeX), dist_choose(ctrs_ch2[1, :], ctrs_ch3[1, :], pix_sizeZ, pix_sizeX)]

            if np.argmin(all_d_es) == 1 or np.argmin(all_d_es) == 2:
                dists_mtx[4:, counts]  =  all_d_es[1], all_d_es[2]
            else:
                dists_mtx[4:, counts]  =  all_d_es[0], all_d_es[3]

        pbar.close()

        workbook  =  xlsxwriter.Workbook(analysis_folder + "/distances_modified.xlsx")
        sheet1    =  workbook.add_worksheet("")

        sheet1.write(0, 0, "Nucs Id")
        sheet1.write(0, 1, "Dist ch1_1-ch2_1")
        sheet1.write(0, 2, "Dist ch1_2-ch2_2")
        sheet1.write(0, 3, "Dist ch1_1-ch3_1")
        sheet1.write(0, 4, "Dist ch1_2-ch3_2")
        sheet1.write(0, 5, "Dist ch2_1-ch3_1")
        sheet1.write(0, 6, "Dist ch2_2-ch3_2")

        for ll in range(dists_mtx.shape[1]):
            sheet1.write(ll + 1, 0, idxs_nucs[ll])
            if dists_mtx[0, ll] < 1000000000:
                sheet1.write(ll + 1, 1, dists_mtx[0, ll])
            else:
                sheet1.write(ll + 1, 1, "----")

            if dists_mtx[1, ll] < 1000000000:
                sheet1.write(ll + 1, 2, dists_mtx[1, ll])
            else:
                sheet1.write(ll + 1, 2, "----")
            if dists_mtx[2, ll] < 1000000000:
                sheet1.write(ll + 1, 3, dists_mtx[2, ll])
            else:
                sheet1.write(ll + 1, 3, "----")
            if dists_mtx[3, ll] < 1000000000:
                sheet1.write(ll + 1, 4, dists_mtx[3, ll])
            else:
                sheet1.write(ll + 1, 4, "----")
            if dists_mtx[4, ll] < 1000000000:
                sheet1.write(ll + 1, 5, dists_mtx[4, ll])
            else:
                sheet1.write(ll + 1, 5, "----")
            if dists_mtx[5, ll] < 1000000000:
                sheet1.write(ll + 1, 6, dists_mtx[5, ll])
            else:
                sheet1.write(ll + 1, 6, "----")

        if spts_clusters_flag[0] == 2:
            sheet2  =  workbook._add_sheet("Overlap")
            sheet2.write(0, 0, "CH1_Id")
            sheet2.write(0, 1, "Volume")
            sheet2.write(0, 2, "z centroid")
            sheet2.write(0, 3, "x centroid")
            sheet2.write(0, 4, "y centroid")
            sheet2.write(0, 5, "CH2_Id")
            sheet2.write(0, 6, "Volume")
            sheet2.write(0, 7, "z centroid")
            sheet2.write(0, 8, "x centroid")
            sheet2.write(0, 9, "y centroid")
            sheet2.write(0, 10, "CH3_Id")
            sheet2.write(0, 11, "Volume")
            sheet2.write(0, 12, "z centroid")
            sheet2.write(0, 13, "x centroid")
            sheet2.write(0, 14, "y centroid")

            sheet2.write(1, 16, "CH2 on Clstrs")
            sheet2.write(2, 16, "CH2 not on Clstrs")
            sheet2.write(3, 16, "CH3 on Clstrs")
            sheet2.write(4, 16, "CH3 not on Clstrs")

            sheet2.write(0, 17, "Numb")
            sheet2.write(0, 18, "%")

            numb_spts_ch2    =  np.unique(spts_ch2_lbls[spts_ch2_lbls != 0]).size
            numb_spts_ch3    =  np.unique(spts_ch3_lbls[spts_ch3_lbls != 0]).size
            ch2_on_clstr     =  spts_ch2_lbls * np.sign(spts_ch1_lbls)
            ch2_on_clstr     =  np.unique(ch2_on_clstr[ch2_on_clstr != 0]).size
            ch2_noton_clstr  =  numb_spts_ch2 - ch2_on_clstr
            ch3_on_clstr     =  spts_ch3_lbls * np.sign(spts_ch1_lbls)
            ch3_on_clstr     =  np.unique(ch3_on_clstr[ch3_on_clstr != 0]).size
            ch3_noton_clstr  =  numb_spts_ch3 - ch3_on_clstr

            sheet2.write(1, 17, np.int64(ch2_on_clstr))
            sheet2.write(2, 17, np.int64(ch2_noton_clstr))
            sheet2.write(3, 17, np.int64(ch3_on_clstr))
            sheet2.write(4, 17, np.int64(ch3_noton_clstr))
            sheet2.write(1, 18, np.float64(100 * ch2_on_clstr / numb_spts_ch2))
            sheet2.write(2, 18, np.float64(100 * ch2_noton_clstr / numb_spts_ch2))
            sheet2.write(3, 18, np.float64(100 * ch3_on_clstr / numb_spts_ch3))
            sheet2.write(4, 18, np.float64(100 * ch3_noton_clstr / numb_spts_ch3))

            row_idx  =  0
            for bb in range(len(rgp_ch1["label"])):
                sheet2.write(1 + row_idx, 0, "Clst_" + str(rgp_ch1["label"][bb]))
                sheet2.write(1 + row_idx, 1, rgp_ch1["area"][bb])
                sheet2.write(1 + row_idx, 2, rgp_ch1["centroid-0"][bb])
                sheet2.write(1 + row_idx, 3, rgp_ch1["centroid-1"][bb])
                sheet2.write(1 + row_idx, 4, rgp_ch1["centroid-2"][bb])

                bff_ch2  =  spts_ch2_lbls[rgp_ch1["coords"][bb][:, 0], rgp_ch1["coords"][bb][:, 1], rgp_ch1["coords"][bb][:, 2]]
                bff_ch2  =  np.unique(bff_ch2[bff_ch2 != 0])
                for dd in range(bff_ch2.size):
                    iidd  =  np.where(rgp_ch2["label"] == bff_ch2[dd])[0][0]
                    sheet2.write(1 + row_idx + dd, 5, "Spts_" + str(rgp_ch2["label"][iidd]))
                    sheet2.write(1 + row_idx + dd, 6, rgp_ch2["area"][iidd])
                    sheet2.write(1 + row_idx + dd, 7, rgp_ch2["centroid-0"][iidd])
                    sheet2.write(1 + row_idx + dd, 8, rgp_ch2["centroid-1"][iidd])
                    sheet2.write(1 + row_idx + dd, 9, rgp_ch2["centroid-2"][iidd])

                bff_ch3  =  spts_ch3_lbls[rgp_ch1["coords"][bb][:, 0], rgp_ch1["coords"][bb][:, 1], rgp_ch1["coords"][bb][:, 2]]
                bff_ch3  =  np.unique(bff_ch3[bff_ch3 != 0])
                for pp in range(bff_ch3.size):
                    ddii  =  np.where(rgp_ch3["label"] == bff_ch3[pp])[0][0]
                    sheet2.write(1 + row_idx + pp, 10, "Spts_" + str(rgp_ch3["label"][ddii]))
                    sheet2.write(1 + row_idx + pp, 11, rgp_ch3["area"][ddii])
                    sheet2.write(1 + row_idx + pp, 12, rgp_ch3["centroid-0"][ddii])
                    sheet2.write(1 + row_idx + pp, 13, rgp_ch3["centroid-1"][ddii])
                    sheet2.write(1 + row_idx + pp, 14, rgp_ch3["centroid-2"][ddii])

                row_idx  +=  np.max([bff_ch2.size, bff_ch3.size])

        workbook.close()

        mtx2show              =  np.zeros(spts_ch2.shape + (3,), dtype=np.uint8)
        mtx2show[:, :, :, 0]  =  125 * np.sign(nucs_dapi) * (1 - spts_ch1_fin) * (1 - spts_ch2_fin) * (1 - spts_ch3_fin)
        mtx2show[:, :, :, 1]  =  125 * np.sign(nucs_dapi) * (1 - spts_ch1_fin) * (1 - spts_ch2_fin) * (1 - spts_ch3_fin)
        mtx2show[:, :, :, 2]  =  125 * np.sign(nucs_dapi) * (1 - spts_ch1_fin) * (1 - spts_ch2_fin) * (1 - spts_ch3_fin)
        mtx2show[:, :, :, 1] +=  255 * spts_ch1_fin
        mtx2show[:, :, :, 0] +=  255 * spts_ch2_fin
        mtx2show[:, :, :, 2] +=  255 * spts_ch3_fin
        # pg.image(mtx2show)

        self.mtx2show   =  mtx2show
        self.nucs_dapi  =  nucs_dapi


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

    def update_progressbar1(self, val1):
        """First progress bar updater"""
        self.progressbar1.setValue(val1)
        QtWidgets.qApp.processEvents()




