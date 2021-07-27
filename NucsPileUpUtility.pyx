import numpy as np
cimport numpy as np
#import pyximport
#pyximport.install(reload_support=True)

DTYPE = np.int
ctypedef np.int_t DTYPE_t

cpdef unique_cython_int_dim3(np.ndarray[DTYPE_t, ndim=3] a):
    cdef int z, x, y
    cdef int zlen = a.shape[0]
    cdef int xlen = a.shape[1]
    cdef int ylen = a.shape[2]
    cdef set s = set()
    cdef list s_list = []
    
    for z in range(zlen):
        for x in range(xlen):
            for y in range(ylen):
                s.add(a[z, x, y])
    for s_el in s:
        s_list.append(s_el)
    return s_list[1:]


cpdef unique_cython_int(np.ndarray[DTYPE_t, ndim=2] a):
    cdef int x, y
    cdef int xlen = a.shape[0]
    cdef int ylen = a.shape[1]
    cdef set s = set()
    cdef list s_list = []
    
    for x in range(xlen):
        for y in range(ylen):
            s.add(a[x, y])
    for s_el in s:
        s_list.append(s_el)
    return s_list[1:]


cpdef median_filter(np.ndarray[DTYPE_t, ndim=2] a):
    cdef list ll = []
    cdef int x, y
    cdef int xlen = a.shape[0]
    cdef int ylen = a.shape[1]
    for x in range(xlen):
        for y in range(ylen):
            if a[x, y] != 0:
                ll.append(a[x, y])
    ll.sort()                
    return ll[len(ll) // 2]



cpdef NucsPileUpUtility(list frames2work, np.ndarray[DTYPE_t, ndim=3] nucs_lbls):

    cdef np.ndarray[DTYPE_t, ndim=3] nucs_lbls_piled = np.zeros([nucs_lbls.shape[0], nucs_lbls.shape[1], nucs_lbls.shape[2]], dtype=np.int)
    cdef np.ndarray[DTYPE_t, ndim=2] mask_bff        = np.zeros([nucs_lbls.shape[1], nucs_lbls.shape[2]], dtype=np.int)
    cdef np.ndarray[DTYPE_t, ndim=2] mtx_ref         = np.zeros([nucs_lbls.shape[1], nucs_lbls.shape[2]], dtype=np.int)

    cdef int new_tag  =  1
    cdef int zlen     =  nucs_lbls.shape[0]
    cdef int xlen     =  nucs_lbls.shape[1]
    cdef int ylen     =  nucs_lbls.shape[2]
    cdef int f_len    =  len(frames2work)
    cdef int z, zz, x, y, j, idx_ref, idxs_lght, zf
    cdef list idxs, idxs2
    
    for z in range(f_len):
        zf         =  frames2work[z]
        idxs       =  unique_cython_int(nucs_lbls[zf, :, :])                              # list of nuclei tags in the considered frame
        idxs_lght  =  len(idxs)
        for j in range(idxs_lght):
            
            #print(j, idxs)
#            print(j)
            for x in range(xlen):
                for y in range(ylen):
                    if nucs_lbls[zf, x, y] == idxs[j]:
                        #print()
                        mask_bff[x, y] = 1
                        nucs_lbls_piled[zf, x, y] = new_tag
                        nucs_lbls[zf, x, y]       = 0
            #print(mask_bff.sum())
            for zz in range(zf + 1, zlen):
                idxs2  =  unique_cython_int(mask_bff * nucs_lbls[zz])
                #print(mask_bff.sum())

                if len(idxs2) > 0:
                    for x in range(xlen):
                        for y in range(ylen):
                            if mask_bff[x, y] == 1:
                                mtx_ref[x, y] = nucs_lbls[zz, x, y]
                                
                    idx_ref    =  median_filter(mtx_ref) 
                    mask_bff  *=  0
                    mtx_ref   *=  0
                    for x in range(xlen):
                        for y in range(ylen):
                            if nucs_lbls[zz, x, y] == idx_ref:
                                mask_bff[x, y]            = 1
                                nucs_lbls_piled[zz, x, y] = new_tag
                                nucs_lbls[zz, x, y]       = 0
                                
                else:
                    break
            
            #nucs_lbls  *=  1 - np.sign(nucs_lbls_piled)                          # remove the well reconstructed nucleus
            new_tag    +=  1                                                        # update the tag
            
    return  nucs_lbls, nucs_lbls_piled


# import numpy as np
# cimport numpy as np

# DTYPE = np.int
# ctypedef np.int_t DTYPE_t


# cpdef spts_int_vol(np.ndarray[DTYPE_t, ndim=3] a, np.ndarray[DTYPE_t, ndim=3] raw, list i_in, list i_out):

#     cdef np.ndarray[DTYPE_t, ndim=2] d = np.zeros([a.shape[1], a.shape[2]], dtype=np.int)
#     cdef np.ndarray[DTYPE_t, ndim=2] s = np.zeros([a.shape[1], a.shape[2]], dtype=np.int)

#     cdef int j_in
#     cdef int j_out
#     cdef int zlen  =  a.shape[0]
#     cdef int xlen  =  a.shape[1]
#     cdef int ylen  =  a.shape[2]
#     cdef int x, y, z

#     for j_in in i_in:
#         for z in range(zlen):
#             for x in range(xlen):
#                 for y in range(ylen):
#                     if a[z, x, y] == j_in:
#                         s[x, y] += 1
#                         d[x, y] += raw[z, x, y]
         
#     return s, d
