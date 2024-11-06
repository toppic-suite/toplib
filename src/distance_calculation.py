import numpy as np
import numba as nb


@nb.njit
def peak_50_cosine_mass_ppm(sorted_mz1, sorted_int1, sorted_ch1, sorted_mz2, sorted_int2, sorted_ch2, err, ppm):
    cos_unembed = 0
    for i in range(len(sorted_mz1)):
        for j in range(len(sorted_mz2)):
            tol2 = max(err, ppm * sorted_mz1[i] * 1e-6)
            th1 = 0.0 + tol2
            th2 = 1.00235 + tol2
            th3 = 1.00235 - tol2
            if sorted_ch1[i] == sorted_ch2[j]:
                if abs(sorted_mz2[j] - sorted_mz1[i]) <= th1 or (
                        (abs(sorted_mz2[j] - sorted_mz1[i]) <= th2) and (abs(sorted_mz2[j] - sorted_mz1[i]) >= th3)):
                    cos_unembed += sorted_int1[i] * sorted_int2[j]
                    if cos_unembed > 1:
                        return 1.0 - (cos_unembed - sorted_int1[i] * sorted_int2[j])

    return 1.0 - cos_unembed

    
@nb.njit
def peak_50_ed_mass_ppm(sorted_mz1, sorted_int1, sorted_ch1, sorted_mz2, sorted_int2, sorted_ch2, err, ppm):
    ed_sum = 0
    idx1 = []
    idx2 = []
    for i in range(len(sorted_mz1)):
        for j in range(len(sorted_mz2)):
            tol2 = max(err, ppm * sorted_mz1[i] * 1e-6)
            th1 = 0.0 + tol2
            th2 = 1.00235 + tol2
            th3 = 1.00235 - tol2
            if sorted_ch1[i] == sorted_ch2[j]:
                if abs(sorted_mz2[j] - sorted_mz1[i]) <= th1 or (
                        (abs(sorted_mz2[j] - sorted_mz1[i]) <= th2) and (abs(sorted_mz2[j] - sorted_mz1[i]) >= th3)):
                    ed_sum = ed_sum + (sorted_int1[i] - sorted_int2[j]) ** 2
                    idx1.append(i)
                    idx2.append(j)

    for i in range(len(sorted_mz1)):
        if i not in idx1:
            ed_sum = ed_sum + (sorted_int1[i] ** 2)

    for j in range(len(sorted_mz2)):
        if j not in idx2:
            ed_sum = ed_sum + (sorted_int2[j] ** 2)

    return np.sqrt(ed_sum)


@nb.njit
def peak_50_entropy_mass_ppm(sorted_mz1, sorted_int1, sorted_ch1, sorted_mz2, sorted_int2, sorted_ch2, err, ppm):
    data_unembed1_normalized = (0.5) * sorted_int1 / np.sum(sorted_int1)
    data_unembed2_normalized = (0.5) * sorted_int2 / np.sum(sorted_int2)

    entropy_unembed = 0
    for k in range(len(sorted_mz1)):
        for j in range(len(sorted_mz2)):
            tol2 = max(err, ppm * sorted_mz1[k] * 1e-6)
            th1 = 0.0 + tol2
            th2 = 1.00235 + tol2
            th3 = 1.00235 - tol2
            if sorted_ch1[k] == sorted_ch2[j]:
                if abs(sorted_mz2[j] - sorted_mz1[k]) <= th1 or (
                        (abs(sorted_mz2[j] - sorted_mz1[k]) <= th2) and (abs(sorted_mz2[j] - sorted_mz1[k]) >= th3)):
                    x_sum = data_unembed1_normalized[k] + data_unembed2_normalized[j]
                    fx = x_sum * np.log2(x_sum)
                    fx1 = data_unembed1_normalized[k] * np.log2(data_unembed1_normalized[k])
                    fx2 = data_unembed2_normalized[j] * np.log2(data_unembed2_normalized[j])
                    entropy_unembed += fx - fx1 - fx2
                    if entropy_unembed > 1:
                        return 1.0 - (entropy_unembed - (fx - fx1 - fx2))                                        
                    
    return 1.0 - entropy_unembed


def dist_func_cal(sorted_mz_all, sorted_int_all, sorted_ch_all, err, ppm, metric):
    # calculate distance of different functions
    func_dist = np.zeros((len(sorted_int_all),len(sorted_int_all)))
    for i in range(len(sorted_int_all)):
        for j in range(i + 1, len(sorted_int_all)):  
            if metric=='euclidean':
                func_dist[i,j] = peak_50_ed_mass_ppm(sorted_mz_all[i], sorted_int_all[i], sorted_ch_all[i], sorted_mz_all[j], sorted_int_all[j], sorted_ch_all[j], err, ppm)
            elif metric =='cosine':
                func_dist[i,j] = peak_50_cosine_mass_ppm(sorted_mz_all[i], sorted_int_all[i], sorted_ch_all[i], sorted_mz_all[j], sorted_int_all[j], sorted_ch_all[j], err, ppm)
            else:
                func_dist[i,j] = peak_50_entropy_mass_ppm(sorted_mz_all[i], sorted_int_all[i], sorted_ch_all[i], sorted_mz_all[j], sorted_int_all[j], sorted_ch_all[j], err, ppm)
            func_dist[j,i] = func_dist[i,j]

    return func_dist[np.triu_indices(len(sorted_int_all),k=1)]
    
