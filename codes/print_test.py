key_list = ['sc_r_iqr_50', 'bm_iqr_10', 'bm_iqr_25', 'bm_iqr_50', 'bm_iqr_75',
       'bm_iqr_90', 'np_iqr_10', 'np_iqr_25', 'np_iqr_50', 'np_iqr_75',
       'np_iqr_90', 'vp_m_iqr_10', 'vp_m_iqr_25', 'vp_m_iqr_50', 'vp_m_iqr_75',
       'vp_m_iqr_90', 'Tp_iqr_10', 'Tp_iqr_25', 'Tp_iqr_50', 'Tp_iqr_75',
       'Tp_iqr_90', 'proton_beta_iqr_10', 'proton_beta_iqr_25',
       'proton_beta_iqr_50', 'proton_beta_iqr_75', 'proton_beta_iqr_90',
       'vA_iqr_10', 'vA_iqr_25', 'vA_iqr_50', 'vA_iqr_75', 'vA_iqr_90',
       'alfven_ratio_iqr_10', 'alfven_ratio_iqr_25', 'alfven_ratio_iqr_50',
       'alfven_ratio_iqr_75', 'alfven_ratio_iqr_90']

# ind1 = 2
ind1 = input('Enter the index: ')
ind1 = int(ind1)
ind2 = ind1 + 1
for key in key_list[ind1:ind2]:
    print(key)
    for xx in dfn[key]:
        #print(f"{np.round(xx, 5)} \n \n")
        print(('{:.2e}\n\n'.format(xx)))
