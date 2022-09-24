def POLICY(VF_final,  dVF_final,  save_final,  VF,  dVF,  save,  Portfolio,  K,  Omega,  wagerate):

    # INITIALIZATION #

    VFendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                      size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice))
    cohendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                       size_laborincome, size_shock, size_shock, size_shock, size_shock, size_portfoliochoice))

    VF_final_old = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                            size_laborincome, size_shock, size_shock, size_shock, size_shock))

    iter = 0

    critV = 10000.0
    # std::cout << "iter\t"
    #           << "critV\n"

    while critV > epsV and iter < 250:
        # we need copy to make a separate object

        VFendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                          size_laborincome, size_shock, size_shock, size_portfoliochoice))
        cohendo = np.zeros((size_asset, size_risk, size_risk, size_laborincome,
                           size_laborincome, size_shock, size_shock, size_portfoliochoice))
        VF_final_old = VF_final

        # std::cout << std::setprecision(16) << VF[(5, 5)] << "\n"
        # std::cout << std::setprecision(16) << VFnew[(5, 5)] << "\n"

        # main EGM computation
        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):

                                                tempnext = 0
                                                dtempnext = 0

                                                for risk_indexnext in list(range(size_risk)):
                                                    for laborincome_indexnext in list(range(size_laborincome)):
                                                        tempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * \
                                                            VF[asset_index, risk_indexnext, risk_index,
                                                                laborincome_indexnext, laborincome_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index]
                                                        dtempnext += risk_trans[risk_index][risk_indexnext] * laborincome_trans[laborincome_index][laborincome_indexnext] * \
                                                            dVF[asset_index, risk_indexnext, risk_index,
                                                                laborincome_indexnext, laborincome_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index]

                                                cohendo[asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] = K[asset_index] + inv_MU(betapar * dtempnext)
                                                VFendo[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] = U(
                                                    cohendo[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index] - K[asset_index]) + betapar * tempnext

        # std::cout << "EGM done\n"

        # rescalings

        for portfoliochoice_index in list(range(size_portfoliochoice+1)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            shockstate_current = riskshock_index * size_shock + laborshock_index

                                            threshold_ii = 0

                                            for asset_index in list(range(size_asset)):

                                                # method 1: cash on hand
                                                cohexo = (1.0 + (r_f + pi + risk_states[risk_index])*riskshock_states[shockstate_current] * Omega[portfoliochoice_index] + r_f * (
                                                    1 - Omega[portfoliochoice_index])) * K[asset_index] + wagerate * laborincome_states[laborincome_index]*laborshock_states[shockstate_current]

                                                if cohexo < cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                     laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]:
                                                    save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = K[0]
                                                    VF[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = U(cohexo - save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) + (
                                                        VFendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - U((cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                           laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - K[0])))

                                                if cohexo >= cohendo[(0, risk_index, risk_pre_index, laborincome_index,
                                                                     laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]:

                                                    itest = threshold_ii

                                                    while (itest < size_asset) and cohexo > cohendo[((itest, risk_index, risk_pre_index, laborincome_index,
                                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index))]:
                                                        itest += 1

                                                    if (itest == size_asset):
                                                        # extrapolation
                                                        vfweight = (cohexo - cohendo[(size_asset - 2, risk_index, risk_pre_index, laborincome_index,
                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[(
                                                                                          size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[(size_asset - 2, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                                                  laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                        igridL = size_asset - 2
                                                        igridH = size_asset - 1

                                                    else:
                                                        # standard interior
                                                        vfweight = (cohexo - cohendo[(itest - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                      laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (cohendo[(
                                                                                          itest, risk_index, risk_pre_index, laborincome_index,
                                                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - cohendo[(itest - 1, risk_index, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                                                  laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                        igridL = itest - 1
                                                        igridH = itest - 0

                                                    VF[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                        laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, VFendo[(
                                                            igridL, risk_index, risk_pre_index, laborincome_index,
                                                            laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VFendo[(igridH, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                                                  laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)])
                                                    save[(asset_index, risk_index, risk_pre_index, laborincome_index,
                                                          laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = inter1d(vfweight, K[igridL], K[igridH])

                                                    threshold_ii = min(
                                                        size_asset - 2, itest)

        # std::cout << "rescaling done\n"

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            shockstate_pre = riskshock_index * size_shock + laborshock_index

                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):
                                                tempnext = 0.0

                                                for risk_pre_indexnext in list(range(size_risk)):
                                                    for laborincome_pre_indexnext in list(range(size_laborincome)):
                                                        for riskshock_pre_indexnext in list(range(size_shock)):
                                                            for laborshock_pre_indexnext in list(range(size_shock)):

                                                                shockstate_current = riskshock_pre_indexnext * \
                                                                    size_shock + laborshock_pre_indexnext
                                                                tempnext += risk_trans[risk_pre_index][risk_pre_indexnext]*laborincome_trans[laborincome_pre_index][laborincome_pre_indexnext] * risklaborshock_trans[shockstate_current] * VF[(
                                                                    asset_index, risk_pre_indexnext, risk_pre_index, laborincome_pre_indexnext, laborincome_pre_index, riskshock_pre_indexnext, riskshock_pre_index, laborshock_pre_indexnext, laborshock_pre_index, portfoliochoice_index)]

                                                # std::cout << tempnext << "\n"
                                                if (portfoliochoice_index == 0):
                                                    temp = tempnext
                                                    itemp = 0

                                                if (tempnext > temp):

                                                    temp = tempnext
                                                    itemp = portfoliochoice_index

                                            VF_final[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = VF[
                                                asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp]
                                            save_final[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = save[
                                                asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, itemp]
                                            Portfolio[asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                      riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index] = Omega[itemp]

        # std::cout << "port done\n"
        # std::cout << std::setprecision(16) << VF[(5, 5)] << "\n"

        # computing new derivatives and convergence
        critV = 0.0

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):

                                            if (asset_index >= 2):
                                                dVF_final[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = nderiv(VF_final[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                             riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[(
                                                                                                                                                                 asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                 riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                            critV = max(critV, abs(VF_final[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                             riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)
                                                                            ] - VF_final_old[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]))

                                            # left corner
                                            dVF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                       riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(
                                                           1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                           riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                      riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[1] - K[0])
                                            # right corner
                                            dVF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                       riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] = (VF_final[(size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                   laborincome_pre_index,
                                                                                                                                                   riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] - VF_final[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)]) / (K[size_asset - 1] - K[size_asset - 2])

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = relaxVF * VF_final[(asset_index, risk_index, risk_pre_index,
                                                                                                                                                                            laborincome_index, laborincome_pre_index,
                                                                                                                                                                            riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index)] + (1 - relaxVF) * VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                                 riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]

        for asset_index in list(range(size_asset)):
            for risk_index in list(range(size_risk)):
                for risk_pre_index in list(range(size_risk)):
                    for laborincome_index in list(range(size_laborincome)):
                        for laborincome_pre_index in list(range(size_laborincome)):
                            for riskshock_index in list(range(size_shock)):
                                for riskshock_pre_index in list(range(size_shock)):
                                    for laborshock_index in list(range(size_shock)):
                                        for laborshock_pre_index in list(range(size_shock)):
                                            for portfoliochoice_index in list(range(size_portfoliochoice+1)):

                                                if (asset_index >= 2):
                                                    dVF[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                         riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = nderiv(VF[(asset_index - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(asset_index - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], VF[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index, riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)], K[asset_index - 2], K[asset_index - 1], K[asset_index])

                                                # left corner
                                                dVF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                     riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                  laborincome_pre_index,
                                                                                                                                                                  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(0, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[1] - K[0])
                                                # right corner
                                                dVF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                     riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] = (VF[(size_asset - 1, risk_index, risk_pre_index, laborincome_index,
                                                                                                                                                                  laborincome_pre_index,
                                                                                                                                                                  riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)] - VF[(size_asset - 2, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index,
                                                                                                                                                                                                                                                                              riskshock_index, riskshock_pre_index, laborshock_index, laborshock_pre_index, portfoliochoice_index)]) / (K[size_asset - 1] - K[size_asset - 2])

        iter += 1
        print("iteration={:d}, critV={:d}".format(iter, critV))
    return VF_final, dVF_final, save_final, VF, dVF, save, Portfolio
