
# def SIMULATION(* save,  * dist,  * capitalout,  K[size_asset]):


# {
#     * distold, critdist, distverif, weight

#     distold = (*)calloc((ARRLL_dim), sizeof())
#     null(distold, ARRLL_dim)

#     int isave, asset_index, risk_index, risk_indexnext, risk_pre_index, risk_pre_indexnext, laborincome_index, laborincome_indexnext, laborincome_pre_index, laborincome_pre_indexnext, iter, state_current_pre, state_next_pre, state_current, state_next

#     critdist = 1.0
#     iter = 0
#     # save[(0, 0, 0, 0)] = 0.001
#     while (critdist > epsdist & & iter < 300)
#     {
#         copy(dist, distold, ARRLL_dim)
#         null(dist, ARRLL_dim)

#         # std::cout << critdist << "\n"

#         # distribution dynamics
#         for (asset_index=0 asset_index < size_asset asset_index++)
#         {
#             for (risk_index=0 risk_index < size_risk risk_index++)
#             {
#                 for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#                 {
#                     for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#                     {
#                         for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                         {
#                             if (distold[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)] > 0)
#                             {
#                                 state_current = risk_index * size_risk + laborincome_index
#                                 state_current_pre = risk_pre_index * size_risk + laborincome_pre_index

#                                 isave = min((int)(floor(getgrid(save[(
#                                     asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]))), size_asset - 2)
#                                 weight = (save[(asset_index, risk_index, risk_pre_index,
#                                           laborincome_index, laborincome_pre_index)] - K[isave]) / (K[isave + 1] - K[isave])
#                                 # std::cout << "weight=" << weight << "\n"

#                                 for (risk_indexnext=0 risk_indexnext < size_risk risk_indexnext++)
#                                 {
#                                     for (risk_pre_indexnext=0 risk_pre_indexnext < size_risk risk_pre_indexnext++)
#                                     {
#                                         for (laborincome_indexnext=0 laborincome_indexnext < size_laborincome laborincome_indexnext++)
#                                         {
#                                             for (laborincome_pre_indexnext=0 laborincome_pre_indexnext < size_laborincome laborincome_pre_indexnext++)
#                                             {
#                                                 state_next_pre = risk_pre_indexnext * size_risk + laborincome_pre_indexnext
#                                                 state_next = risk_indexnext * size_risk + laborincome_indexnext
#                                                 dist[(isave, risk_indexnext, risk_pre_indexnext, laborincome_indexnext, laborincome_pre_indexnext)] += (1.0 - weight) * risk_labor_trans[state_current][state_next] *
#                                                 risk_labor_pre_trans[state_current_pre][state_next_pre] * distold[(
#                                                     asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                                                 dist[(min(isave + 1, size_asset - 1), risk_indexnext, risk_pre_indexnext, laborincome_indexnext, laborincome_pre_indexnext)] += (
#                                                     weight)*risk_labor_trans[state_current][state_next] * risk_labor_pre_trans[state_current_pre][state_next_pre] * (distold[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)])
#                                             }
#                                         }
#                                     }
#                                 }
#                             }
#                         }
#                     }
#                 }
#             }
#         }

#         # convergence
#         critdist = 0.0
#         distverif = 0.0

#         for (asset_index=0 asset_index < size_asset asset_index++)
#         {
#             for (risk_index=0 risk_index < size_risk risk_index++)
#             {
#                 for (risk_pre_index=0 risk_pre_index < size_risk risk_pre_index++)
#                 {
#                     for (laborincome_index=0 laborincome_index < size_laborincome laborincome_index++)
#                     {
#                         for (laborincome_pre_index=0 laborincome_pre_index < size_laborincome laborincome_pre_index++)
#                         {
#                             critdist = (max(critdist, abs(dist[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)
#                                                                ] - distold[(asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)])))
#                             distverif += dist[(
#                                 asset_index, risk_index, risk_pre_index, laborincome_index, laborincome_pre_index)]
#                         }
#                     }
#                 }
#             }
#         }

#         iter++
#     }
#     std: : cout << "iteration=" << iter << ", critdist=" << critdist << ", distverify=" << distverif << "\n"

#     # *capitalout = 0.0

#     # for (asset_index = 0 asset_index < size_asset asset_index++)
#     # {
#     #     for (risk_index = 0 risk_index < size_risk risk_index++)
#     #     {
#     #         for (risk_pre_index = 0 risk_pre_index < size_risk risk_pre_index++)
#     #         {
#     #             for (laborincome_index = 0 laborincome_index < size_laborincome laborincome_index++)
#     #             {
#     #                 *capitalout += dist[(asset_index, risk_index, risk_pre_index, laborincome_index)] * K[asset_index]
#     #             }
#     #         }
#     #     }
#     # }
# }
