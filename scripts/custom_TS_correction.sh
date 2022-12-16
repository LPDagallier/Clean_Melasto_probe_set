#!/bin/bash

###################################################
#### 7 Individual stop codons removal
###################################################
cat TS_to_correct.FNA > TS_to_correct_CORRECTED.FNA

## WWQZ_GW_QUTB-4992
sed -i '/>WWQZ_GW_QUTB-4992/,+1s/TAAAAATCGTCAAGGATG//g' TS_to_correct_CORRECTED.FNA

## NIJU_GW_CLMX-5271
sed -i '/>NIJU_GW_CLMX-5271/,+1s/TAGAGCGGCAACCCGCATGCTGTTGGATTT//g' TS_to_correct_CORRECTED.FNA

## SWGX-5326
sed -i '/>SWGX-5326/,+1s/TGAAGG$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-5326
sed -i '/>WWQZ-5326/,+1s/TGAAGA$//g' TS_to_correct_CORRECTED.FNA

## UAXP-5339
sed -i '/>UAXP-5339/,+1s/TGAATGGAGGAA$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-5355
sed -i '/>WWQZ-5355/,+1s/TGACCGGGAAACCCTTTGTTTCTCTA//g' TS_to_correct_CORRECTED.FNA

## SWGX-5578
sed -i '/>SWGX-5578/,+1s/TAG$//g' TS_to_correct_CORRECTED.FNA

## BZDF-5421
sed -i '/>BZDF-5421/,+1s/TAGTCCCCT$//g' TS_to_correct_CORRECTED.FNA

## NIJU-5594
sed -i '/>NIJU-5594/,+1s/TGA$//g' TS_to_correct_CORRECTED.FNA

## UAXP-5594
sed -i '/>UAXP-5594/,+1s/TAA$//g' TS_to_correct_CORRECTED.FNA

## YNUE-5596
sed -i '/>YNUE-5596/,+1s/A$//g' TS_to_correct_CORRECTED.FNA

## SWGX_GW_BERS-5634
sed -i '/>SWGX_GW_BERS-5634/,+1s/CT$//g' TS_to_correct_CORRECTED.FNA

## UAXP_EG_Ambtr-5634
sed -i '/>UAXP_EG_Ambtr-5634/,+1s/CT$//g' TS_to_correct_CORRECTED.FNA

## WWQZ_EG_Ambtr-5634
sed -i '/>WWQZ_EG_Ambtr-5634/,+1s/CT$//g' TS_to_correct_CORRECTED.FNA

## NIJU_GW_Orysa-5857
sed -i '/>NIJU_GW_Orysa-5857/,+1s/TAGAGCAAGAAC$//g' TS_to_correct_CORRECTED.FNA

## RJNQ-5910
sed -i '/>RJNQ-5910/,+1s/G$//g' TS_to_correct_CORRECTED.FNA

## FGDU-5941
sed -i '/>FGDU-5941/,+1s/TA$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-5941
sed -i '/>WWQZ-5941/,+1s/TAAAGCAATTGGCACCCACACCA$//g' TS_to_correct_CORRECTED.FNA

## BZDF-5981
sed -i '/>BZDF-5981/,+1s/TAAAGAAAAAAAAAATGTCGGTTGTATCAT//g' TS_to_correct_CORRECTED.FNA

## RJNQ_GW_Orysa-5981
sed -i '/>RJNQ_GW_Orysa-5981/,+1s/TGAATG$//g' TS_to_correct_CORRECTED.FNA

## YNUE-5981
sed -i '/>YNUE-5981/,+1s/TGAATG$//g' TS_to_correct_CORRECTED.FNA

## FEDW_GW_Orysa-6048
sed -i '/>FEDW_GW_Orysa-6048/,+1s/TGATATGTTCAA$//g' TS_to_correct_CORRECTED.FNA

## AYMT_GW_MTHW-6110 
sed -i '/>AYMT_GW_MTHW-6110/,+1s/T$//g' TS_to_correct_CORRECTED.FNA

## SWGX_GW_MTHW-6110
sed -i '/>SWGX_GW_MTHW-6110/,+1s/T$//g' TS_to_correct_CORRECTED.FNA

## BZDF-6366
sed -i '/>BZDF-6366/,+1s/A$//g' TS_to_correct_CORRECTED.FNA

## YNUE-6393
sed -i '/>YNUE-6393/,+1s/TAG$//g' TS_to_correct_CORRECTED.FNA

## BZDF-6401
sed -i '/>BZDF-6401/,+1s/T$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-6450
sed -i '/>WWQZ-6450/,+1s/TGAGGC$//g' TS_to_correct_CORRECTED.FNA

## BZDF_GW_Orysa-6506
sed -i '/>BZDF_GW_Orysa-6506/,+1s/A$//g' TS_to_correct_CORRECTED.FNA

## FGDU_GW_Orysa-6506
sed -i '/>FGDU_GW_Orysa-6506/,+1s/A$//g' TS_to_correct_CORRECTED.FNA

## RJNQ_GW_Orysa-6506
sed -i '/>RJNQ_GW_Orysa-6506/,+1s/A$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-6506
sed -i '/>WWQZ-6506/,+1s/T$//g' TS_to_correct_CORRECTED.FNA

## WWQZ_GW_Ambtr-6527
sed -i '/>WWQZ_GW_Ambtr-6527/,+1s/TT$//g' TS_to_correct_CORRECTED.FNA

## NIJU_GW_UMUL-6528
sed -i '/>NIJU_GW_UMUL-6528/,+1s/TAG$//g' TS_to_correct_CORRECTED.FNA

## AYMT-6538
sed -i '/>AYMT-6538/,+1s/TGAGAGCAATGCTGGATTCTTTTCAAGTT$//g' TS_to_correct_CORRECTED.FNA

## SWGX-6620
sed -i '/>SWGX-6620/,+1s/TAG$//g' TS_to_correct_CORRECTED.FNA

## FGDU_GW_Ambtr-6667
sed -i '/>FGDU_GW_Ambtr-6667/,+1s/AG$//g' TS_to_correct_CORRECTED.FNA

## RJNQ_GW_Ambtr-6667
sed -i '/>RJNQ_GW_Ambtr-6667/,+1s/TAATGTTGGCTATTGCCAGAGGATCAAGCA$//g' TS_to_correct_CORRECTED.FNA

## AYMT-6780
sed -i '/>AYMT-6780/,+1s/TAGCGAGCGGGGAGCAAGTCAAT$//g' TS_to_correct_CORRECTED.FNA

## SWGX-6780
sed -i '/>SWGX-6780/,+1s/TAA$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-6780
sed -i '/>WWQZ-6780/,+1s/TGAGGAACAAGTTTAG$//g' TS_to_correct_CORRECTED.FNA

## WWQZ_GW_Ambtr-6864
sed -i '/>WWQZ_GW_Ambtr-6864/,+1s/G$//g' TS_to_correct_CORRECTED.FNA

## UAXP-6933
sed -i '/>UAXP-6933/,+1s/G$//g' TS_to_correct_CORRECTED.FNA

## NIJU-6962
sed -i '/>NIJU-6962/,+1s/TGAAGG$//g' TS_to_correct_CORRECTED.FNA

## AJJE-6969
sed -i '/>AJJE-6969/,+1s/TGAAATGAAGTAGATCAATCCA$//g' TS_to_correct_CORRECTED.FNA

## AYMT-6969
sed -i '/>AYMT-6969/,+1s/TGAGTTCGTGACAGAGACGTTC$//g' TS_to_correct_CORRECTED.FNA

## FGDU-6969
sed -i '/>FGDU-6969/,+1s/TGAGTTCGCAATAAAGCCATTT$//g' TS_to_correct_CORRECTED.FNA

## UAXP-6969
sed -i '/>UAXP-6969/,+1s/AC$//g' TS_to_correct_CORRECTED.FNA

## FEDW-7111
sed -i '/>FEDW-7111/,+1s/TAA$//g' TS_to_correct_CORRECTED.FNA

## WWQZ-7111
sed -i '/>WWQZ-7111/,+1s/TGA$//g' TS_to_correct_CORRECTED.FNA

## FEDW_GW_Orysa-7241
sed -i '/>FEDW_GW_Orysa-7241/,+1s/TAGAGAAAC$//g' TS_to_correct_CORRECTED.FNA

## YNUE-7336
sed -i '/>YNUE-7336/,+1s/TAG$//g' TS_to_correct_CORRECTED.FNA

