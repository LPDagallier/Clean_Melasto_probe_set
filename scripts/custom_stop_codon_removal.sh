#!/bin/bash

###################################################
#### 7 Individual stop codons removal
###################################################

## Brachyotum__5699__hc-5699
sed -i '/>Brachyotum__5699__hc-5699/,+1s/TTGCCCCTCAAATGGTAGTATCCC/TTGCCCCTCAAATGGNANTATCCC/g' final_sequences_set.FNA

## Tibouchina__AT3G59040__hc-LOC00667
sed -i '/>Tibouchina__AT3G59040__hc-LOC00667/,+1s/CGATCTGAAAKCTGGAGGAACATGAGATAA/CGATCTGAAAKCTGGAGGAACATGAGANAA/g' final_sequences_set.FNA

## Affzelli__5449__hit-LOC00776
sed -i '/>Affzelli__5449__hit-LOC00776/,+1s/ACTTGGGCTAAATTGTTGTGTTTTTGA/ACTTGGGCTAAATTGTTGTGTTTTNNN/g' final_sequences_set.FNA

## Memecylon_tr_combo__AT1G02020__hc-LOC01351
sed -i '/>Memecylon_tr_combo__AT1G02020__hc-LOC01351/,+1s/TGCAACAAAAGTAAATCCTGGTTGTATTTTATGTAG//g' final_sequences_set.FNA

## Memecylon_tr_combo__AT4G39520__hc-LOC01520
sed -i '/>Memecylon_tr_combo__AT4G39520__hc-LOC01520/,+1s/AGGCTAAACAAGTAGCCACCT/AGGCTAAACAAGNAGCCACCT/g' final_sequences_set.FNA

## Tibouchina__TC05G024510__hc-LOC04226
sed -i '/>Tibouchina__TC05G024510__hc-LOC04226/,+1s/GGKTAGTTTATTCATTGAATTKKMSYT//g' final_sequences_set.FNA

## Tibouchina__AT4G00740__hc-LOC04935
sed -i '/>Tibouchina__AT4G00740__hc-LOC04935/,+1s/GCATGGTAACTTCTGATGAAACATTTTGTA/GCATGGTANCTTCTGATGAAACATTTTGTA/g' final_sequences_set.FNA

## Tibouchina__5639__hc-LOC05244
sed -i '/>Tibouchina__5639__hc-LOC05244/,+1s/TTKSTCTAGGGCAAT/TTKSTCTNGGGCAAT/g' final_sequences_set.FNA

## Affzelli__5426__hit-LOC06188
sed -i '/>Affzelli__5426__hit-LOC06188/,+1s/ATGTGAAAAAATGTTGGTAGT/ATGTGGAAAAATGTTGGTAGT/g' final_sequences_set.FNA

## Tibouchina__NM_179236.3__hc-LOC06614
sed -i '/>Tibouchina__NM_179236.3__hc-LOC06614/,+1s/TCAAGGTAACTTTTTGTGCTTTTTGCA/TCAAGGTNNCTTTTTGTGCTTTTTGCA/g' final_sequences_set.FNA
sed -i '/>Tibouchina__NM_179236.3__hc-LOC06614/,+1s/YTSCTGCTGCTATGGTGACAG//g' final_sequences_set.FNA

## Tibouchina__5339__hc-LOC07010
sed -i '/>Tibouchina__5339__hc-LOC07010/,+1s/CAKACCAGGKATKKWTYTCMMCTTRAAAAAATCTAAACGAGGCAGTTG//g' final_sequences_set.FNA

## Memecylon_tr_combo__TC04G019470__hc-LOC07251
sed -i '/>Memecylon_tr_combo__TC04G019470__hc-LOC07251/,+1s/TCTATGTAGAGATGGGTTGCAGTTGTTTTTGTATTC/TCTATGTNNAGATGGGTTGCAGTTGTTTTTGTATTC/g' final_sequences_set.FNA

## Tibouchina__6563__hc-LOC07413
sed -i '/>Tibouchina__6563__hc-LOC07413/,+1s/TCCAAGCTTCCAAAGTAAAGGACTTCTGTGGTGACAATATTG/TCCAAGCTTCCAAAGNANAGGACTTCTGTGGTGACAATATTG/g' final_sequences_set.FNA

## Torricelli_tr_combo__TC01G021540__hc-LOC07589
sed -i '/>Torricelli_tr_combo__TC01G021540__hc-LOC07589/,+1s/TCATGACTTTTCTGTGCATCCATCAACTTCCTA/TCANGNCTTTTCTGTGCATCCATCAACTTCCTA/g' final_sequences_set.FNA

## Tibouchina__NM_001198514.1__hc-LOC08266
sed -i '/>Tibouchina__NM_001198514.1__hc-LOC08266/,+1s/TTYWYSSYSTYYYSTKYKTTCTGA//g' final_sequences_set.FNA

## Memecylon_tr_combo__TC09G010970__hc-LOC08580
sed -i '/>Memecylon_tr_combo__TC09G010970__hc-LOC08580/,+1s/GTAAGTGTTCATTATTTCACACGCCTTCTGAGTTGATTTTGTTGATTCTGT//g' final_sequences_set.FNA

## Tibouchina__TC05G027690__hc-LOC09163
sed -i '/>Tibouchina__TC05G027690__hc-LOC09163/,+1s/ATGATCGGTAGAGCACACGTCTGAACTCCAGTCACGAGATTC//g' final_sequences_set.FNA

## Tibouchina__7572__hc-LOC09834
sed -i '/>Tibouchina__7572__hc-LOC09834/,+1s/CAAGGAGTTATAACGTGACTCAAG//g' final_sequences_set.FNA

## Tibouchina__5942__hc-LOC14673
sed -i '/>Tibouchina__5942__hc-LOC14673/,+1s/RGCMRRRRKCMSYRSSSSYRGMARRARRCMSKWMGMCCAGAAGAAGACCGGTAG//g' final_sequences_set.FNA

## Tibouchina__AT5G13630__hc-LOC15950
sed -i '/>Tibouchina__AT5G13630__hc-LOC15950/,+1s/GTATAGCCGGCTGMSMKRGYSRTSAARRTKGTK//g' final_sequences_set.FNA

## Brachyotum__AT3G02760__hc-LOC24270
sed -i '/>Brachyotum__AT3G02760__hc-LOC24270/,+1s/AATTCCTAGAATTCATKSKKKCKYTAWKRWRRCSWKMWA//g' final_sequences_set.FNA



