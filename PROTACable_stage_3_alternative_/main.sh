#/bin/bash



python -W ignore generate_with_protein.py --linker_size 5,20 --model models/pockets_difflinker_full.ckpt --anchors 0,0 --protein test/pv_out_nhyd99_mae_complex_6vs3_11_1ERR_14_pp_model_xx14-out-out-2.pdb --fragments test/ligands/frags_new.pdb

