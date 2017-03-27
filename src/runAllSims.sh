 #!/bin/bash
mkdir cavity_re100 cylinder_re10 foil_neutral poiseuille cavity_re1000 cylinder_re100 foil_pitched couette cylinder_re50 foil_stall
./lbm ../input/couette_31.in
./lbm ../input/poiseuille_31.in
./lbm ../input/foil_neutral.in
./lbm ../input/foil_pitched.in
./lbm ../input/foil_stall.in
./lbm ../input/cylinder_401_re10.in
./lbm ../input/cylinder_401_re50.in
./lbm ../input/cylinder_401_re100.in
./lbm ../input/cavity_501_re100.in
./lbm ../input/cavity_501_re1000.in