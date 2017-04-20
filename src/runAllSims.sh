#!/bin/bash
mkdir ../res/cavity_re100 ../res/cylinder_re1 ../res/foil_neutral ../res/poiseuille ../res/cavity_re1000 ../res/cylinder_re100 ../res/foil_pitched ../res/couette ../res/cylinder_re20 ../res/foil_stall
#./lbm ../input/couette_31.in
#./lbm ../input/poiseuille_31.in
./lbm ../input/foil_neutral.in
./lbm ../input/foil_pitched.in
./lbm ../input/foil_stall.in
./lbm ../input/cylinder_401_re1.in
./lbm ../input/cylinder_401_re20.in
./lbm ../input/cylinder_401_re100.in
#./lbm ../input/cavity_501_re100.in
#./lbm ../input/cavity_501_re1000.in
