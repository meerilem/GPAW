#!/bin/bash

cation=( Pyr14 BMIm EMIm TEPA BPy )
anion=( PF6 TFSI Br Cl I BF4 BCN4 FSI )
#cation=( EMIm )
#anion=( Cl )

home=$(pwd)
for i in ${cation[@]}; do
    for j in ${anion[@]}; do
        mkdir $i$j/
	mv $i$j\_1s.out $home/$i$j/
    done
done

