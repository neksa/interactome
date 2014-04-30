#!/bin/bash

mkdir -p plots/contacts_in_charged_aa

for loc in backbone sidechain mixed
do
    for var in ARG_ARG ARG_LYS ASP_ARG ASP_ASP ASP_GLU ASP_LYS GLU_ARG GLU_GLU GLU_LYS LYS_LYS
    do
        mkdir -p plots/contacts_in_charged_aa/$loc
        cp plots/contact_probabilities/$loc/$var.png   plots/contacts_in_charged_aa/$loc/$var.png
    done
done
