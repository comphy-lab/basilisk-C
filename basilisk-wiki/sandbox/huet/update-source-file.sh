#!/bin/bash

file=$1
source_path="lagrangian_caps/"

sed -i "s|lag-mesh.h|capsule-ft.h|g" $file
sed -i "s|neo-hookean.h|neo-hookean-ft.h|g" $file
sed -i "s|skalak.h|skalak-ft.h|g" $file
sed -i "s|caps-viscosity.h|viscosity-ft.h|g" $file
sed -i "s|caps_viscosity.h|viscosity-ft.h|g" $file
sed -i "s|common-shapes.h|common-shapes-ft.h|g" $file

sed -i "s|\.nlp|\.nln|g" $file
sed -i "s|->nlp|->nln|g" $file
sed -i "s|MB(|CAPS(|g" $file
sed -i "s|mbs.nbmb|NCAPS|g" $file
sed -i "s|mbs.mb|allCaps.caps|g" $file
sed -i "s|mbs.mb[i]|CAPS(i)|g" $file
sed -i "s|mbs.mb[j]|CAPS(j)|g" $file
sed -i "s|mbs.mb[k]|CAPS(k)|g" $file

sed -i "s|initialize_membranes(|initialize_capsules(|g" $file
sed -i "s|initialize_spherical_mb(|activate_spherical_capsule(|g" $file
sed -i "s|initialize_biconcave_mb(|activate_biconcave_capsule(|g" $file
sed -i "s|initialize_rbc_mb(|activate_biconcave_capsule(|g" $file
sed -i "s|initialize_membrane_stencils(|initialize_capsule_stencils(|g" $file
sed -i "s|initialize_membranes_stencils(|initialize_all_capsules_stencils(|g" $file

sed -i "s|dump_membranes(|dump_capsules(|g" $file
sed -i "s|restore_membranes(|restore_capsules(|g" $file