#!/bin/bash

./generate_volume_image_script_custom.py \
-v --host 172.23.162.70 --username sberry --password sheeP6ep \
--experiment \
20171130-kim2-cytoo-scaling-FISH-Y-custom \
20171130-kim2-cytoo-scaling-FISH-Y-L \
20171130-kim2-cytoo-scaling-FISH-Y-M \
20171130-kim2-cytoo-scaling-FISH-DC-L \
20171130-kim2-cytoo-scaling-FISH-DC-S \
20171130-kim2-cytoo-scaling-FISH-DC-custom \
20171130-kim2-control-scaling-FISH-sorted \
20171130-kim2-control-scaling-FISH-unsorted \
--input_path \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y_custom/STACK/ \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y-L/STACK/ \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y-M/STACK/ \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC-L/STACK/ \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC-S/STACK/ \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC_custom/STACK/ \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH-unpatterned/TIFF/CYCLE_01/glass-sorted/STACK/ \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH-unpatterned/TIFF/CYCLE_01/glass-unsorted/STACK/ \
--output_path \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y_custom/VOL \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y-L/VOL \
~/pelkmanslab-share1/Data/Users/Scott/20171201-Kim2-MicropatternFISH/TIFF/Y-M/VOL \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC-L/VOL \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC-S/VOL \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH/TIFF/CYCLE_01/DC_custom/VOL \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH-unpatterned/TIFF/CYCLE_01/glass-sorted/VOL \
~/pelkmans-samba2/20171130-kim2-cytoo-scaling-FISH-unpatterned/TIFF/CYCLE_01/glass-unsorted/VOL \
--fname_stem \
20171130-kim2-cytoo-scaling-FISH-customY \
No\ Response\ Detected \
No\ Response\ Detected \
No\ Response\ Detected \
20171130-kim2-cytoo-scaling-FISH-DC-S \
No\ Response\ Detected \
20171130-kim2-cytoo-scaling-FISH-unpatterned-unsorted \
20171130-kim2-cytoo-scaling-FISH-unpatterned-unsorted \
--watch 10 -s multi-plate --plate plate01
