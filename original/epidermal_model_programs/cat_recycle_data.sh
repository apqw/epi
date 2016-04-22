#!/bin/bash

cat last_data_uvp last_data_a last_data_B last_data_w > cat_recycle_data_tmp.dat
awk ' NF > 0 { print $0 } ' cat_recycle_data_tmp.dat > recycled_data
rm -f cat_recycle_data_tmp.dat
