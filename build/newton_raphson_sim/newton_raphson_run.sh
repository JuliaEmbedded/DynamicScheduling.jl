#! /bin/sh

rm -r work
vlib work

vlog -work work +acc=blnr -sv -noincr -timescale 1ns/1ps nr_sim_tb.sv

vcom -2008 /home/bb2515/dynamatic/components/elastic_components.vhd
vcom -2008 /home/bb2515/dynamatic/components/delay_buffer.vhd
vcom -2008 /home/bb2515/dynamatic/components/multipliers.vhd
vcom -2008 /home/bb2515/dynamatic/components/mul_wrapper.vhd

vcom -2008 intop_sdiv_32ns_32ns_32_36_1.vhd
vcom -2008 arithmetic_units_dss.vhd

vcom -2008 /home/bb2515/dynamatic/components/arithmetic_units.vhd
vcom -2008 /home/bb2515/dynamatic/components/MemCont.vhd
vcom -2008 newton_raphson_graph.vhd

vopt -work work nr_sim_tb_top -o work_opt
vsim -modelsimini comp_sim_libs/modelsim.ini nr_sim_tb_top -do doit.tcl
