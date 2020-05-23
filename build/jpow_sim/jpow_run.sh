#! /bin/sh

rm -r work
vlib work

vlog -work work +acc=blnr -sv -noincr -timescale 1ns/1ps jpow_sim_tb.sv

vcom -2008 /home/bb2515/dynamatic/components/elastic_components.vhd
vcom -2008 /home/bb2515/dynamatic/components/delay_buffer.vhd
vcom -2008 /home/bb2515/dynamatic/components/multipliers.vhd
vcom -2008 /home/bb2515/dynamatic/components/mul_wrapper.vhd
vcom -2008 /home/bb2515/dynamatic/components/arithmetic_units.vhd
vcom -2008 /home/bb2515/dynamatic/components/MemCont.vhd
vcom -2008 jpow_graph.vhd

vopt -work work jpow_sim_tb_top -o work_opt
vsim jpow_sim_tb_top -do doit.tcl