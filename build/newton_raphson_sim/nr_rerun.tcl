vlog -work work +acc=blnr -sv -noincr -timescale 1ns/1ps nr_sim_tb.sv
do newton_raphson_graph_modelsim.tcl
vopt -work work nr_sim_tb_top -o work_opt
vsim nr_sim_tb_top

add wave *
run -all
