vlog -work work +acc=blnr -sv -noincr -timescale 1ns/1ps if_else_sim_tb.sv
do if_else_graph_modelsim.tcl
vopt -work work if_else_sim_tb_top -o work_opt
vsim if_else_sim_tb_top

add wave *
run -all
