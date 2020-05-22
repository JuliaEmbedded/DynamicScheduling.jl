vlog -work work +acc=blnr -sv -noincr -timescale 1ns/1ps mul_sim_tb.sv
do mul_graph_modelsim.tcl
vopt -work work mul_sim_tb_top -o work_opt
vsim mul_sim_tb_top

add wave *
run -all
