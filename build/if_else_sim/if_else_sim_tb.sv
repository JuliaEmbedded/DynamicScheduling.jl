module if_else_sim_tb_top();

logic clk, rst;

//dut ctrl inputs
logic start_in, start_valid;
logic end_ready;

//dut ctrl outputs
logic start_ready;
logic end_valid;
logic [31:0] end_out;

//dut arg inputs
logic [31:0] a_din;
logic [31:0] b_din;
logic a_valid_in, b_valid_in; //not used in dut

//dut arg outputs
logic a_ready_out, b_ready_out; //not used in dut

int a,b, res;
int test_total;
int dut_res;

function int gen_result;
    input int a, b;
    automatic int k = a - b;
    if (a > b) begin
        k = k * a;
        return k;
    end else if (b > a) begin
        k = k * b;
        return k;
    end else
        return k;
endfunction

//tb program
initial begin
    clk = 1;
    rst = 1;
    #15 end_ready = 1;
    
    test_total = 100;
    for (int test_it=0; test_it < test_total; test_it++) begin
        a=$random();
        b=$random();
        res = gen_result(a, b);
        rst = 0;
        
        //drive dut
        @(posedge clk);
        $display("Test %d: a: %d, b: %d", test_it, a, b);
        a_din = a;
        b_din = b;
        start_in = 1;
        start_valid = 1;
        
        @(posedge clk); //wait two cycles for dut to get info
        @(posedge clk);
        a_din = 0;
        b_din = 0;
        start_in = 0;
        start_valid = 0;
        
        @(posedge end_valid);
        dut_res = end_out;
        $display("Test %d: Predicted: %x, Actual: %x, Decimal: %d", test_it, res, dut_res, res);
        if (res != dut_res) $error("TEST FAIL");

        @(posedge clk);
        rst = 1;
        @(posedge clk);
        
    end

    $display("TEST SUCCESS cool off.");
    @(posedge clk);
    @(posedge clk);
    
    $stop;
end


//tb clock
always #5 clk = ~clk;


if_else_graph ie1(
    .clk        (clk        ), 
    .rst        (rst        ), 
    .start_in   (start_in   ),
    .start_valid(start_valid),
    .start_ready(start_ready),
    .end_out    (end_out    ),
    .end_valid  (end_valid  ),
    .end_ready  (end_ready  ),
    .a_din      (a_din      ),
    .a_valid_in (a_valid_in ),
    .a_ready_out(a_ready_out),
    .b_din      (b_din      ),
    .b_valid_in (b_valid_in ),
    .b_ready_out(b_ready_out)
);

endmodule
