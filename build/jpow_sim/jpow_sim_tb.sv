module jpow_sim_tb_top();

logic clk, rst;

//dut ctrl inputs
logic start_in, start_valid;
logic end_ready;

//dut ctrl outputs
logic start_ready;
logic end_valid;
logic [31:0] end_out;

//dut arg inputs
logic [31:0] x_din;
logic [31:0] n_din;
logic x_valid_in, n_valid_in; //not used in dut

//dut arg outputs
logic x_ready_out, n_ready_out; //not used in dut

int x,n, res;
int test_total;
int dut_res;

function int gen_result;
    input int x, n;
    automatic int r = 1;
    while (n > 0) begin
        n--;
        r = r*x;
    end
    return r;
/*    automatic int k = a - b;
    if (a > b) begin
        k = k * a;
        return k;
    end else if (b > a) begin
        k = k * b;
        return k;
    end else
        return k;
*/
endfunction

//tb program
initial begin
    clk = 1;
    rst = 1;
    #15 end_ready = 1;
    //targeted
    x=2;
    n=3;
    res = gen_result(x, n);
    rst = 0;
    
    //drive dut
    @(posedge clk);
    $display("Test %d: x: %d, n: %d", -1, x, n);
    x_din = x;
    n_din = n;
    start_in = 1;
    start_valid = 1;
    
    @(posedge clk); //wait two cycles for dut to get info
    @(posedge clk);
    x_din = 0;
    n_din = 0;
    start_in = 0;
    start_valid = 0;
    
    @(posedge end_valid);
    dut_res = end_out;
    $display("Test %d: Predicted: %x, Actual: %x, Decimal: %d", -1, res, dut_res, res);
    if (res != dut_res) $error("TEST FAIL");

    @(posedge clk);
    rst = 1;
    @(posedge clk);
     
    test_total = 10;
    for (int test_it=0; test_it < test_total; test_it++) begin
        x=$urandom_range(1,200);
        n=$urandom_range(1,20);
        res = gen_result(x, n);
        rst = 0;
        
        //drive dut
        @(posedge clk);
        $display("Test %d: x: %d, n: %d", test_it, x, n);
        x_din = x;
        n_din = n;
        start_in = 1;
        start_valid = 1;
        
        @(posedge clk); //wait two cycles for dut to get info
        @(posedge clk);
//        x_din = 0;
//        n_din = 0;
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


jpow_graph jp1(
    .clk        (clk        ), 
    .rst        (rst        ), 
    .start_in   (start_in   ),
    .start_valid(start_valid),
    .start_ready(start_ready),
    .end_out    (end_out    ),
    .end_valid  (end_valid  ),
    .end_ready  (end_ready  ),
    .x_din      (x_din      ),
    .x_valid_in (x_valid_in ),
    .x_ready_out(x_ready_out),
    .n_din      (n_din      ),
    .n_valid_in (n_valid_in ),
    .n_ready_out(n_ready_out)
);

endmodule

