module nr_sim_tb_top();

logic clk, rst;

//dut ctrl inputs
logic start_in, start_valid;
logic end_ready;

//dut ctrl outputs
logic start_ready;
logic end_valid;
logic [31:0] end_out;

//dut arg inputs
logic [31:0] rts_din;
logic [31:0] x1_din;
logic [31:0] xh_din;
logic rts_valid_in, x1_valid_in, xh_valid_in; //not used in dut

//dut arg outputs
logic rts_ready_out, x1_ready_out, xh_ready_out; //not used in dut

int rts,x1,xh, res;
int test_total;
int dut_res;

function int gen_result;
    input int rts, x1, xh;

    automatic int i = 0;
    automatic int dx;
    while (i < 300) begin
        automatic int df = 4 * rts;
        automatic int f  = 2 * rts * rts - 100;
        i++;
        if (f < 0)
            x1 = rts;
        else
            xh = rts;

        if (((rts - x1) * df - f) * ((rts - xh) * df - f) <= 0) begin

            dx  = (xh - x1) / 2;
            rts = x1 + dx;
        end else begin
            dx = rts / 4;
            rts -= dx;
        end
    end
    return rts;
endfunction

//tb program
initial begin
    clk = 1;
    rst = 1;
    #15 end_ready = 1;
    //targeted
    rts=11;
    x1=11;
    xh=11;
    res = gen_result(rts, x1, xh);
    rst = 0;
    
    //drive dut
    @(posedge clk);
    $display("Test %d: rts: %d, x1: %d, xh: %d", 0, rts, x1, xh);
    rts_din = rts;
    x1_din = x1;
    xh_din = xh;
    start_in = 1;
    start_valid = 1;
    
    @(posedge clk); //wait two cycles for dut to get info
    @(posedge clk);
    //rts_din = 0;
    //x1_din = 0;
    //xh_din = 0;
    start_in = 0;
    start_valid = 0;
    
    @(posedge end_valid);
    dut_res = end_out;
    $display("Test %d: Predicted: %x, Actual: %x, Decimal: %d", 0, res, dut_res, res);
    if (res != dut_res) $error("TEST FAIL");

    @(posedge clk);
    rst = 1;
    @(posedge clk);
     
    //targeted 2
    rts=-7;
    x1=-7;
    xh=-7;
    res = gen_result(rts, x1, xh);
    rst = 0;
    
    //drive dut
    @(posedge clk);
    $display("Test %d: rts: %d, x1: %d, xh: %d", 1, rts, x1, xh);
    rts_din = rts;
    x1_din = x1;
    xh_din = xh;
    start_in = 1;
    start_valid = 1;
    
    @(posedge clk); //wait two cycles for dut to get info
    @(posedge clk);
    //rts_din = 0;
    //x1_din = 0;
    //xh_din = 0;
    start_in = 0;
    start_valid = 0;
    
    @(posedge end_valid);
    dut_res = end_out;
    $display("Test %d: Predicted: %x, Actual: %x, Decimal: %d", 1, res, dut_res, res);
    if (res != dut_res) $error("TEST FAIL");

    @(posedge clk);
    rst = 1;
    @(posedge clk);

    test_total = 10;
    for (int test_it=2; test_it < test_total; test_it++) begin
        rts=$random();
        x1=$random();
        xh=$random();
        res = gen_result(rts, x1, xh);
        rst = 0;
        
        //drive dut
        @(posedge clk);
        $display("Test %d: rts: %d, x1: %d, xh: %d", test_it, rts, x1, xh);
        rts_din = rts;
        x1_din = x1;
        xh_din = xh;
        start_in = 1;
        start_valid = 1;
        
        @(posedge clk); //wait two cycles for dut to get info
        @(posedge clk);
        rts_din = 0;
        x1_din = 0;
        xh_din = 0;
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


newton_raphson_graph jp1(
    .clk            (clk          ), 
    .rst            (rst          ), 
    .start_in       (start_in     ),
    .start_valid    (start_valid  ),
    .start_ready    (start_ready  ),
    .end_out        (end_out      ),
    .end_valid      (end_valid    ),
    .end_ready      (end_ready    ),
    .rts_din        (rts_din      ),
    .rts_valid_in   (rts_valid_in ),
    .rts_ready_out  (rts_ready_out),
    .x1_din         (x1_din       ),
    .x1_valid_in    (x1_valid_in  ),
    .x1_ready_out   (x1_ready_out ),
    .xh_din         (xh_din       ),
    .xh_valid_in    (xh_valid_in  ),
    .xh_ready_out   (xh_ready_out )
);

endmodule

