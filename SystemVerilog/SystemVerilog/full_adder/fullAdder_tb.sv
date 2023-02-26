module full_adder_tb();

    timeunit 1ns; // define unit time
    timeprecision 1ps; // 1/1000

    localparam CLK_PERIOD = 10; //clock period
    logic clk;

    //clock signal defining          
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk <= ~clk; //clk!=clk for every 5ns
    end

    logic a     = 0; //using logic can find errors easily (undefiend, infinete impedence,1,0)
    logic b        ; //if use bit there is only 1s' and 0s'
    logic c_in  = 0;
    logic sum   ;
    logic c_out ;


    full_adder fa (.*);
    //full_adder fa(.clk(clk), .a(a), .b(b), .sum(sum), .c_in(c_in), c_out(c_out));

    initial begin
        @(posedge clk);
        #(CLK_PERIOD*3);//delay clock cycle

        @(posedge clk);
        a    <=0;
        b    <=0;
        c_in <=0;
        
        @(posedge clk);
        a    <=0;
        b    <=0;
        c_in <=1;

        @(posedge clk);
        a    <=1;
        b    <=1;
        c_in <=0;

        @(posedge clk);
        a    <=1;
        b    <=1;
        c_in <=1;
    end

endmodule