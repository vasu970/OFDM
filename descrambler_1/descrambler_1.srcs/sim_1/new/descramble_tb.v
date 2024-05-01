`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 17.02.2024 07:14:27
// Design Name: 
// Module Name: descramble_tb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module descramble_tb( );
 reg clock=0;
    reg enable;
    reg reset;
    reg in_bit=0;
    reg input_strobe=0;
wire out_bit;
     wire output_strobe;
     integer i;
     descramble f1(clock,enable,reset,in_bit,input_strobe,out_bit,output_strobe);
     initial
     begin
     forever
     #5 clock=~clock;
     end
     initial 
     begin
     enable=1;
     input_strobe=1;
     reset=0;
     #1 reset=1;
     #6 reset=0;
     end
     initial
     begin
     for (i=0;i<199;i=i+1)
    #5 in_bit =$random;
     #400
     $finish;
     end
endmodule
