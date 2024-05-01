/*============================================================================
   module tb_b213.v

   Test bench for (2,1,3) backward label  Viterbi Decoder.
 
===========================================================================*/

`timescale 1 ns/1 ns

module tb_b213;


   reg [`n-1:0] Rx_data [0:37];

   reg [`n-1:0] 	Rx;
   reg                  seq_ready;
   reg 			clock;
   reg 			reset;

   wire [`k-1:0] 	Dx;
   wire [63:0] 		HD;
   wire 		sync_error;
   wire 		oe;

   integer 		i;
   
   
/*===  Instantiate efficient (2,1,3) decoder ===*/
bVITERBI_213 TB_U2 (.Dx(Dx), .oe(oe), .sync_error(sync_error), .Rx(Rx), .seq_ready(seq_ready), .clock(clock), .reset(reset));


/*=== Create Free Running Clock ===*/

always
  begin
     #(0.5*`CLOCK_PERIOD);
     clock=~clock;
  end

/*=== Load Test Vectors ===*/
initial
  begin
     Rx_data[0]  = 2'b00;
     Rx_data[1]  = 2'b11;
     Rx_data[2]  = 2'b10;
     Rx_data[3]  = 2'b10;
     Rx_data[4]  = 2'b11;
     Rx_data[5]  = 2'b01;
     Rx_data[6]  = 2'b10;
     Rx_data[7]  = 2'b00;
     Rx_data[8]  = 2'b00;
     Rx_data[9]  = 2'b01;
     Rx_data[10] = 2'b00;
     Rx_data[11] = 2'b10;
     Rx_data[12] = 2'b11;
     Rx_data[13] = 2'b11;
     Rx_data[14] = 2'b11;
     Rx_data[15] = 2'b10;
     Rx_data[16] = 2'b10;
     Rx_data[17] = 2'b00;
     Rx_data[18] = 2'b00;
     Rx_data[19] = 2'b01;
     Rx_data[20] = 2'b11;
     Rx_data[21] = 2'b11;
     Rx_data[22] = 2'b00;
     Rx_data[23] = 2'b00;
     Rx_data[24] = 2'b11;
     Rx_data[25] = 2'b10;
     Rx_data[26] = 2'b10;
     Rx_data[27] = 2'b11;
     Rx_data[28] = 2'b01;
     Rx_data[29] = 2'b10;
     Rx_data[30] = 2'b00;
     Rx_data[31] = 2'b00;
     Rx_data[32] = 2'b01;
     Rx_data[33] = 2'b00;
     Rx_data[34] = 2'b10;
     Rx_data[35] = 2'b11;
     Rx_data[36] = 2'b11;
     Rx_data[37] = 2'b11;
     
   end

/*=== Clock data into the decoder ===*/
initial
  begin
     clock    = 0;
     reset    = 1;            //Apply reset
     #`CLOCK_PERIOD reset =0; //Disable reset after one clock period

     seq_ready =1;

     for (i=0; i<15; i=i+1)
       begin
               Rx=Rx_data[i];
               #(`CLOCK_PERIOD*3);

        end

     // After 1st Traceback, clock in remaining bits after each traceback
     for (i=15; i<33; i=i+1)
       begin
               Rx=Rx_data[i];
               #(`CLOCK_PERIOD*19);

        end

     seq_ready =0;
     
  end
   

/*=== Create Simulation Waveforms ===*/

//initial
//  begin
//     $fsdbDumpfile("verilog.fsdb");// debussy dump format
//     $fsdbDumpStrength;
//     $fsdbDumpvars;
//  end
   

initial
   #27000 $finish;


/*=== Monitor Viterbi Decoder Responses ===*/
initial
   $monitor("time=%4g, clock=%b, Rx=%b, we=%b, te=%b, oe=%b, wptr=%h, tptr=%h, min_state=%b, tb_reg=%b,A_in=%h,P_in=%b,plabels=%b,Dx=%b, ns=%b",
             $time,
             clock,
	     Rx,
	     TB_U2.U3.we,
	     TB_U2.U3.te,
	     TB_U2.U3.oe,	    
	     TB_U2.U3.write_ptr,
    	     TB_U2.U3.trace_ptr,
             TB_U2.U3.min_state,
	    
             TB_U2.U3.tb_reg,

            {TB_U2.U3.A0_in, 
             TB_U2.U3.A1_in,
             TB_U2.U3.A2_in,
             TB_U2.U3.A3_in,
             TB_U2.U3.A4_in,
             TB_U2.U3.A5_in,
             TB_U2.U3.A6_in,
             TB_U2.U3.A7_in},
	     	     
             {TB_U2.U3.P0_in, 
             TB_U2.U3.P1_in,
             TB_U2.U3.P2_in,
             TB_U2.U3.P3_in,
             TB_U2.U3.P4_in,
             TB_U2.U3.P5_in,
             TB_U2.U3.P6_in,
             TB_U2.U3.P7_in},
	    

            {TB_U2.U3.P0[`T-1], 
             TB_U2.U3.P1[`T-1],
             TB_U2.U3.P2[`T-1],
             TB_U2.U3.P3[`T-1],
             TB_U2.U3.P4[`T-1],
             TB_U2.U3.P5[`T-1],
             TB_U2.U3.P6[`T-1],
             TB_U2.U3.P7[`T-1]},
	    
             TB_U2.U3.Dx,
	     TB_U2.U3.NEXT_STATE
             );


endmodule





