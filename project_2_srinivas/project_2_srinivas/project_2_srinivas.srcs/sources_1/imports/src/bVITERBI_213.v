/*============================================================================
   module bVITERBI_213.v

   Top Level Module for the (2,1,3) backward label  Viterbi Decoder.

   This Module Contains:-
    - Branch Metric Unit (BMU)
    - Add Compare Select Unit (ACSU)
    - Metric and path Memory Units
    - Synchronization and Control Unit
   
  ===========================================================================*/


`timescale 1 ns/1 ns

module bVITERBI_213 (
,           // Decoded output signal
                     oe,           // Dx output enable
                     sync_error,   // Sync error output signal
		     
                     Rx,           // Input bits
		     seq_ready,    // Ready signal to decode next sequence
                     clock,
                     reset);


reg input_strobe=1;

   output [`k-1:0] Dx;
   output          oe;
   output          sync_error;
   
   input [`n-1:0]  Rx;
   input           seq_ready;
   input           reset;
   input           clock;
   
   wire [`W-1:0] A0_out;
   wire [`W-1:0] A1_out;
   wire [`W-1:0] A2_out;
   wire [`W-1:0] A3_out;
   wire [`W-1:0] A4_out;
   wire [`W-1:0] A5_out;
   wire [`W-1:0] A6_out;
   wire [`W-1:0] A7_out;   
		

   wire [`W-1:0] acs0_ppm_out;
   wire [`W-1:0] acs1_ppm_out;
   wire [`W-1:0] acs2_ppm_out;
   wire [`W-1:0] acs3_ppm_out;
   wire [`W-1:0] acs4_ppm_out;
   wire [`W-1:0] acs5_ppm_out;   
   wire [`W-1:0] acs6_ppm_out;   
   wire [`W-1:0] acs7_ppm_out;
   
   
   wire [`k-1:0] acs0_Bx_out;
   wire [`k-1:0] acs1_Bx_out;
   wire [`k-1:0] acs2_Bx_out;
   wire [`k-1:0] acs3_Bx_out;
   wire [`k-1:0] acs4_Bx_out;
   wire [`k-1:0] acs5_Bx_out;
   wire [`k-1:0] acs6_Bx_out;
   wire [`k-1:0] acs7_Bx_out;   

   wire [1:0] 	 HD1;
   wire [1:0] 	 HD2;
   wire [1:0] 	 HD3;
   wire [1:0] 	 HD4;
   wire [1:0] 	 HD5;
   wire [1:0] 	 HD6;
   wire [1:0] 	 HD7;
   wire [1:0] 	 HD8;
   wire [1:0] 	 HD9;
   wire [1:0] 	 HD10;
   wire [1:0] 	 HD11;
   wire [1:0] 	 HD12;
   wire [1:0] 	 HD13;
   wire [1:0] 	 HD14;
   wire [1:0] 	 HD15;
   wire [1:0] 	 HD16;  

/*=== Instantiate the BMU === */
bBMU_213 U1 (.HD1(HD1),   .HD2(HD2),   .HD3(HD3),   .HD4(HD4),   
             .HD5(HD5),   .HD6(HD6),   .HD7(HD7),   .HD8(HD8),   
             .HD9(HD9),   .HD10(HD10), .HD11(HD11), .HD12(HD12), 
             .HD13(HD13), .HD14(HD14), .HD15(HD15), .HD16(HD16), 
             .Rx(Rx), .le(le), .reset(reset), .clock(clock)
            );

/*=== Instantiate ACSU ===*/


bACSU_213 U2 (.acs0_ppm_out(acs0_ppm_out), .acs1_ppm_out(acs1_ppm_out), 
              .acs2_ppm_out(acs2_ppm_out), .acs3_ppm_out(acs3_ppm_out), 
              .acs4_ppm_out(acs4_ppm_out), .acs5_ppm_out(acs5_ppm_out), 
              .acs6_ppm_out(acs6_ppm_out), .acs7_ppm_out(acs7_ppm_out),

              .acs0_Bx_out(acs0_Bx_out),   .acs1_Bx_out(acs1_Bx_out), 
              .acs2_Bx_out(acs2_Bx_out),   .acs3_Bx_out(acs3_Bx_out), 
              .acs4_Bx_out(acs4_Bx_out),   .acs5_Bx_out(acs5_Bx_out), 
              .acs6_Bx_out(acs6_Bx_out),   .acs7_Bx_out(acs7_Bx_out),

              .acs0_ppm_ina(A0_out),  .acs0_ppm_inb(A4_out), .HD0_ina(HD1), .HD0_inb(HD9),
              .acs1_ppm_ina(A0_out),  .acs1_ppm_inb(A4_out), .HD1_ina(HD2), .HD1_inb(HD10),
              .acs2_ppm_ina(A1_out),  .acs2_ppm_inb(A5_out), .HD2_ina(HD3), .HD2_inb(HD11),
              .acs3_ppm_ina(A1_out),  .acs3_ppm_inb(A5_out), .HD3_ina(HD4), .HD3_inb(HD12),
              .acs4_ppm_ina(A2_out),  .acs4_ppm_inb(A6_out), .HD4_ina(HD5), .HD4_inb(HD13),
              .acs5_ppm_ina(A2_out),  .acs5_ppm_inb(A6_out), .HD5_ina(HD6), .HD5_inb(HD14),
              .acs6_ppm_ina(A3_out),  .acs6_ppm_inb(A7_out), .HD6_ina(HD7), .HD6_inb(HD15),
              .acs7_ppm_ina(A3_out),  .acs7_ppm_inb(A7_out), .HD7_ina(HD8), .HD7_inb(HD16),
	      .ae(ae), .reset(reset), .clock(clock)
              );


/*=== Instantiate Control Unit  ===*/

bCONTROL_213 U3 (.Dx(Dx), .oe(oe), .le(le), .ae(ae), .sync_error(sync_error), 
                 .A0_out(A0_out), .A1_out(A1_out), .A2_out(A2_out), .A3_out(A3_out), 
                 .A4_out(A4_out), .A5_out(A5_out), .A6_out(A6_out), .A7_out(A7_out),
                 
                 .A0_in(acs0_ppm_out), .A1_in(acs1_ppm_out), .A2_in(acs2_ppm_out), .A3_in(acs3_ppm_out),  
                 .A4_in(acs4_ppm_out), .A5_in(acs5_ppm_out), .A6_in(acs6_ppm_out), .A7_in(acs7_ppm_out),
                 
                 .P0_in(acs0_Bx_out), .P1_in(acs1_Bx_out), .P2_in(acs2_Bx_out), .P3_in(acs3_Bx_out), 
                 .P4_in(acs4_Bx_out), .P5_in(acs5_Bx_out), .P6_in(acs6_Bx_out), .P7_in(acs7_Bx_out),
                 .seq_ready(seq_ready), .reset(reset), .clock(clock)
                );
descramble f1(
    clock,
    oe,
    reset,

    Dx,
     input_strobe,

   out_bit,
     output_strobe
);

endmodule

