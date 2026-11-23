module ks_present (input clk, input load, input [79:0] key, output [63:0] rk);
  function [3:0] sbox(input [3:0] x);
    case (x)
      4'h0: sbox = 4'hc;
      4'h1: sbox = 4'h5;
      4'h2: sbox = 4'h6;
      4'h3: sbox = 4'hb;
      4'h4: sbox = 4'h9;
      4'h5: sbox = 4'h0;
      4'h6: sbox = 4'ha;
      4'h7: sbox = 4'hd;
      4'h8: sbox = 4'h3;
      4'h9: sbox = 4'he;
      4'ha: sbox = 4'hf;
      4'hb: sbox = 4'h8;
      4'hc: sbox = 4'h4;
      4'hd: sbox = 4'h7;
      4'he: sbox = 4'h1;
      4'hf: sbox = 4'h2;
    endcase
  endfunction

  reg [79:0] st; reg [4:0] i;
  wire [79:0] r61 = {st[18:0], st[79:19]};
  wire [79:0] s1 = {sbox(r61[79:76]), r61[75:0]};
  wire [79:0] nx = {s1[79:20], s1[19:15] ^ (i + 5'd1), s1[14:0]};
  assign rk = st[79:16];
  always @(posedge clk) begin
    if (load) begin st <= key; i <= 5'd0; end
    else begin st <= nx; i <= i + 5'd1; end
  end
endmodule
