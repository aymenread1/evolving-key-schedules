module gaks_des_hw_step (input [55:0] si, input [4:0] t, output reg [55:0] so);
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

  reg [27:0] ci, di, c1, d1, dm, c2, d2, c3, d3; reg [55:0] s; reg two; wire [4:0] t1 = t + 5'd1;
  always @(*) begin
    case (t)
      5'd0: two = 1'b0;
      5'd1: two = 1'b0;
      5'd2: two = 1'b0;
      5'd3: two = 1'b0;
      5'd4: two = 1'b0;
      5'd5: two = 1'b0;
      5'd6: two = 1'b0;
      5'd7: two = 1'b0;
      5'd8: two = 1'b0;
      5'd9: two = 1'b1;
      5'd10: two = 1'b1;
      5'd11: two = 1'b1;
      5'd12: two = 1'b1;
      5'd13: two = 1'b1;
      5'd14: two = 1'b1;
      5'd15: two = 1'b0;
      5'd16: two = 1'b1;
      5'd17: two = 1'b1;
      5'd18: two = 1'b1;
      5'd19: two = 1'b1;
      5'd20: two = 1'b1;
      5'd21: two = 1'b1;
      5'd22: two = 1'b0;
      default: two = 1'b0;
    endcase
    ci = si[55:28]; di = si[27:0];
    c1 = two ? {ci[25:0], ci[27:26]} : {ci[26:0], ci[27:27]};
    d1 = two ? {di[25:0], di[27:26]} : {di[26:0], di[27:27]};
    dm = 0; c2 = 0; d2 = 0; c3 = 0; d3 = 0;
    case (t[0:0])
      1'd0: begin
        dm = (c1 ^ d1) & 28'h3f05edf; c2 = c1 ^ dm; d2 = d1 ^ dm;
        d3 = d2 ^ {c2[16:0], c2[27:17]}; c3 = c2 ^ {d3[6:0], d3[27:7]};
      end
      1'd1: begin
        dm = (c1 ^ d1) & 28'hbc2d763; c2 = c1 ^ dm; d2 = d1 ^ dm;
        d3 = d2 ^ {c2[22:0], c2[27:23]}; c3 = c2 ^ {d3[19:0], d3[27:20]};
      end
    endcase
    s = {c3, d3};
    s[23:20] = sbox(s[23:20]);
    s[27:24] = sbox(s[27:24]);
    so = s ^ {23'd0, t1, 23'd0, t1};
  end
endmodule

module ks_gaks_des_hw (input clk, input load, input [63:0] key, output [47:0] rk);
  reg [55:0] st; reg [4:0] t; wire [55:0] k56, nx;
  assign k56[55] = key[7];
  assign k56[54] = key[15];
  assign k56[53] = key[23];
  assign k56[52] = key[31];
  assign k56[51] = key[39];
  assign k56[50] = key[47];
  assign k56[49] = key[55];
  assign k56[48] = key[63];
  assign k56[47] = key[6];
  assign k56[46] = key[14];
  assign k56[45] = key[22];
  assign k56[44] = key[30];
  assign k56[43] = key[38];
  assign k56[42] = key[46];
  assign k56[41] = key[54];
  assign k56[40] = key[62];
  assign k56[39] = key[5];
  assign k56[38] = key[13];
  assign k56[37] = key[21];
  assign k56[36] = key[29];
  assign k56[35] = key[37];
  assign k56[34] = key[45];
  assign k56[33] = key[53];
  assign k56[32] = key[61];
  assign k56[31] = key[4];
  assign k56[30] = key[12];
  assign k56[29] = key[20];
  assign k56[28] = key[28];
  assign k56[27] = key[1];
  assign k56[26] = key[9];
  assign k56[25] = key[17];
  assign k56[24] = key[25];
  assign k56[23] = key[33];
  assign k56[22] = key[41];
  assign k56[21] = key[49];
  assign k56[20] = key[57];
  assign k56[19] = key[2];
  assign k56[18] = key[10];
  assign k56[17] = key[18];
  assign k56[16] = key[26];
  assign k56[15] = key[34];
  assign k56[14] = key[42];
  assign k56[13] = key[50];
  assign k56[12] = key[58];
  assign k56[11] = key[3];
  assign k56[10] = key[11];
  assign k56[9] = key[19];
  assign k56[8] = key[27];
  assign k56[7] = key[35];
  assign k56[6] = key[43];
  assign k56[5] = key[51];
  assign k56[4] = key[59];
  assign k56[3] = key[36];
  assign k56[2] = key[44];
  assign k56[1] = key[52];
  assign k56[0] = key[60];
  assign rk[47] = st[42];
  assign rk[46] = st[39];
  assign rk[45] = st[45];
  assign rk[44] = st[32];
  assign rk[43] = st[55];
  assign rk[42] = st[51];
  assign rk[41] = st[53];
  assign rk[40] = st[28];
  assign rk[39] = st[41];
  assign rk[38] = st[50];
  assign rk[37] = st[35];
  assign rk[36] = st[46];
  assign rk[35] = st[33];
  assign rk[34] = st[37];
  assign rk[33] = st[44];
  assign rk[32] = st[52];
  assign rk[31] = st[30];
  assign rk[30] = st[48];
  assign rk[29] = st[40];
  assign rk[28] = st[49];
  assign rk[27] = st[29];
  assign rk[26] = st[36];
  assign rk[25] = st[43];
  assign rk[24] = st[54];
  assign rk[23] = st[15];
  assign rk[22] = st[4];
  assign rk[21] = st[25];
  assign rk[20] = st[19];
  assign rk[19] = st[9];
  assign rk[18] = st[1];
  assign rk[17] = st[26];
  assign rk[16] = st[16];
  assign rk[15] = st[5];
  assign rk[14] = st[11];
  assign rk[13] = st[23];
  assign rk[12] = st[8];
  assign rk[11] = st[12];
  assign rk[10] = st[7];
  assign rk[9] = st[17];
  assign rk[8] = st[0];
  assign rk[7] = st[22];
  assign rk[6] = st[3];
  assign rk[5] = st[10];
  assign rk[4] = st[14];
  assign rk[3] = st[6];
  assign rk[2] = st[20];
  assign rk[1] = st[27];
  assign rk[0] = st[24];
  gaks_des_hw_step u (.si(st), .t(t), .so(nx));
  always @(posedge clk) begin
    if (load) begin st <= k56; t <= 5'd0; end
    else if (t != 5'd23) begin st <= nx; t <= t + 5'd1; end
  end
endmodule
