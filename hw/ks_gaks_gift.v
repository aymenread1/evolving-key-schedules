module gaks_gift_step (input [127:0] si, input [5:0] t, output reg [127:0] so);
  function [3:0] sbox(input [3:0] x);
    case (x)
      4'h0: sbox = 4'h1;
      4'h1: sbox = 4'ha;
      4'h2: sbox = 4'h4;
      4'h3: sbox = 4'hc;
      4'h4: sbox = 4'h6;
      4'h5: sbox = 4'hf;
      4'h6: sbox = 4'h3;
      4'h7: sbox = 4'h9;
      4'h8: sbox = 4'h2;
      4'h9: sbox = 4'hd;
      4'ha: sbox = 4'hb;
      4'hb: sbox = 4'h7;
      4'hc: sbox = 4'h5;
      4'hd: sbox = 4'h0;
      4'he: sbox = 4'h8;
      4'hf: sbox = 4'he;
    endcase
  endfunction
  function [127:0] upd(input [127:0] k);
    upd = {k[31:16] >> 2 | k[31:16] << 14, k[15:0] >> 12 | k[15:0] << 4, k[127:32]};
  endfunction
  reg [127:0] u, s; reg [63:0] h0, l0, dm, h1, l1, h2, l2; wire [5:0] t1 = t + 6'd1;
  always @(*) begin
    u = upd(si); h0 = u[127:64]; l0 = u[63:0]; dm = 0; h1 = 0; l1 = 0; h2 = 0; l2 = 0;
    case (t[0])
      1'd0: begin
        dm = (h0 ^ l0) & 64'h74324e75105437e8; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[2:0], h1[63:3]}; h2 = h1 ^ {l2[37:0], l2[63:38]};
      end
      1'd1: begin
        dm = (h0 ^ l0) & 64'h97a2da2a17ba8f04; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[17:0], h1[63:18]}; h2 = h1 ^ {l2[52:0], l2[63:53]};
      end
    endcase
    s = {h2, l2};
    s[31:28] = sbox(s[31:28]);
    s[67:64] = sbox(s[67:64]);
    s[71:68] = sbox(s[71:68]);
    s[75:72] = sbox(s[75:72]);
    s[91:88] = sbox(s[91:88]);
    so = s ^ {58'd0, t1, 58'd0, t1};
  end
endmodule

module ks_gaks_gift (input clk, input load, input [127:0] key, output [31:0] rk);
  reg [127:0] st; reg [5:0] t; wire [127:0] nx;
  gaks_gift_step u0 (.si(st), .t(t), .so(nx));
  assign rk = st[31:0];
  always @(posedge clk) begin
    if (load) begin st <= key; t <= 6'd0; end
    else if (t != 6'd36) begin st <= nx; t <= t + 6'd1; end
  end
endmodule
