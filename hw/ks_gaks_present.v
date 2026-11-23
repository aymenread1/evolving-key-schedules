module gaks_present_step (input [79:0] si, input [5:0] t, output reg [79:0] so);
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

  reg [79:0] r61, rc; reg [39:0] h0, l0, dm, h1, l1, h2, l2;
  always @(*) begin
    r61 = {si[18:0], si[79:19]};
    h0 = r61[79:40]; l0 = r61[39:0];
    dm = 0; h1 = 0; l1 = 0; h2 = 0; l2 = 0;
    case (t[1:0])
      2'd0: begin
        dm = (h0 ^ l0) & 40'hbe3ced35a8; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[35:0], h1[39:36]}; h2 = h1 ^ {l2[19:0], l2[39:20]};
      end
      2'd1: begin
        dm = (h0 ^ l0) & 40'hfabb061e50; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[29:0], h1[39:30]}; h2 = h1 ^ {l2[30:0], l2[39:31]};
      end
      2'd2: begin
        dm = (h0 ^ l0) & 40'hd3fc7220bf; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[37:0], h1[39:38]}; h2 = h1 ^ {l2[23:0], l2[39:24]};
      end
      2'd3: begin
        dm = (h0 ^ l0) & 40'h3ce448da4c; h1 = h0 ^ dm; l1 = l0 ^ dm;
        l2 = l1 ^ {h1[35:0], h1[39:36]}; h2 = h1 ^ {l2[19:0], l2[39:20]};
      end
    endcase
    h2[3:0] = sbox(h2[3:0]);
    h2[39:36] = sbox(h2[39:36]);
    l2[27:24] = sbox(l2[27:24]);
    case (t)
      6'd0: rc = 80'h50df1e641a53336f5fa1;
      6'd1: rc = 80'h99ae30e9a35002428f8d;
      6'd2: rc = 80'hd526939617b0ae2f6fe6;
      6'd3: rc = 80'hfba8f36a1a0437ac6ad6;
      6'd4: rc = 80'h25101e4607a912a97ad3;
      6'd5: rc = 80'h504220484c610a6abc87;
      6'd6: rc = 80'hd8edea53894685c571b4;
      6'd7: rc = 80'h1ef80ad049bcc327515f;
      6'd8: rc = 80'h8d9626f56d938b2668c3;
      6'd9: rc = 80'hcfe1cf25476129203acf;
      6'd10: rc = 80'h4b8951aa8d596d3917d3;
      6'd11: rc = 80'h4de2bc3e09cb85ccf2ab;
      6'd12: rc = 80'h87fdfc2764cce71e2e52;
      6'd13: rc = 80'h835d247098e6c317bb33;
      6'd14: rc = 80'h3a8d494c812700d969b6;
      6'd15: rc = 80'haee4b7f29e1bf4855b0b;
      6'd16: rc = 80'h1692f7122ffa1860c5c8;
      6'd17: rc = 80'h69f1bd440fb2535d8830;
      6'd18: rc = 80'h005aa944a7784bef434e;
      6'd19: rc = 80'h63cfe2e819005cfbfcf7;
      6'd20: rc = 80'h3efe85120bf122d0c318;
      6'd21: rc = 80'hc37d2739405f26fc77c5;
      6'd22: rc = 80'h0127f1076b2ad0b969e5;
      6'd23: rc = 80'h5fbbfe707caf7a08a60b;
      6'd24: rc = 80'h2124d8a94da646867457;
      6'd25: rc = 80'h55943d388a8477d48dd5;
      6'd26: rc = 80'heb69dee134c23cbd8cdb;
      6'd27: rc = 80'h70be4f64ef1b338673a8;
      6'd28: rc = 80'h219028e7e5127bfa3089;
      6'd29: rc = 80'hb8bdc3e4ee49afff5590;
      6'd30: rc = 80'h24cd2f360fbcef93ec69;
      6'd31: rc = 80'ha45c774865009c1be6e4;
      6'd32: rc = 80'h9bcce322b19fbc10d256;
      6'd33: rc = 80'hc7d4f3ffa7997d15953b;
      6'd34: rc = 80'hb88e8ead02050456c86f;
      6'd35: rc = 80'hd2282d369a2f0302bcfb;
      6'd36: rc = 80'h3c423722db550bec70ff;
      6'd37: rc = 80'ha673da8d850baf820c1e;
      default: rc = 80'h0;
    endcase
    so = {h2, l2} ^ rc;
  end
endmodule

module ks_gaks_present (input clk, input load, input [79:0] key, output [63:0] rk);
  reg [79:0] st; reg [5:0] t; wire [79:0] nx;
  gaks_present_step u (.si(st), .t(t), .so(nx));
  assign rk = st[79:16];
  always @(posedge clk) begin
    if (load) begin st <= key; t <= 6'd0; end
    else if (t != 6'd38) begin st <= nx; t <= t + 6'd1; end
  end
endmodule
