module ks_gift (input clk, input load, input [127:0] key, output [31:0] rk);
  function [127:0] upd(input [127:0] k);
    upd = {k[31:16] >> 2 | k[31:16] << 14, k[15:0] >> 12 | k[15:0] << 4, k[127:32]};
  endfunction
  reg [127:0] st;
  assign rk = st[31:0];
  always @(posedge clk) st <= load ? key : upd(st);
endmodule
