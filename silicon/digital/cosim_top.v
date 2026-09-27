// cosim_top -- fdvpd_ctrl as ngspice's d_cosim sees it, plus a cycle log.
//
// Same ports as fdvpd_ctrl (d_cosim maps them in declaration order, inputs
// then outputs, vectors MSB first).  Per detector cycle it writes one line to
// the file named by +log=<path> (default cosim_log.txt):
//
//   n  dac  t_win_ps  dt_ps  t_cmp_ps  code  railed
//
// dac is the code the DAC held when S1 last opened, dt the S2 window, t_cmp
// the first comparator clock after it -- where the testbench reads vn - vp.
`timescale 1ps/1fs

module cosim_top (
    input  wire        rst_n,
    input  wire        ref_clk,
    input  wire        ckv,
    input  wire [7:0]  fcw_int,
    input  wire [15:0] fcw_frac,
    input  wire [6:0]  trim_in,
    input  wire        cmp_p,
    input  wire        cmp_n,
    output wire        ckvd,
    output wire [20:0] dsw,
    output wire [20:0] dswb,
    output wire [6:0]  trim,
    output wire        s1_n,
    output wire        s2_n,
    output wire        cmp_clk,
    output wire [5:0]  bpos,
    output wire [5:0]  bneg,
    output wire [5:0]  bmid,
    output wire [6:0]  code,
    output wire        code_valid,
    output wire        railed,
    output wire [7:0]  count_ref
);
    // CKV and REF are generated here when det_cosim.py defines their times
    // (ps; the edges land where the analog sources cross 0.9 V).  Bridged
    // from ngspice, d_cosim dropped a CKV edge that came 28 ps before a
    // comparator edge -- once in 500 -- and the divider slipped a phase for
    // good.  Nothing analog uses either clock, so nothing is lost.
    wire ckv_u, ref_u;
`ifdef TCKV
    reg ckv_i = 1'b0, ref_i = 1'b0;
    real tn, tr;
    integer k, j;
    initial for (k = 0; 1; k = k + 1) begin
        tn = `TCKV0 + k * `TCKV;
        #(tn - $realtime) ckv_i = 1'b1;
        #(tn + `TCKV / 2.0 - $realtime) ckv_i = 1'b0;
    end
    initial for (j = 0; 1; j = j + 1) begin
        tr = `TREF0 + j * `TREF;
        #(tr - $realtime) ref_i = 1'b1;
        #(tr + `TREF / 2.0 - $realtime) ref_i = 1'b0;
    end
    assign ckv_u = ckv_i, ref_u = ref_i;
`else
    assign ckv_u = ckv, ref_u = ref_clk;
`endif

    fdvpd_ctrl u (.rst_n(rst_n), .ref_clk(ref_u), .ckv(ckv_u),
        .fcw_int(fcw_int), .fcw_frac(fcw_frac), .trim_in(trim_in),
        .cmp_p(cmp_p), .cmp_n(cmp_n), .ckvd(ckvd), .dsw(dsw), .dswb(dswb),
        .trim(trim), .s1_n(s1_n), .s2_n(s2_n), .cmp_clk(cmp_clk),
        .bpos(bpos), .bneg(bneg), .bmid(bmid), .code(code), .code_valid(code_valid),
        .railed(railed), .count_ref(count_ref));

    function integer dac_of(input [20:0] d);
        integer j, n;
        begin
            n = 0;
            for (j = 0; j < 15; j = j + 1) n = n + 64 * d[j];
            for (j = 0; j < 6; j = j + 1) n = n + (d[20 - j] << j);
            dac_of = n;
        end
    endfunction

    integer fd, n = 0, dac = 0;
    real t_win = 0.0, dt = 0.0, t_cmp = 0.0;
    reg first = 0;
    reg [1023:0] path;
    initial begin
        if (!$value$plusargs("log=%s", path)) path = "cosim_log.txt";
        fd = $fopen(path, "w");
    end
    always @(posedge s1_n) dac = dac_of(dsw);
    always @(negedge s2_n) t_win = $realtime;
    always @(posedge s2_n) begin dt = $realtime - t_win; first = 1; end
    always @(posedge cmp_clk) if (first) begin t_cmp = $realtime; first = 0; end
    always @(posedge code_valid) begin
        n = n + 1;
        $fdisplay(fd, "%0d %0d %.3f %.3f %.3f %0d %0d", n, dac, t_win, dt, t_cmp,
                  code, railed);
        $fflush(fd);
    end
endmodule
