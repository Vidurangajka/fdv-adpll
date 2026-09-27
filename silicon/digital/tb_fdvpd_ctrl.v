// tb_fdvpd_ctrl -- behavioural check of the control block against an ideal
// analog model of the detector.
//
//   iverilog -o tb tb_fdvpd_ctrl.v fdvpd_ctrl.v && vvp tb
//
// The testbench stands in for the analog blocks: it times the ramp window from
// s2_n, decodes the DAC code from dsw when S1 opens, and builds the residue
//     vd = (c + 154) * 781.25 uV - 0.4 mV/ps * dt - sum(switched CDAC weights)
// that the comparator sees.  REF is 50 MHz; CKV runs at the bit-5 fractional
// channel's frequency, phased (+ PHASE_PS) so the gated edge lands at the lock
// point.  For each cycle it checks: the window stayed in 0.3 .. 2.3 ns, and the
// SAR code equals an ideal 7 b quantisation of the residue it was handed.
`timescale 1ps/1fs

module tb;
    parameter real T_REF   = 20000.0;
    parameter real FCW_PD  = 10.0 - 1.0 / 32.0;           // fractional bit 5
    parameter real PHASE_PS = 0.0;                         // phase error
    parameter integer N_CYC = 40;
    localparam real T_CKVD = T_REF / FCW_PD;
    localparam real T_CKV  = T_CKVD / 2.0;
    localparam real LSB_DAC = 0.78125;                     // mV per code
    localparam real SR      = 0.4;                         // mV/ps
    localparam real LSB_ADC = 0.1953125;                   // mV

    reg rst_n = 1, ref_clk = 0, ckv = 0;   // reset is a pulse: flops reset on its falling edge
    wire ckvd, s1_n, s2_n, cmp_clk, code_valid, railed;
    wire [20:0] dsw, dswb;
    wire [6:0] trim, code;
    wire [5:0] bpos, bneg;
    wire [7:0] count_ref;
    reg cmp_p = 0, cmp_n = 0;

    fdvpd_ctrl dut (.rst_n(rst_n), .ref_clk(ref_clk), .ckv(ckv),
        .fcw_int(8'd9), .fcw_frac(16'd63488), .trim_in(7'd64),
        .cmp_p(cmp_p), .cmp_n(cmp_n), .ckvd(ckvd), .dsw(dsw), .dswb(dswb),
        .trim(trim), .s1_n(s1_n), .s2_n(s2_n), .cmp_clk(cmp_clk),
        .bpos(bpos), .bneg(bneg), .code(code), .code_valid(code_valid),
        .railed(railed), .count_ref(count_ref));

    // ---- clocks ---------------------------------------------------------------
    // ckvd rises on odd ckv edges; the gate passes the edge where cnt becomes
    // 10, i.e. ckv edge 19.  REF is placed 300 ps before it (code 0 at reset).
    real t_ckv0 = 1000.0;
    real t_ref0;
    initial begin
        t_ref0 = t_ckv0 + 18.0 * T_CKV - 300.0 + PHASE_PS;
        #10 rst_n = 0;
        #490 rst_n = 1;
    end
    initial begin : ckv_gen
        #(t_ckv0);
        forever begin ckv = 1; #(T_CKV / 2.0); ckv = 0; #(T_CKV / 2.0); end
    end
    initial begin : ref_gen
        #(t_ref0);
        forever begin ref_clk = 1; #(T_REF / 2.0); ref_clk = 0; #(T_REF / 2.0); end
    end

    // ---- the analog stand-in ----------------------------------------------------
    real t_win, dt, vd0;
    integer dac_code, n_ok = 0, n_bad = 0, n_cyc = 0, k;
    function integer dac_of(input [20:0] d);
        integer j, n;
        begin
            n = 0;
            for (j = 0; j < 15; j = j + 1) n = n + 64 * d[j];
            for (j = 0; j < 6; j = j + 1) n = n + (d[20 - j] << j);
            dac_of = n;
        end
    endfunction
    // the DAC code the pair holds is whatever dsw was when S1 opened
    always @(posedge s1_n) dac_code = dac_of(dsw);
    always @(negedge s2_n) t_win = $realtime;
    always @(posedge s2_n) begin
        dt  = $realtime - t_win;
        vd0 = (dac_code + 154) * LSB_DAC - SR * dt;
    end
    function real vd_now(input dummy);
        real v; integer j;
        begin
            v = vd0;
            for (j = 0; j < 6; j = j + 1)
                // real arithmetic: as 1-bit unsigned, 0 - 1 would be +1
                v = v - (1.0 * bpos[j] - 1.0 * bneg[j]) * (1 << j) * LSB_ADC;
            vd_now = v;
        end
    endfunction
    always @(posedge cmp_clk) begin
        #150;
        if (vd_now(0) > 0.0) cmp_p = 1; else cmp_n = 1;
    end
    always @(negedge cmp_clk) begin #100; cmp_p = 0; cmp_n = 0; end

    // ---- checks -------------------------------------------------------------------
    integer ideal;
    always @(posedge code_valid) begin
        n_cyc = n_cyc + 1;
        // ideal offset-binary: residue ~ (code - 63.5) LSB, saturating
        ideal = $rtoi(vd0 / LSB_ADC + 64.0 + 1000.0) - 1000;
        if (ideal < 0) ideal = 0;
        if (ideal > 127) ideal = 127;
        if (n_cyc > 2) begin
            if ((code == ideal || code == ideal - 1 || code == ideal + 1)
                && dt >= 250.0 && dt <= 2400.0) n_ok = n_ok + 1;
            else begin
                n_bad = n_bad + 1;
                $display("cycle %0d: dac %0d dt %0.1f ps vd %0.3f mV -> code %0d, ideal %0d",
                         n_cyc, dac_code, dt, vd0, code, ideal);
            end
        end
        if (n_cyc <= 6 || n_cyc % 8 == 0)
            $display("cycle %0d: dac %4d  dt %7.1f ps  residue %7.3f mV  code %3d  count %0d",
                     n_cyc, dac_code, dt, vd0, code, count_ref);
        if (n_cyc == N_CYC) begin
            $display("RESULT %0d ok, %0d bad of %0d checked", n_ok, n_bad, N_CYC - 2);
            $finish;
        end
    end
    initial begin #(T_REF * (N_CYC + 10)); $display("TIMEOUT"); $finish; end
endmodule
