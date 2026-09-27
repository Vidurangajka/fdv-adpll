// fdvpd_ctrl -- timing and control for the sky130 FDVPD.
//
// Everything the detector needs besides its analog blocks:
//
//   * CKVd = CKV / 2, and a free-running CKVd edge counter
//   * the reference accumulator, advanced once per cycle at the gated edge
//   * the ramp window: set by REF, ended by the CKVd edge the accumulator
//     names -- by ABSOLUTE COUNT, not by a time window.  A time-threshold
//     selection races at every fractional wrap: the same edge is the target
//     for one code and "too early" for the next, and whichever way the
//     threshold falls, the DAC and the gate disagree and the SAR rails once
//     per sawtooth period.  The count is exact; only acquisition can put the
//     gated edge outside the window, and a timeout covers that.
//   * a self-timed SAR: the comparator's own "decided" signal clocks the next
//     step (a 1 GHz synchronous SAR does not close in sky130 hd at ss)
//   * the encode window for the DAC (S1) and the DAC switch drives
//   * the CKVd count sampled at REF, for the frequency-lock path
//
// Detector cycle, times after REF (dt = the ramp window, 0.3 .. 2.3 ns):
//   0            win rises: the ramp current steers into vn
//   dt           the gated CKVd edge: win falls, the current steers back to VDD
//   dt + 2       seq stage 0, the next CKVd edge; the SAR runs (7
//                comparisons, ~5.3 ns)
//   stage 1      the DAC takes the next cycle's code (vp is disconnected)
//   stage 3      the SAR's code is taken; S1 closes: vp settles to the DAC,
//                vn to VDD, the CDAC plates return to mid
//   stage 7      S1 opens, dt + 16 ns: at least 1.7 ns before the next
//                REF.  8 ns of encode is ~6 tau of the post-layout DAC.
//
// Only the CKVd divider runs at CKV.  The sequencer is a one-hot shift
// register on CKVd: nothing between flops, and 2 ns a stage, so it closes at
// ss 1.60 V 100 C where a 1 GHz one does not (a hd flop's clk-to-Q and setup
// alone are ~1.2 ns there).
//
// Conventions shared with the analog blocks:
//   s1_n, s2_n  PMOS switch gates, low = closed
//   cmp_clk     comparator clock, high = evaluate.  cmp_p / cmp_n are its
//               outputs after inverters: both low in reset, one rises when it
//               decides; cmp_p = 1 means vn > vp (the residue is positive)
//   bpos[k]     weight 2^k: p plate to H, n plate to L (lowers vn - vp)
//   bneg[k]     weight 2^k: p plate to L, n plate to H
//   bmid[k]     neither: both plates to mid.  Decoded here rather than by a
//               NOR per weight in the analog cell, which then holds only
//               switches, capacitors and the comparator
//   dsw[20:0]   DAC switches t0..t14, b5..b0: 1 = the segment feeds outp
//
// `define SIM uses behavioural models for the three hard cells (a clock-gate
// latch and two delay lines); otherwise the sky130 cells are instantiated.

`timescale 1ps/1fs

module fdvpd_ctrl (
    input  wire        rst_n,
    input  wire        ref_clk,
    input  wire        ckv,
    input  wire [7:0]  fcw_int,      // integer part of FCW / 2 (9 or 10 here)
    input  wire [15:0] fcw_frac,     // fractional part, 2^-16 of a CKVd cycle
    input  wire [6:0]  trim_in,
    input  wire        cmp_p,
    input  wire        cmp_n,
    output wire        ckvd,
    output reg  [20:0] dsw,
    output reg  [20:0] dswb,
    output wire [6:0]  trim,
    output wire        s1_n,
    output wire        s2_n,
    output wire        cmp_clk,
    output wire [5:0]  bpos,
    output wire [5:0]  bneg,
    output wire [5:0]  bmid,
    output reg  [6:0]  code,
    output reg         code_valid,
    output reg         railed,
    output reg  [7:0]  count_ref
);
    assign trim = trim_in;

    // ---- CKVd and its free-running counter ------------------------------------
    // The one flop at CKV.  At ss 1.60 V 100 C a reset flop and an inverter
    // take 1.04 ns; without the reset the loop fits, just: OpenSTA with ideal
    // wiring gives dfxbp_1 Q_N -> D +10 ps and dfxtp_2 -> inv_1 +20 ps, the
    // best of the library.  Its phase is arbitrary and nothing depends on it,
    // so it has no reset -- the behavioural model keeps one so that RTL
    // simulation starts known.
`ifdef SIM
    reg ckvd_r;
    always @(posedge ckv or negedge rst_n)
        if (!rst_n) ckvd_r <= 1'b0; else ckvd_r <= ~ckvd_r;
    assign ckvd = ckvd_r;
`else
    // The divider's loop is Q -> inverter -> D, and Q drives nothing else: at
    // ss 1.60 V 100 C clk-to-Q and setup alone nearly fill the 1 ns.  CKVd is
    // buffered off the inverter's output (either phase is a /2 clock), and
    // the generated clock is declared on that buffer (fdvpd_ctrl.sdc), which
    // keeps CTS out of the loop -- declared on Q, CTS buffered the feedback
    // and the divider missed by 0.28 ns.  The resizer is kept off both nets
    // (config.json RSZ_DONT_TOUCH_RX).  Even so a hd flop does not quite
    // toggle at 1 GHz there; the real fix is a custom TSPC stage in the DCO.
    wire ckvd_q, ckvd_n;
    sky130_fd_sc_hd__dfxtp_2 u_div (.CLK(ckv), .D(ckvd_n), .Q(ckvd_q));
    sky130_fd_sc_hd__inv_2 u_divi (.A(ckvd_q), .Y(ckvd_n));
    sky130_fd_sc_hd__clkbuf_2 u_ckvd_buf (.A(ckvd_n), .X(ckvd));
`endif

    // cnt counts CKVd edges, in four 2 b digits, each stepping on a carry
    // registered one edge early.  At ss 1.60 V 100 C a reset flop's clk-to-Q
    // is 0.83 ns of the 2 ns, and an 8 b increment (or even a 4 b one with its
    // enable) does not fit in the rest.  Same count.
    reg [1:0] c0, c1, c2, c3;
    reg       cy1, cy2, cy3, f1, f12;
    wire [7:0] cnt = {c3, c2, c1, c0};
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) begin
            c0 <= 2'd0; c1 <= 2'd0; c2 <= 2'd0; c3 <= 2'd0;
            cy1 <= 1'b0; cy2 <= 1'b0; cy3 <= 1'b0; f1 <= 1'b0; f12 <= 1'b0;
        end else begin
            c0  <= c0 + 2'd1;
            // c1 and c2 only move when c0 wraps, so where c0 is 2 they have
            // held since the edge before: their all-ones flags are registered
            f1  <= (c1 == 2'd3);
            f12 <= (c1 == 2'd3) & (c2 == 2'd3);
            // c0 is 2 now, so it wraps at the next edge: carry then
            cy1 <= (c0 == 2'd2);
            cy2 <= (c0 == 2'd2) & f1;
            cy3 <= (c0 == 2'd2) & f12;
            // enabled 2 b increments as XORs: written as `if (cy) c <= c + 1`
            // they synthesised to a mux that missed by 4 ps at ss
            c1 <= c1 ^ {cy1 & c1[0], cy1};
            c2 <= c2 ^ {cy2 & c2[0], cy2};
            c3 <= c3 ^ {cy3 & c3[0], cy3};
        end

    // ---- accumulator: the gated edge is number ceil(Phi_R) ----------------------
    // floor(Phi_R), plus one whenever the fractional accumulator is non-zero
    reg  [7:0]  n_floor;
    reg  [15:0] acc;
    reg  [9:0]  dac_next;
    wire [7:0]  n_target = n_floor + {7'd0, acc != 16'd0};
    // en = (cnt == n_target - 1), pipelined so the clock gate's enable comes
    // straight off a flop (the compare in one CKVd cycle missed by 2 ns at
    // ss): n_target - 3 is registered on CKVd, each 2 b digit is matched
    // against it, and the matches are ANDed.  Registering t3 keeps synthesis
    // from folding the count into the accumulator's arithmetic.  n_target
    // moves at the gated edge and may still be settling for the next two
    // CKVd edges (a false path: it has ~18 ns), so t3_r can hold a transient
    // for two edges and m for three: matches there are masked.  A real one is
    // never that close -- FCW / 2 is at least 9.
    wire [7:0]  t3 = n_target - 8'd3;
    reg  [7:0]  t3_r;
    reg  [3:0]  m;
    reg         en;
    reg  [7:0]  seq;
    always @(posedge ckvd or negedge rst_n)
        // t3_r out of reset is n_target - 3 for the reset n_floor (10): any
        // other value is an early match, and a stop there shifts every
        // later target
        if (!rst_n) begin t3_r <= 8'd7; m <= 4'd0; en <= 1'b0; end
        else begin
            t3_r <= t3;
            m    <= {c3 == t3_r[7:6], c2 == t3_r[5:4], c1 == t3_r[3:2],
                     c0 == t3_r[1:0]};
            en   <= (&m) & ~seq[0] & ~seq[1] & ~seq[2];
        end

    // ---- the gate and the ramp window ---------------------------------------------
    wire gclk;
`ifdef SIM
    reg en_l;
    always @(ckvd or en) if (!ckvd) en_l = en;
    assign gclk = ckvd & en_l;
`else
    sky130_fd_sc_hd__dlclkp_1 u_icg (.CLK(ckvd), .GATE(en), .GCLK(gclk));
`endif

    reg  timeout;
    wire stop_async;
`ifdef SIM
    assign stop_async = gclk | timeout;
`else
    // clocks the 34 accumulator flops: yosys's or2_0 gave 2.9 ns edges at ss
    sky130_fd_sc_hd__or2_4 u_stop (.A(gclk), .B(timeout), .X(stop_async));
`endif
    reg  win;
    always @(posedge ref_clk or posedge stop_async or negedge rst_n)
        if (!rst_n)          win <= 1'b0;
        else if (stop_async) win <= 1'b0;
        else                 win <= 1'b1;
    // s2_n steers the ramp current: into vn while low, into VDD otherwise
    // (the steering pair and its drivers are in the front end)
    assign s2_n = ~win;

    // ---- the next cycle's code, at each gated edge ------------------------------
    wire [16:0] acc_sum = {1'b0, acc} + {1'b0, fcw_frac};
    wire [16:0] tfrac   = 17'h10000 - {1'b0, acc_sum[15:0]};   // ceil - frac
    wire [10:0] tfrac_r = (tfrac + 17'd32) >> 6;                 // rounded, /64
    always @(posedge stop_async or negedge rst_n)
        if (!rst_n) begin
            acc <= 16'd0; n_floor <= 8'd10; dac_next <= 10'd0;
        end else begin
            acc      <= acc_sum[15:0];
            n_floor  <= n_floor + fcw_int + {7'd0, acc_sum[16]};
            // round to 10 b; a whole period (tfrac = 1) is code 0
            dac_next <= (acc_sum[15:0] == 16'd0) ? 10'd0
                      : (tfrac_r > 11'd1023 ? 10'd1023 : tfrac_r[9:0]);
        end

    // ---- CKVd domain: synchronise the window, time out, sequence ----------------
    // Everything after the gate runs on CKVd.  At CKV a hd flop's clk-to-Q and
    // setup alone fill the 1 ns period at ss 1.60 V 100 C; at 2 ns per stage
    // the stage times -- the SAR's 6 ns, the 8 ns encode -- are unchanged.
    reg [1:0] wsync;
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) wsync <= 2'b00; else wsync <= {wsync[0], win};
    // The window's FALLING edge is a gated CKVd edge plus the gate and the
    // flop, so win is synchronous when it falls and the stop is seen at the
    // next CKVd edge: wsync[0] still holds the window, win no longer does.
    // Its rising edge is REF, asynchronous, but win can only rise while
    // wsync[0] is low, where stopped ignores it; only the timeout (on
    // wsync[1]) looks at that edge.
    wire stopped = wsync[0] & ~win;

    // no gated edge within 3 CKVd cycles of the window opening (a real window
    // is <= 2.3 ns, so it spans at most two CKVd edges): end it
    reg [2:0] wage;
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) begin wage <= 3'd0; timeout <= 1'b0; end
        else begin
            wage    <= wsync[1] ? {wage[1:0], 1'b1} : 3'd0;
            timeout <= &wage;
        end

    // one-hot: seq[i] is high in the i-th CKVd cycle after the window closed
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) seq <= 8'd0; else seq <= {seq[6:0], stopped};

    reg s1_closed, sar_on;
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) begin s1_closed <= 1'b1; sar_on <= 1'b0; end
        else begin
            if (stopped)  s1_closed <= 1'b0;
            if (seq[3])   s1_closed <= 1'b1;
            if (seq[7])   s1_closed <= 1'b0;
            if (stopped)  sar_on <= 1'b1;
            if (seq[3])   sar_on <= 1'b0;
        end
    assign s1_n = ~s1_closed;

    // DAC code: thermometer + binary, loaded at stage 1
    wire [3:0]  msb = dac_next[9:6];
    wire [20:0] dsw_next = {dac_next[0], dac_next[1], dac_next[2], dac_next[3],
                            dac_next[4], dac_next[5],
                            msb > 4'd14, msb > 4'd13, msb > 4'd12, msb > 4'd11,
                            msb > 4'd10, msb > 4'd9,  msb > 4'd8,  msb > 4'd7,
                            msb > 4'd6,  msb > 4'd5,  msb > 4'd4,  msb > 4'd3,
                            msb > 4'd2,  msb > 4'd1,  msb > 4'd0};
    // a gated clock, not 42 enable muxes on seq[1]: GCLK pulses at the CKVd
    // edge where seq[1] is high, the edge an enable would have loaded on
    wire dclk;
`ifdef SIM
    reg dg_l;
    always @(ckvd or seq[1]) if (!ckvd) dg_l = seq[1];
    assign dclk = ckvd & dg_l;
`else
    sky130_fd_sc_hd__dlclkp_4 u_icg_dsw (.CLK(ckvd), .GATE(seq[1]), .GCLK(dclk));
`endif
    always @(posedge dclk or negedge rst_n)
        if (!rst_n) begin dsw <= 21'd0; dswb <= {21{1'b1}}; end
        else begin dsw <= dsw_next; dswb <= ~dsw_next; end

    // ---- self-timed SAR -----------------------------------------------------------
    // Each comparison: cmp_clk rises, the comparator decides (valid), the
    // decision is stored on valid's rising edge and the CDAC switches; after a
    // settling delay cmp_clk falls, the comparator resets (valid falls), and
    // after the same delay cmp_clk rises again.
    wire valid = cmp_p | cmp_n;
    wire valid_d;
`ifdef SIM
    assign #300 valid_d = valid;
`else
    wire valid_d1;
    sky130_fd_sc_hd__dlygate4sd2_1  u_dly_v1 (.A(valid), .X(valid_d1));
    sky130_fd_sc_hd__dlymetal6s4s_1 u_dly_v2 (.A(valid_d1), .X(valid_d));
`endif
    reg  [6:0] step;                     // one-hot: which comparison is next
    reg  [6:0] dec;
    wire done = (step == 7'd0);          // all seven taken
    assign cmp_clk = sar_on & ~done & ~valid_d;

    always @(posedge valid or negedge sar_on)
        if (!sar_on) begin step <= 7'b0000001; dec <= 7'd0; end
        else begin
            dec  <= dec | (step & {7{cmp_p}});
            step <= step << 1;
        end
    // comparison 0 is the sign; comparison j (1..6 follow) sets weight
    // 2^(5-j) the way its decision says, and holds it to the end
    genvar k;
    generate for (k = 0; k < 6; k = k + 1) begin : g_w
        wire taken = done | (|(step & ~((7'b0000010 << (5 - k)) - 7'd1)));
        assign bpos[k] = sar_on & taken &  dec[5 - k];
        assign bneg[k] = sar_on & taken & ~dec[5 - k];
        assign bmid[k] = ~(bpos[k] | bneg[k]);
    end endgenerate

    // the code is taken at stage 3, through a clock gate like the DAC's
    wire cclk;
`ifdef SIM
    reg cg_l;
    always @(ckvd or seq[3]) if (!ckvd) cg_l = seq[3];
    assign cclk = ckvd & cg_l;
`else
    sky130_fd_sc_hd__dlclkp_2 u_icg_code (.CLK(ckvd), .GATE(seq[3]), .GCLK(cclk));
`endif
    always @(posedge cclk or negedge rst_n)
        if (!rst_n) begin code <= 7'd0; railed <= 1'b0; end
        else begin
            // offset binary: dec[0] is the sign, dec[1..6] the magnitude bits
            code   <= {dec[0], dec[1], dec[2], dec[3], dec[4], dec[5], dec[6]};
            railed <= (dec == 7'h7F) | (dec == 7'h00);
        end
    // high for the CKVd cycle after the code is taken
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) code_valid <= 1'b0; else code_valid <= seq[3];

    // ---- the CKVd count at REF, for the frequency-lock path ----------------------
    reg [7:0] gray;
    // Gray code of the count a cycle late (one XOR level, not an increment);
    // the REF side adds the cycle back
    always @(posedge ckvd or negedge rst_n)
        if (!rst_n) gray <= 8'd0; else gray <= cnt ^ (cnt >> 1);
    reg [7:0] g1, g2;
    integer i;
    reg [7:0] bin;
    always @(posedge ref_clk or negedge rst_n)
        if (!rst_n) begin g1 <= 8'd0; g2 <= 8'd0; count_ref <= 8'd0; end
        else begin
            g1 <= gray; g2 <= g1;
            bin[7] = g2[7];
            for (i = 6; i >= 0; i = i - 1) bin[i] = bin[i + 1] ^ g2[i];
            count_ref <= bin + 8'd1;
        end
endmodule
