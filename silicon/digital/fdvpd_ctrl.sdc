# fdvpd_ctrl timing constraints.
#
#   ckv       1 GHz, the DCO output: only the divide-by-two
#   ckvd      ckv / 2, generated: the edge counter, the Gray copy, the
#             synchroniser, the sequencer and the DAC and code registers
#   ref_clk   50 MHz, asynchronous to ckv
#
# Not timed here, and why:
#   * the accumulator (acc, n_floor, dac_next) is clocked by stop_async -- the
#     gated CKVd edge -- once per 20 ns reference cycle, and the gate compares
#     against it ~18 ns later; its paths are cut and covered by the RTL and
#     gate-level simulations
#   * the outputs -- see below
#   * the SAR's step / dec flops are clocked by the comparator's own "decided"
#     signal (self-timed); their loop runs through the analog comparator and is
#     checked by gate-level simulation with extracted delays
create_clock -name ckv -period 1.0 [get_ports ckv]
# ckvd is declared on the divider's output buffer, an instantiated cell with
# a fixed name.  Two earlier anchors were wrong in opposite ways: "the driver
# of net ckvd" found the buffer the resizer put on the port, so the clock
# reached the port and nothing inside (129 register pins untimed, and the run
# looked clean); u_div/Q made the divider's own feedback a clock net, and CTS
# buffered it.  STA takes the source latency through u_div's clk-to-Q and
# the inverter; -invert because the buffer hangs off the inverter.
create_generated_clock -name ckvd -source [get_ports ckv] -divide_by 2 -invert \
    [get_pins u_ckvd_buf/X]
create_clock -name ref_clk -period 20.0 [get_ports ref_clk]
set_clock_groups -asynchronous -group {ckv ckvd} -group {ref_clk}

set_input_delay 0.2 -clock ckv [get_ports {fcw_int* fcw_frac* trim_in*}]
set_input_delay 0.0 -clock ckvd [get_ports {cmp_p cmp_n}]
set_false_path -from [get_ports rst_n]
set_false_path -from [get_ports {fcw_int* fcw_frac* trim_in*}]
set_false_path -through [get_nets {acc*}]
set_false_path -through [get_nets {n_floor*}]
set_false_path -through [get_nets {dac_next*}]
set_false_path -from [get_ports {cmp_p cmp_n}]
set_false_path -through [get_nets {step*}]
set_false_path -through [get_nets {dec*}]
# timeout only ever ends the window, as the gated edge does, so a glitch where
# the two overlap is another stop: no clock-gating check on it
set_false_path -through [get_pins u_stop/B]

# The outputs drive analog switches and the loop filter, not ckv flops: their
# arrival is checked against the analog timing in the mixed-signal simulation
set_false_path -to [all_outputs]

set_load 0.02 [all_outputs]
set_driving_cell -lib_cell sky130_fd_sc_hd__buf_2 -pin X [all_inputs]
