This a slightly-modified version of a gain-phase analyzer
that ships with the Red Pitaya source code. It's been modified
to provide more useful output (i.e. in CSV format). The user
can also specify amplitude as a function of frequency by passing
a file containining frequency-amplitude pairs.

An example configuration file could contain the following:

	1e2:0.9
	1e3:0.4
	1e4:0.1
	1e5:0.6
	1e6:0.3

Most delimiters are accepted. Amplitudes are calculated by linear
interpolation between these frequency-amplitude pairs.

A complete list of options:
    [c|ch|chan|channel]                  set output channel {1, 2}
    [a|amp|ampl|amplitude]               set output amplitude [0-1]
    [dc|bias|dc_bias|offset]             set output dc bias [0-1]
                                         Note: bias + amplitude < 1.
    [r|range]                            set low and high freq bounds
                                         e.g. ./bode range 1e3 1e5
    [startf|startfreq|fstart|f_lo|f_low] Set low freq bound
    [stopf|stopfreq|fstop|f_hi|f_high]   Set high freq bound
    [lin|linear]                         Linear freq sweep
    [log|logarithmic]                    Logarithmic frequency sweep
    [avg|avgn|average]                   Number of samples to average together
    	     			         at each point
    [n|count|steps]                      Number of points to measure
    [fin|config]                         Use a config file for a frequency:amplitude
    				         curve
    [fout|out|csvout|csv]                Put the output CSV here

Example invocations:
./bode config foo.conf n 50 log csv out.csv avg 3 range 1e2 1e5
./bode count 30 lin csv foo.csv range 10 25e4 bias 0.2 a 0.3 chan 1
