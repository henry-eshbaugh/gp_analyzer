/**
 * $Id: bode.c 1246  $
 *
 * @brief Red Pitaya Bode plotter
 *
 * @Author1 Martin Cimerman (main developer,concept program translation)
 * @Author2 Zumret Topcagic (concept code developer)
 * @Author3 Luka Golinar (functioanlity of web interface)
 * @Author4 Peter Miklavcic (manpage and code review)
 * Contact: <cim.martin@gmail.com>, <luka.golinar@gmail.com>
 *
 * GENERAL DESCRIPTION:
 *
 * The code below defines the Bode analyzer on a Red Pitaya.
 * It uses acquire and generate from the Test/ folder.
 * Data analysis returns frequency, phase and amplitude.
 *
 * VERSION: VERSION defined in Makefile
 *
 * This part of code is written in C programming language.
 * Please visit http://en.wikipedia.org/wiki/C_(programming_language)
 * for more details on the language used herein.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/param.h>

#include "argparse.h"
#include "linterp.h"
#include "main_osc.h"
#include "fpga_osc.h"
#include "fpga_awg.h"
#include "version.h"

#define M_PI 3.14159265358979323846

const char *g_argv0 = NULL; 		/* Program name 		 */
const double c_max_frequency = 62.5e6;  /* Maximal signal frequency [Hz] */
const double c_min_frequency = 0; 	/* Minimal signal frequency [Hz] */
const double c_max_amplitude = 1.0; 	/* Maximal signal amplitude [V]  */

#define AWG_BUF_LEN (16*1024) // AWG buffer length [samples]
int32_t data[AWG_BUF_LEN]; // AWG data buffer

/** Signal types */
typedef enum {
	eSignalSine,     /* Sinusoidal waveform */
	eSignalSquare,   /* Square waveform */
	eSignalTriangle, /* Triangular waveform */
	eSignalSweep,	 /* Sinusoidal frequency sweep */
	eSignalConst     /* Constant signal */
} signal_e;

/** AWG FPGA parameters */
typedef struct {
	int32_t  offsgain;   // AWG offset & gain
	uint32_t wrap;	   // AWG buffer wrap value
	uint32_t step;	   // AWG step interval
} awg_param_t;


/** Oscilloscope module parameters as defined in main module
 * @see rp_main_params
 */
float t_params[PARAMS_NUM] = { 0, 1e6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/** Decimation translation table */
#define DEC_MAX 6 // Max decimation index
static int g_dec[DEC_MAX] = { 1,  8,  64,  1024,  8192,  65536 };

/** Forward declarations */
void write_data_fpga(uint32_t ch, const int32_t *data, const awg_param_t *awg);
int acquire_data(float **s, uint32_t size);
void synthesize_signal(double ampl, double freq, signal_e type, double endfreq,
					   int32_t *data, awg_param_t *params);
int bode_data_analysis(float **s, uint32_t size,
					   double DC_bias, float *Amplitude,
					   float *Phase, double w_out, int f);

/** Print usage information */
#define usage() _usage(__FILE__, __LINE__, __func__)
void _usage(char *file, int line, const char *func) {
	const char *format =
		"%s: Bode analyzer, build %s\n"
		"Usage:\n"
		"\t[c|ch|chan|channel]                  set output channel {1, 2}\n"
		"\t[a|amp|ampl|amplitude]               set output amplitude [0-1]\n"
		"\t[dc|bias|dc_bias|offset]             set output dc bias [0-1]\n"
		"\t                                     Note: bias + amplitude < 1.\n"
		"\t[r|range]                            set low and high freq bounds\n"
		"\t                                     e.g. ./bode range 1e3 1e5\n"
		"\t[startf|startfreq|fstart|f_lo|f_low] Set low freq bound\n"
		"\t[stopf|stopfreq|fstop|f_hi|f_high]   Set high freq bound\n"
		"\t[lin|linear]                         Linear freq sweep\n"
		"\t[log|logarithmic]                    Logarithmic frequency sweep\n"
		"\t[avg|avgn|average]                   Number of samples to average together at each point\n"
		"\t[n|count|steps]                      Number of points to measure\n"
		"\t[fin|config]                         Use a config file for a frequency:amplitude curve\n"
		"\t[fout|out|csvout|csv]                Put the output CSV here\n"
		"----------------------------------------------------------------\n"
		"Example invocations:\n"
		"\t%s config foo.conf n 50 log csv out.csv avg 3 range 1e2 1e5\n"
		"\t%s count 30 lin csv foo.csv range 10 25e4 bias 0.2 a 0.3 chan 1\n"
		"----------------------------------------------------------------\n"
		"If you're seeing this output, it's likely something went wrong.\n"
		"Here's where we were called: %s:(%d) from %s()\n";
	fprintf(stderr, format, g_argv0, __TIMESTAMP__, g_argv0, g_argv0, file,
		line, func);
}

/** Allocates a memory with size num_of_el, memory has 1 dimension */
float *create_table_size(int n)
{
	float *new_table = malloc(n * sizeof *new_table);
	return new_table;
}

/** Allocates a memory with size num_of_cols*num_of_rows */
float **create_2D_table_size(int num_of_rows, int num_of_cols)
{
	float **new_table = malloc(num_of_rows * sizeof *new_table);
	for(int i = 0; i < num_of_rows; i++)
		new_table[i] = create_table_size(num_of_cols);
	return new_table;
}

float max_array(float *arrayptr, int numofelements)
{
	float max = arrayptr[0];
	for (int i = 0; i < numofelements; i++)
		max = max > arrayptr[i] ? max
					: arrayptr[i];
	return max;
}

/** Trapezoidal method for integration */
float trapz(float *arrayptr, float T, int size1)
{
	float result = 0;
	int i;

	for (i = 0; i < size1-1; i++)
		result += arrayptr[i] + arrayptr[i+1];

	result *= T * 0.5f;
	return result;
}

/** Finds a mean value of an array */
float mean_array(float *arrayptr, int numofelements)
{
	float mean = 0;

	for(int i = 0; i < numofelements; i++)
		mean += arrayptr[i];

	mean = mean / numofelements;
	return mean;
}

/** Finds a mean value of an array by columns, acquiring values from rows */
float mean_array_column(float **arrayptr, int length, int column)
{
	float result = 0.0f;

	for(int i = 0; i < length; i++)
		result += arrayptr[i][column];

	result /= length;
	return result;
}

/* TODAY: TODO:
 * redo argument parsing here;
 * include linterp.h stuff.
 */

/* EXISTING ARGUMENTS
 * Channel: {1, 2} : channel to output on
 * amplitude: [0, 1] : amplitude of o/p signal
 * dc bias: [0, 1] : output dc offset
 * 	note dc bias + amplitude < 1
 * averaging
 * number of steps
 */

/* NEW ARGUMENTS
 * File for amplitude/freq mappings
 */

/** Bode analyzer */
int main(int argc, char *argv[])
{
	/** Set program name */
	g_argv0 = argv[0];

	unsigned int ch = 0, scale_type = 1,
		     averaging_num = 2, steps = 10;
	double ampl = 0.5, DC_bias = 0.0;
	double start_frequency, end_frequency;
	FILE *fin = NULL, *fout = NULL;

	LOOP_ARGV {
		IF_ARG("c")              goto channel;
		else IF_ARG("ch")        goto channel;
		else IF_ARG("chan")      goto channel;
		else IF_ARG("channel") {
			channel:
			NEXT_ARG_MUST_EXIST;
			ch = atoi(THIS_ARG) - 1;
		}

		else IF_ARG("a")         goto amplitude;
		else IF_ARG("amp")       goto amplitude;
		else IF_ARG("ampl")      goto amplitude;
		else IF_ARG("amplitude") {
			amplitude:
			NEXT_ARG_MUST_EXIST;
			ampl = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("dc")        goto dc_bias;
		else IF_ARG("bias")      goto dc_bias;
		else IF_ARG("dc_bias")   goto dc_bias;
		else IF_ARG("offset") {
			dc_bias:
			NEXT_ARG_MUST_EXIST;
			DC_bias = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("r")         goto range;
		else IF_ARG("range") {
			range:
			NEXT_ARG_MUST_EXIST;
			start_frequency = strtod(THIS_ARG, NULL);
			NEXT_ARG_MUST_EXIST;
			end_frequency = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("startf")    goto start_freq;
		else IF_ARG("startfreq") goto start_freq;
		else IF_ARG("fstart")    goto start_freq;
		else IF_ARG("f_lo")      goto start_freq;
		else IF_ARG("f_low") {
			start_freq:
			NEXT_ARG_MUST_EXIST;
			start_frequency = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("stopf")     goto stop_freq;
		else IF_ARG("stopfreq")  goto stop_freq;
		else IF_ARG("fstop")     goto stop_freq;
		else IF_ARG("f_hi")      goto stop_freq;
		else IF_ARG("f_high") {
			stop_freq:
			NEXT_ARG_MUST_EXIST;
			end_frequency = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("lin")       goto lin;
		else IF_ARG("linear")
			lin: scale_type = 0;

		else IF_ARG("log")       goto log;
		else IF_ARG("logarithmic")
			log: scale_type = 1;

		else IF_ARG("avg")       goto average;
		else IF_ARG("avgn")      goto average;
		else IF_ARG("average") {
			average:
			NEXT_ARG_MUST_EXIST;
			averaging_num = strtod(THIS_ARG, NULL);
		}

		else IF_ARG("n")         goto steps;
		else IF_ARG("count")     goto steps;
		else IF_ARG("steps") {
			steps:
			NEXT_ARG_MUST_EXIST;
			steps = atoi(THIS_ARG);
		}

		else IF_ARG("fin") goto fin;
		else IF_ARG("config") {
			fin:
			NEXT_ARG_MUST_EXIST;
			fin = fopen(THIS_ARG, "r");
			if (!fin) printf("Warning: no file\n");
		}

		else IF_ARG("fout") goto fout;
		else IF_ARG("out")  goto fout;
		else IF_ARG("csvout") goto fout;
		else IF_ARG("csv") {
			fout:
			NEXT_ARG_MUST_EXIST;
			fout = fopen(THIS_ARG, "w");
			if (!fout) printf("Warning: no file\n");
		}

		else usage();

	}

	maplist_t map;
	if (fin) map = parse_freq_ampl_map(fin);	
	if (!fout) fout = stdout;

	/* basically these are runtime asserts */
	CHECK_OK(ch > 1,                            "Invalid channel value.\n");
	CHECK_OK(ampl < 0,                          "Invalid amplitude.\n");
	CHECK_OK(ampl > c_max_amplitude,            "Amplitude too high\n");
	CHECK_OK(DC_bias < 0,                       "Got negative dc bias.\n");
	CHECK_OK(DC_bias > 1,                       "dc bias too high\n");
	CHECK_OK(ampl + DC_bias > 1,                "Invalid ampltiude/bias combination\n");
	CHECK_OK(ampl + DC_bias <= 0,               "Invalid amplitude/bias combination.\n");
	CHECK_OK(averaging_num < 1,                 "Invalid averaging value.\n");
	CHECK_OK(steps < 2,                         "Invalid number of steps.\n");
	CHECK_OK(start_frequency < c_min_frequency, "Start frequency is too low\n");
	CHECK_OK(start_frequency > c_max_frequency, "Start frequency is too high\n");
	CHECK_OK(end_frequency < c_min_frequency,   "Stop frequency is too low\n");
	CHECK_OK(end_frequency > c_max_frequency,   "Stop frequency is too high\n");
	CHECK_OK(end_frequency < start_frequency,   "Stop frequency is less than the start frequency\n");

	/** Parameters initialization and calculation */
	double          frequency_step;
	double          a, b, c;
	double		endfreq = 0;		/* endfreq set for generate's sweep */
	double		k;
	double		w_out;  		/* angular velocity */
	uint32_t        min_periods = 10;       /* max 20 */
	uint32_t        size; 			/* number of samples varies with number of periods */
	signal_e        type = eSignalSine;
	int             f = 0; 			/* decimation index */
	int             i1, fr; 		/* iterators in for loops */
	int             equal = 0; 		/* parameter initialized for generator functionality */
	int             shaping = 0; 		/* parameter initialized for generator functionality */
	int             transientEffectFlag = 1;
	int             stepsTE = 10; 		/* number of steps for transient effect(TE) elimination */
	int             TE_step_counter;
	int             progress_int;
	char            command[70];
	char            hex[45];

	/* if user sets less than 10 steps than stepsTE is decreased	      *
	 * for transient efect to be eliminated only 10 steps of measurements *
	 * is eliminated						      */
	if (steps < 10)
		stepsTE = steps;

	TE_step_counter = stepsTE;

	/* handle log scale */
	if (scale_type)
		b = log10f(end_frequency),
		a = log10f(start_frequency),
		c = (b - a) / (steps - 1);
		/* ^ this is one expr, note the commas */
	else
		frequency_step = (end_frequency - start_frequency) / (steps - 1);


	float **s                      = create_2D_table_size(SIGNALS_NUM, SIGNAL_LENGTH); // raw data saved to this location
	float *Amplitude               = malloc((averaging_num + 1) * sizeof(float));
	float *Amplitude_output        = malloc((steps + 1) * sizeof(float));
	float *Phase                   = malloc((averaging_num + 1) * sizeof(float));
	float *Phase_output            = malloc((steps + 1) * sizeof(float));
	float **data_for_averaging     = create_2D_table_size((averaging_num + 1), 2);
	float *measured_data_amplitude = malloc((2) * sizeof(float));
	float *measured_data_phase     = malloc((2) * sizeof(float));
	float *frequency               = malloc((steps + 1) * sizeof(float));

	/* Initialization of Oscilloscope application */
	if (rp_app_init() < 0) {
		fprintf(stderr, "rp_app_init() failed!\n");
		return -1;
	}

	/* ----- generate/acquire loop ----- */
	for (fr = 0; fr < steps; fr++) {
		/* set scale type */
		if (scale_type) {
			k = powf(10, (c * (float)fr) + a);
			frequency[fr] = k;
		} else {
			frequency[fr] = start_frequency + (frequency_step * fr);
		}

		/* account for transient effects by outputting low freqs */
		if (transientEffectFlag == 1){

			if (TE_step_counter > 0) {
				int k = TE_step_counter/stepsTE - 1;
				frequency[fr] = start_frequency;
				frequency[fr] += k * start_frequency / 2;
				TE_step_counter--;
			}

			if (TE_step_counter == 0){
				fr = 0;
				frequency[fr] = start_frequency;
				transientEffectFlag = 0;
			}

		}

		/* progress bar has transient effect accounting      */
		/* TODO: get ncurses on RP, have a real progress bar */
		progress_int = transientEffectFlag == 1 ?
				  stepsTE - TE_step_counter
				: stepsTE + fr;
		progress_int *= 100;
		progress_int /= steps + stepsTE - 1;

		if (progress_int <= 100) {
			sprintf(hex, "%x", (int)(255 - (255*progress_int/100)));
			strcpy(command, "/opt/redpitaya/bin/monitor 0x40000030 0x" );
			strcat(command, hex);
			system(command);
		}

		/* ----- signal synthesis ----- */
		w_out = frequency[fr] * 2 * M_PI; // w = omega = angular freq
		awg_param_t params;
		/* if given a config file, use it to calculate amplitude */
		if (fin)
			ampl = linterp(map, frequency[fr]);
		/* prepare the signal buffer */
		synthesize_signal(ampl, frequency[fr], type, endfreq, data, &params);
		/* Write the data to the FPGA and set FPGA AWG state machine */
		write_data_fpga(ch, data, &params);
		/* sleep 1ms for synchronization*/
		usleep(1000);

		/* ----- signal acquisition ----- */
		for (i1 = 0; i1 < averaging_num; i1++) {
			/* find appropriate decimation index */
			     if	(frequency[fr] >= 160000) f = 0;
			else if (frequency[fr] >= 20000)  f = 1;
			else if (frequency[fr] >= 2500)   f = 2;
			else if (frequency[fr] >= 160)    f = 3;
			else if (frequency[fr] >= 20)	  f = 4;
			else if (frequency[fr] >= 2.5)    f = 5;

			/* set decimtion index */
			if (f != DEC_MAX) {
				t_params[TIME_RANGE_PARAM] = f;
			} else {
				fprintf(stderr, "Invalid decimation %d\n", f);
				usage();
				return -1;
			}

			/* calculate number of samples */
			size = round(  (min_periods*125e6)
				     / (frequency[fr] * g_dec[f]));

			/* set scope module parameters for signal acqusition */
			t_params[EQUAL_FILT_PARAM] = equal;
			t_params[SHAPE_FILT_PARAM] = shaping;
			if(rp_set_params((float *)&t_params, PARAMS_NUM) < 0) {
				fprintf(stderr, "rp_set_params() failed!\n");
				return -1;
			}

			/* acquire the signal */
			if (acquire_data(s, size) < 0) {
				printf("error acquiring data @ acquire_data\n");
				return -1;
			}

			/* analysis */
			if (bode_data_analysis(s, size, DC_bias, Amplitude,
					       Phase, w_out, f) < 0) {
				printf("error in bode_data_analysis()\n");
				return -1;
			}

			/* Saving data */
			data_for_averaging[i1][1] = *Amplitude;
			data_for_averaging[i1][2] = *Phase;
		} /* averaging loop end */

		/* Calculating and saving mean values */
		measured_data_amplitude[1] = mean_array_column(data_for_averaging, averaging_num, 1);
		measured_data_phase[1]	   = mean_array_column(data_for_averaging, averaging_num, 2 );

		if (transientEffectFlag == 0) {
			Amplitude_output[fr] = measured_data_amplitude[1];
			Phase_output[fr] = measured_data_phase[1];
		}
	}
	
	/* ----- done generating/acquiring, turn off the output ----- */
	awg_param_t params;
	synthesize_signal(0, 1000, type, endfreq, data, &params);
	write_data_fpga(ch, data, &params);


	fprintf(fout, "Freq,Phi,A");
	for (int po = 0; po < steps; po++)
		fprintf(fout, "%.2f,%.5f,%.5f\n", frequency[po], Phase_output[po], Amplitude_output[po]);

	return 1;
}

/* A lot of this code touches the FPGA directly. That's okay. It doesn't
 * use the normal rp_*() family of functions. That's also okay. This works,
 * and there's enough uncertainty in the way the RP folks do things that
 * we don't really want to question *why* something works. This is okay.
 */

/**
 * Synthesize a desired signal.
 *
 * Generates/synthesized  a signal, based on three pre-defined signal
 * types/shapes, signal amplitude & frequency. The data[] vector of
 * samples at 125 MHz is generated to be re-played by the FPGA AWG module.
 *
 * @param ampl  Signal amplitude [V].
 * @param freq  Signal Frequency [Hz].
 * @param type  Signal type/shape [Sine, Square, Triangle, Constant].
 * @param data  Returned synthesized AWG data vector.
 * @param awg   Returned AWG parameters.
 *
 */
void synthesize_signal(double ampl, double freq,
		       signal_e type, double endfreq,
		       int32_t *data, awg_param_t *awg)
{

	uint32_t i;

	/* Various locally used constants - HW specific parameters */
	const int dcoffs = -155;
	const int trans0 = 30;
	const int trans1 = 300;
	const double tt2 = 0.249;

	/* This is where frequency is used... */
	awg->offsgain = (dcoffs << 16) + 0x1fff;
	awg->step = round(65536 * freq/c_awg_smpl_freq * AWG_BUF_LEN);
	awg->wrap = round(65536 * AWG_BUF_LEN - 1);

	int trans = freq / 1e6 * trans1; /* 300 samples at 1 MHz */
	uint32_t amp = ampl * 4000.0;	/* 1 V ==> 4000 DAC counts */
	if (amp > 8191) {
		/* Truncate to max value if needed */
		amp = 8191;
	}

	if (trans <= 10) {
		trans = trans0;
	}


	/* Fill data[] with appropriate buffer samples */
	for(i = 0; i < AWG_BUF_LEN; i++) {

		/* Sine */
		if (type == eSignalSine) {
			data[i] = round(amp * cos(2*M_PI*(double)i/(double)AWG_BUF_LEN));
		}

		/* Square */
		if (type == eSignalSquare) {
			data[i] = round(amp * cos(2*M_PI*(double)i/(double)AWG_BUF_LEN));
			if (data[i] > 0)
				data[i] = amp;
			else
				data[i] = -amp;

			/* Soft linear transitions */
			double mm, qq, xx, xm;
			double x1, x2, y1, y2;

			xx = i;
			xm = AWG_BUF_LEN;
			mm = -2.0*(double)amp/(double)trans;
			qq = (double)amp * (2 + xm/(2.0*(double)trans));

			x1 = xm * tt2;
			x2 = xm * tt2 + (double)trans;

			if ( (xx > x1) && (xx <= x2) ) {

				y1 = (double)amp;
				y2 = -(double)amp;

				mm = (y2 - y1) / (x2 - x1);
				qq = y1 - mm * x1;

				data[i] = round(mm * xx + qq);
			}

			x1 = xm * 0.75;
			x2 = xm * 0.75 + trans;

			if ( (xx > x1) && (xx <= x2)) {

				y1 = -(double)amp;
				y2 = (double)amp;

				mm = (y2 - y1) / (x2 - x1);
				qq = y1 - mm * x1;

				data[i] = round(mm * xx + qq);
			}
		}

		/* Triangle */
		if (type == eSignalTriangle) {
			data[i] = round(-1.0*(double)amp*(acos(cos(2*M_PI*(double)i/(double)AWG_BUF_LEN))/M_PI*2-1));
		}

		/* Sweep */
		/* Loops from i = 0 to n = 16*1024. Generates a sine wave signal that
		   changes in frequency as the buffer is filled. */
		double start = 2 * M_PI * freq;
		double end = 2 * M_PI * endfreq;
		if (type == eSignalSweep) {
			double sampFreq = c_awg_smpl_freq; // 125 MHz
			double t = i / sampFreq; // This particular sample
			double T = AWG_BUF_LEN / sampFreq; // Wave period = # samples / sample frequency
			/* Actual formula. Frequency changes from start to end. */
			data[i] = round(amp * (sin((start*T)/log(end/start) * ((exp(t*log(end/start)/T)-1)))));
		}

		/* Constant */
		if (type == eSignalConst) data[i] = amp;

		/* TODO: Remove, not necessary in C/C++. */
		if(data[i] < 0)
			data[i] += (1 << 14);
	}
}

/**
 * Write synthesized data[] to FPGA buffer.
 *
 * @param ch	Channel number [0, 1].
 * @param data  AWG data to write to FPGA.
 * @param awg   AWG paramters to write to FPGA.
 */
void write_data_fpga(uint32_t ch,
					 const int32_t *data,
					 const awg_param_t *awg) {

	uint32_t i;

	fpga_awg_init();

	if(ch == 0) {
		/* Channel A */
		g_awg_reg->state_machine_conf = 0x000041;
		g_awg_reg->cha_scale_off	  = awg->offsgain;
		g_awg_reg->cha_count_wrap	 = awg->wrap;
		g_awg_reg->cha_count_step	 = awg->step;
		g_awg_reg->cha_start_off	  = 0;

		for(i = 0; i < AWG_BUF_LEN; i++) {
			g_awg_cha_mem[i] = data[i];
		}
	} else {
		/* Channel B */
		g_awg_reg->state_machine_conf = 0x410000;
		g_awg_reg->chb_scale_off	  = awg->offsgain;
		g_awg_reg->chb_count_wrap	 = awg->wrap;
		g_awg_reg->chb_count_step	 = awg->step;
		g_awg_reg->chb_start_off	  = 0;

		for(i = 0; i < AWG_BUF_LEN; i++) {
			g_awg_chb_mem[i] = data[i];
		}
	}

	/* Enable both channels */
	/* TODO: Should this only happen for the specified channel?
	 *	   Otherwise, the not-to-be-affected channel is restarted as well
	 *	   causing unwanted disturbances on that channel.
	 */
	g_awg_reg->state_machine_conf = 0x110011;

	fpga_awg_exit();
}

/**
 * Acquire data from FPGA to memory (s).
 *
 * @param **s   Points to a memory where data is saved.
 * @param size  Size of data.
 */
int acquire_data(float **s ,uint32_t size) {

	int retries = 150000;
	int j, sig_num, sig_len;
	int ret_val;
	usleep(50000);
	while(retries >= 0) {
		if((ret_val = rp_get_signals(&s, &sig_num, &sig_len)) >= 0) {
			/* Signals acquired in s[][]:
			 * s[0][i] - TODO
			 * s[1][i] - Channel ADC1 raw signal
			 * s[2][i] - Channel ADC2 raw signal
			 */
			for(j = 0; j < MIN(size, sig_len); j++) {
				//printf("%7d, %7d\n",(int)s[1][j], (int)s[2][j]);
			}
			break;
		}
		if(retries-- == 0) {
			fprintf(stderr, "Signal scquisition was not triggered!\n");
			break;
		}
		usleep(1000);
	}
	usleep(30000); // delay for pitaya to operate correctly
	return 1;
}

/**
 * Acquired data analysis function for Bode analyzer.
 *
 * @param s		  Pointer where data is read from.
 * @param size	   Size of data.
 * @param DC_bias	DC component.
 * @param Amplitude  Pointer where to write amplitude data.
 * @param Phase	  Pointer where to write phase data.
 * @param w_out	  Angular velocity (2*pi*freq).
 * @param f		  Decimation selector index.
 */
int
bode_data_analysis(float **s,
		   uint32_t size,
		   double DC_bias,
		   float *Amplitude,
		   float *Phase,
		   double w_out,
		   int f)
{
	int i2, i3;
	float **U_acq = create_2D_table_size(SIGNALS_NUM, SIGNAL_LENGTH);
	/* Signals multiplied by the reference signal (sin) */
	float *U1_sampled_X = (float *) malloc( size * sizeof( float ) );
	float *U1_sampled_Y = (float *) malloc( size * sizeof( float ) );
	float *U2_sampled_X = (float *) malloc( size * sizeof( float ) );
	float *U2_sampled_Y = (float *) malloc( size * sizeof( float ) );
	/* Signals return by trapezoidal method in complex */
	float *X_component_lock_in_1 = (float *) malloc( size * sizeof( float ) );
	float *X_component_lock_in_2 = (float *) malloc( size * sizeof( float ) );
	float *Y_component_lock_in_1 = (float *) malloc( size * sizeof( float ) );
	float *Y_component_lock_in_2 = (float *) malloc( size * sizeof( float ) );
	/* Voltage, current and their phases calculated */
	float U1_amp;
	float Phase_U1_amp;
	float U2_amp;
	float Phase_U2_amp;
	float Phase_internal;
	//float Z_phase_deg_imag;  // may cuse errors because not complex
	float T; // Sampling time in seconds
	float *t = create_table_size(16384);

	T = ( g_dec[f] / 125e6 );

	for(i2 = 0; i2 < (size - 1); i2++)
		t[i2] = i2;

	/* Transform signals from  AD - 14 bit to voltage [ ( s / 2^14 ) * 2 ] */
	for (i2 = 0; i2 < SIGNALS_NUM; i2++) // only the 1 and 2 are used for i2
		for(i3 = 0; i3 < size; i3++)
			/* sequencing is important here;     */
			/* multiply before dividing to avoid */
			/* losing sigfigs */
			U_acq[i2][i3] =
				((s[i2][i3]) * (float)(2 - DC_bias))
				/ 16384.0f;

	/* Acquired signals must be multiplied by the reference signals, used for lock in metod */
	float ang;
	for( i2 = 0; i2 < size; i2++) {
		ang = (i2 * T * w_out);
		//printf("ang(%d) = %f \n", (i2+1), ang);
		U1_sampled_X[i2] = U_acq[1][i2] * sin( ang );
		U1_sampled_Y[i2] = U_acq[1][i2] * sin( ang+ (M_PI/2) );

		U2_sampled_X[i2] = U_acq[2][i2] * sin( ang );
		U2_sampled_Y[i2] = U_acq[2][i2] * sin( ang +(M_PI/2) );
	}

	/* Trapezoidal method for calculating the approximation of an integral */
	X_component_lock_in_1[1] = trapz( U1_sampled_X, (float)T, size );
	Y_component_lock_in_1[1] = trapz( U1_sampled_Y, (float)T, size );

	X_component_lock_in_2[1] = trapz( U2_sampled_X, (float)T, size );
	Y_component_lock_in_2[1] = trapz( U2_sampled_Y, (float)T, size );

	/* Calculating voltage amplitude and phase */
	U1_amp = (float)2 * (sqrtf( powf( X_component_lock_in_1[ 1 ] , (float)2 ) + powf( Y_component_lock_in_1[ 1 ] , (float)2 )));
	Phase_U1_amp = atan2f( Y_component_lock_in_1[ 1 ], X_component_lock_in_1[ 1 ] );

	/* Calculating current amplitude and phase */
	U2_amp = (float)2 * (sqrtf( powf( X_component_lock_in_2[ 1 ], (float)2 ) + powf( Y_component_lock_in_2[ 1 ] , (float)2 ) ) );
	Phase_U2_amp = atan2f( Y_component_lock_in_2[1], X_component_lock_in_2[1] );

	Phase_internal = Phase_U2_amp - Phase_U1_amp;

	if (Phase_internal <=  (-M_PI) )
	{
		Phase_internal = Phase_internal +(2*M_PI);
	}
	else if ( Phase_internal >= M_PI )
	{
		Phase_internal = Phase_internal -(2*M_PI);
	}
	else
	{
		Phase_internal = Phase_internal;
	}

	*Amplitude = 10*log( U2_amp / U1_amp );
	*Phase = Phase_internal * ( 180/M_PI );

	return 1;
}
