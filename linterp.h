#ifndef _linterp_h
#define _linterp_h

#define CONFIG_DELIMS ",/|\\\n\t !@#$%^&*()-_=+`~"

/* relates freq:amplitude */
typedef struct {
	double freq;
	double ampl;
} freq_ampl_map_t;

/* frequency amplitude map list head */
typedef struct {
	size_t n;
	freq_ampl_map_t *maps;
} maplist_t;

extern maplist_t parse_freq_ampl_map(FILE *);
extern double linterp(maplist_t, double);

#endif /* !defined _linterp_h */
