#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "linterp.h"

/* callback for qsort() use in parse_file() */
int _cmp_fa_maps(const void *a, const void *b)
{
	const freq_ampl_map_t *map1 = a,
		              *map2 = b;
	if (map1->freq < map2->freq) return -1;
	if (map1->freq > map2->freq) return  1;

	/* error: duplicate freqs */
	fprintf(stderr, "FIXME: Config file has duplicate points, "
			"(%lf,%lf) and (%lf, %lf)\n",
			map1->freq, map1->ampl,
			map2->freq, map2->ampl);
	return 0;
}

maplist_t parse_freq_ampl_map(FILE *f)
{
	maplist_t ret = {0, NULL};
	size_t len;
	char *buf, *tok;
	int i;

	if (!f) goto buferr;

	/* calculate file length and reset read head */
	fseek(f, 0, SEEK_END);
	len = ftell(f);
	rewind(f);

	/* read the file into memory */
	buf = malloc(len + 1);
	if (!buf) goto buferr;
	fread(buf, len, 1, f);
	buf[len] = 0;

	/* set up the list head */
	ret.n = 8;
	ret.maps = malloc(ret.n * sizeof *ret.maps);
	if (!ret.maps) goto maperr;
	
	/* finally begin parsing */
	tok = strtok(buf, CONFIG_DELIMS);
	if (!tok) goto tokerr;	
	
	for (i = 0; tok; i++, tok = strtok(NULL, CONFIG_DELIMS)) {
	
		/* grow the list if we need to */
		if (i >= ret.n) {
			ret.n *= 2;
			ret.maps = realloc(ret.maps, ret.n * sizeof *ret.maps);
			if (!ret.maps) goto maperr;
		}

		if (sscanf(tok, "%lf", &ret.maps[i].freq) <= 0)
			goto tokerr;
		
		tok = strtok(NULL, CONFIG_DELIMS);	
		/* must have an ampl for each freq */
		if (!tok) goto tokerr;

		if (sscanf(tok, "%lf", &ret.maps[i].ampl) <= 0)
			goto tokerr;

	}

	/* downsize */
	ret.n = i;
	ret.maps = realloc(ret.maps, ret.n * sizeof *ret.maps);
	if (!ret.maps) goto maperr;

	/* sort ascending */
	qsort(ret.maps, ret.n, sizeof *ret.maps, _cmp_fa_maps);
	
	/* no errors -> safe to hit normal exit path */
	goto exitok;

tokerr:	free(ret.maps);
maperr:	ret = (maplist_t) {0, NULL};
exitok: free(buf);
buferr:	return ret; /* = {0, NULL} from initialization */
}

double linterp(maplist_t head, double freq)
{
	/* head.maps is now guaranteed to be sorted *
	 * ascending by frequency, with head.n-1 as *
	 * the highest valid index into head.maps.  *
	 * This makes this function trivial.        */

	/* 1 data point -> constant amplitude */
	if (head.n == 1) return head.maps[0].ampl;

	/* freq out of bounds */
	if (freq < head.maps[0].freq)
		return head.maps[0].ampl;
	
	if (freq > head.maps[head.n-1].freq)
		return head.maps[head.n-1].ampl;

	/* in bounds - interpolate */
	for (int i = 0; i < head.n-1; i++)
		if (head.maps[i].freq < freq
		    && freq < head.maps[i+1].freq) {
			
			double m = (head.maps[i+1].ampl - head.maps[i].ampl)
				 / (head.maps[i+1].freq - head.maps[i].freq);
			double x = freq - head.maps[i].freq;
			double c = head.maps[i].ampl;

			/* and now, the equation of a line */
			return m*x + c;

		}

	/* eh */
	return NAN;
}
