#ifndef MINDISTANCE_H
#define MINDISTANCE_H
/*
int standard_min_gp(double,
		    double,
		    double, 
		    int, struct grid*);
*/
int standard_min(double, double*,
		 double, double*,
		 double, double*);
int index_min(double, double*, int);
int find_id_min(double, double*,
		double, double*,
		double, double*);
#endif
