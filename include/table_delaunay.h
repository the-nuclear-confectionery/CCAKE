#ifndef TABLE_DELAUNAY_H
#define TABLE_DELAUNAY_H

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>

using namespace std;

int generate_delaunay_2d( const vector<vector<double> > & vertices,
                          vector<vector<double> > & faces );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2, 
             double x3, double y3 );
int dtris2 ( int point_num, int base, double point_xy[], int *tri_num, 
             int tri_vert[], int tri_nabe[] );
int file_column_count ( string input_file_name );
int file_row_count ( string input_file_name );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_sign ( int i );
int i4_wrap ( int ival, int ilo, int ihi );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
                                  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
int *i4vec_indicator_new ( int n );
int i4vec_min ( int n, int a[] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2, 
             double yv2, double dv );
bool perm_check ( int n, int p[], int base );
void perm_inverse ( int n, int p[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r82vec_permute ( int n, int p[], int base, double a[] );
int *r82vec_sort_heap_index_a ( int n, int base, double a[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
                                  int ihi, int jhi, string title );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
int swapec ( int i, int *top, int *btri, int *bedg, int point_num, 
             double point_xy[], int tri_num, int tri_vert[], int tri_nabe[], 
             int stack[] );
void timestamp ( );
void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num, 
             int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );

#endif
