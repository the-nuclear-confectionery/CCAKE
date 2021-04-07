#include <stdio.h>
#include <math.h>
#include "eostables.h"
#include <iostream>
#include <fstream>

using namespace std;

double spline ( void (*r8)(double,int&,int&,double&,double&) , void (*y)(double&,double&,int,int) ,double xval ,int &pos)

{
  double x1,x2,y1,y2,yval;
  int left,right;

	
  (*r8)(xval,left, right,x1,x2 ); //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.

  (*y)(y1,y2,left,right);
  
  yval=y1+(xval-x1)*(y2-y1)/(x2-x1);
  
  pos=left;
  

  return yval;
}

double splinepos ( void (*r8)(double&,double&,int,int) , void (*y)(double&,double&,int,int) ,double xval,int pos)

{
  double x1,x2,y1,y2,yval;
  int left,right;
  left=pos;
  right=pos+1;

  (*r8)(x1,x2,left, right); //  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.

  (*y)(y1,y2,left,right);
  
  yval=y1+(xval-x1)*(y2-y1)/(x2-x1);
  
	

  return yval;
}

void y_p( double &y1, double &y2,int left,int right)
{
	y1=ETH[left].p;
	y2=ETH[right].p;

}

void y_plow( double &y1, double &y2,int left,int right)
{
	y1=ETL[left].p;
	y2=ETL[right].p;

}

void y_T( double &y1, double &y2,int left,int right)
{
	y1=ETH[left].T;
	y2=ETH[right].T;
	

}

void y_e( double &y1, double &y2,int left,int right)
{
	y1=ETH[left].e;
	y2=ETH[right].e;

}

void y_elow( double &y1, double &y2,int left,int right)
{
	y1=ETL[left].e;
	y2=ETL[right].e;

}

void y_s( double &y1, double &y2,int left,int right)
{
	y1=ETH[left].s;
	y2=ETH[right].s;

}

void y_Tlow( double &y1, double &y2,int left,int right)
{
	y1=ETL[left].T;
	y2=ETL[right].T;

}

void y_slow( double &y1, double &y2,int left,int right)
{
	y1=ETL[left].s;
	y2=ETL[right].s;

}

void y_dtds( double &y1, double &y2,int left,int right)
{
	y1=ETH[left].dtds;
	y2=ETH[right].dtds;

}

void y_dtdslow( double &y1, double &y2,int left,int right)
{
	y1=ETL[left].dtds;
	y2=ETL[right].dtds;
	

}


void r8_s (  double xval, int &left, int &right,double &x1, double&x2 )
{

  //cout << xval << " " << left << " " << right << " " << x1 << " " << x2 << endl;

  for ( int i = 10; i <nETH; i=i+10 ) 
  {
    if ( xval < ETH[i].s ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  int lnew=left,rnew=right;
  for ( int j= lnew; j <=rnew; j++ ) 
  {
    if ( xval < ETH[j].s ) 
    {
      left = j - 1;
      right = j;
      
      
      break;
    }
  }

  x1=ETH[left].s;
  x2=ETH[right].s;
}


void r8_T (  double xval, int &left, int &right,double &x1, double&x2 )
{

  for ( int i = 10; i <nETH; i=i+10 ) 
  {
    if ( xval < ETH[i].T ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  for ( int i = left; i <right; i++ ) 
  {
    if ( xval < ETH[i].T ) 
    {
      left = i - 1;
      right = i;
      break;
    }
  }

  x1=ETH[left].T;
  x2=ETH[right].T;
}

void r8_slow (  double xval, int &left, int &right,double &x1, double&x2 )
{

  for ( int i = 10; i <nETL; i=i+10 ) 
  {
    if ( xval < ETL[i].s ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  for ( int j = left; j <right; j++ ) 
  {
    if ( xval < ETL[j].s ) 
    {
      left = j - 1;
      right = j;
      break;
    }
  }

  x1=ETL[left].s;
  x2=ETL[right].s;
}

void r8_Tlow (  double xval, int &left, int &right,double &x1, double&x2 )
{

  for ( int i = 10; i <nETL; i=i+10 ) 
  {
    if ( xval < ETL[i].T ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  for ( int i = left; i <right; i++ ) 
  {
    if ( xval < ETL[i].T ) 
    {
      left = i - 1;
      right = i;
      break;
    }
  }

  x1=ETL[left].T;
  x2=ETL[right].T;
}

void r8_e (  double xval, int &left, int &right,double &x1, double&x2 )
{

  for ( int i = 10; i <nETH; i=i+10 ) 
  {
    if ( xval < ETH[i].e ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  for ( int i = left; i <right; i++ ) 
  {
    if ( xval < ETH[i].e ) 
    {
      left = i - 1;
      right = i;
      break;
    }
  }

  x1=ETH[left].e;
  x2=ETH[right].e;
}

void r8_elow (  double xval, int &left, int &right,double &x1, double&x2 )
{

  for ( int i = 10; i <nETL; i=i+10 ) 
  {
    if ( xval < ETL[i].e ) 
    {
      left = i - 10;
      right = i;
      break;
      
    }    
  }
  
  for ( int i = left; i <right; i++ ) 
  {
    if ( xval < ETL[i].e ) 
    {
      left = i - 1;
      right = i;
      break;
    }
  }

  x1=ETL[left].e;
  x2=ETL[right].e;
}

