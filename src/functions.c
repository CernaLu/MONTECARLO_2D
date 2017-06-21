#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "declare.h"

void data(float *L, int *N, float *n,
		float *bvol, float *sigma, float *delta) {
        
  	//#### BOX SIZE #####
	float l;
  	printf("Box size: ");
	  scanf("%f",&l);
	*L=l;
	*bvol = l*l;
	//###################	
	
		//########## Nº OF PARTICLES ######
		float nn, nnn, mod;
		mod = 1.0;
		while( mod != 0.0 ){
		printf("Set number of particles: ");
		  scanf("%f",&nn);
		nnn = sqrt(nn);
		mod = fmod(nnn,1.0);
		}
		*N = (int)nn;
		*n = nnn;
		//################################

			//########### PARTICLE DIM #####
			float sig;
			printf("Particle diameter: ");
			  scanf("%f",&sig);
			*sigma=sig;
			//##############################
	
				//##### SET DELTA ###
				*delta = ( l / nnn );
				//###################
}

void zero_init(int *i, float *L, float *bVol, float *sigma,\
	       int *N, float *n, float *delta, float *l,\
       int *m_rej, int *m_acc)
{
	*i = 0;
	*L = 0.0;
	*bVol = 0.0;
	*sigma = 0.0;
	*N = 0.0;
	*n = 0.0;
	*delta = 0.0;
	*l = 0.0;
	*m_rej = 0;
	*m_acc = 0;
}

float x_i(int i, int N, float delta)
{
  float x_i;
	x_i = (float)(i / ((int)sqrtf((float)N)));
  	x_i *= delta;
	return x_i;
}

float y_i(int i, int N, float delta)
{
  float y_i;
	y_i = fmodf(((float)i),  sqrtf((float)N));
	y_i *= delta;
	return y_i;
}

void rand_gen(void)
{
	srand(time(NULL));
	int r = rand();
	r = ( (2*r) + 1 );
	init_genrand(r);
}

float mc_dist(float x, float y, float x_k,\
		float y_k, float L)
{
	float x2, y2, dist;
	x2 = fabsf(x-x_k);
	while(x2>0.5*L) x2=L-x2;
	y2 = fabsf(y-y_k);
	while(y2>0.5*L) y2=L-y2;
	
	dist = sqrtf( x2*x2 + y2*y2 );
	
	return dist;

}

void mc_move(int i, int n, int *m_rej, int *m_acc, float *x, float *y,\
		 float x_i, float y_i, float sigma,\
		  float l, float L)
{
  int k; 
  int rej, gre, gac;
  float xau, yau, rn1, rn2, dlx, dly, dist;
	xau = x_i; 
	yau = y_i; 
	gre = *m_rej; 
	gac = *m_acc; 
	
	rn1 = genrand_real1();
	rn2 = genrand_real1();
	
	dlx = ( xau + (2*(rn1 - 0.5)*l) );
	dly = ( yau + (2*(rn2 - 0.5)*l) );
	
	while(dlx>=L) dlx-=L;
	while(dlx<0.0)dlx+=L;
	while(dly>=L) dly-=L;
	while(dly<0.0)dly+=L;
	
	rej = 0;
	for(k=0; ((k<n)&&(rej==0)); k++){
		if( k != i ){
			dist = mc_dist(dlx,dly,x[k],y[k],L);
			if( dist <= sigma ){ rej++; }
		}
	}
	if( rej==0 ){
		x[i] = dlx;
		y[i] = dly;
		gac++;
	}
	else gre++;

	*m_acc = gac;
	*m_rej = gre;
}


void mc_hist(float rmax, float L, int N,\
		int *histogram, float deltaR, float *x,\
		float *y)
{
	int i, j, bin;
	float dist, x2, y2;

	for(i=0; i<N; i++){
	 for(j=0; j<i; j++){
		x2 = fabsf( x[i] - x[j] );
		while(x2>(0.5*L)) x2-=L-x2;
		y2 = fabsf( y[i] - y[j] );
		while(y2>(0.5*L)) y2-=L-y2;
		dist=sqrtf(x2*x2 + y2*y2);
		bin = (int)( (dist/rmax)* ((float)B) );
		if( bin < B ) histogram[bin] += 1;
	 }
	}
}
