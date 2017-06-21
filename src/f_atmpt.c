#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declare.h"
#define pi 4.0*atanf(1.0)

int main ()
{
  char opc;
  int N, i, j, mcm, m_rej, m_acc;
  float L, bVol, sigma, Nfloat, delta, n, l, rmax, deltaR;

	zero_init(&i,&L,&bVol,&sigma,\
		&N,&n,&delta,&l,&m_rej,&m_acc);

// STEP 1 ####################################################
	
	data(&L,&N,&n,&bVol,&sigma,&delta);
	
	int *histogram = malloc( B*sizeof(int) );
	float *x = malloc( N*sizeof(float) );
	float *y = malloc( N*sizeof(float) );

	// POSITIONS INITIALIZATION ###########
	FILE *fp;
	fp=fopen("initial_snapshot.dat","w");
	for(i=0; i<N; i++){
		x[i] = x_i(i,N,delta);
		y[i] = y_i(i,N,delta);
		fprintf(fp,"%lf %lf \n",x[i],y[i]);
	}

	printf("Verbose particle initialization? (y/n): ");
	  scanf(" %c",&opc);
	if( opc=='y'){
		for( i=0; i<N; i++){
			printf("x[%d] = %lf \t y[%d] = %lf\n",
					i,x[i],i,y[i]);
		}
	}
//############################################################

	rand_gen();
	l = ( L / 10.0 );


// STEP 2 ####################################################

	rmax = L/2.0;


	printf("Number of epochs: ");
	  scanf("%d",&mcm);
		m_acc=0;
		m_rej=0;
	for(j=0; j<mcm; j++){
		for(i=0; i<N; i++)
			mc_move(i,N,&m_rej,&m_acc,x,y,
				x[i],y[i],sigma,l,L);
		if( (j%1000) == 0 ){
			printf("mv accepted = %d \n",m_acc);
			printf("mv rejected = %d \n\n",m_rej);
		}
		mc_hist(rmax, L, N, histogram, deltaR, x, y);
	}
	// Data to file ###############
	FILE *hst;
	hst=fopen("histogram.dat","w");
	float lol, nono;
	float dr = rmax/(float)B;
	for(i=1; i<B; i++){
		lol = ((float)i) * dr;
		nono = (float)histogram[i];
		nono *= (2.0*L*L)/((N-1.0)*N*2.0*pi*(float)i*dr*dr*(float)mcm);
		fprintf(hst,"%lf %lf \n",lol,nono);
	}
	fclose(hst);

	FILE *mv;
	mv=fopen("mc_moves.dat","w");
	for(i=0; i<N; i++)
		fprintf(mv,"%lf %lf \n",x[i],y[i]);
//############################################################

	// CLOSING AND FREEING MEMO ###########################
	free(histogram);
	free(x);
	free(y);
	fclose(fp);
	fclose(mv);
	fclose(hst);
	//##########
}

