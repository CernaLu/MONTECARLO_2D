#define B 248
// //########################################### FUNCTIONS ########
void data(float *L, int *N, float *n,
		float *bvol, float *sigma, float *delta);
float x_i(int i, int N, float delta);
float y_i(int i, int N, float delta);
float mc_dist(float x, float y, float x_k,\
		float y_k, float L);

void zero_init(int *i, float *L, float *bVol,\
		float *sigma, int *N,\
		float *n, float *delta, float *l,\
		int *m_rej, int *m_acc);

void set_delta(float *delta, float L, int N);
void mc_move(int i, int n, int *m_rej, int *m_acc, float *x, float *y,\
		float x_i, float y_i, float sigma,\
		float l, float L);
void mc_hist(float rmax, float L, int N,\
		int *histogram, float deltaR, float *x,\
		float *y);
// //##################################################################


//################################## RANDGEN ########
void rand_gen(void);
void init_genrand(unsigned long s);
float genrand_real1(void);
//############################################

