/*
  Write a program to solve the follow PDE numerically
  △u=f(x,y)
  where (x,y) on [0,1]x[0,1], and boundary value are all zero.
  Suppose we have a uniform mesh h=1/N, and we give a f(x,y)=(-5.0)*(M_PI*M_PI)*(sin(M_PI*x))*(sin(2.0*M_PI*y)).
  By finite difference method, we can derive a matrix A such that
  Au=b.
  Please, we want to get the result PDE numerically by Multigrid.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

int dec_prime(int N);
int Segmentation(int N, double *b);
int initial(int NK, double *u);
int GSM_smooth(double *b ,double *u, double *x, int NC);
int Scratch_1(double *y, double *u, int NC ,int k);
int creat_residual(double *b, double *u, double *r, int NC);
int Newsol(int n, double *u);
int Restrict(double *r, double *x, int NK, int N);
int Interpolate(double *u,double *x, int N);
int MG_RES(double *u, double *b, double *r,double *y, int N, double *x, int k, double *w, int t);
int MG_INT(double *u, double *b, double *y, int N , double *x, double *z, int k1, double *w, int k3);
int Scratch_1(double *y, double *u, int NC ,int k);
int Scratch_2(double *w, double *r, int t, int NK);
int sum(double *u, double *x, int N, int k1);
int Post_smooth(double *b ,double *u, double *x, int NC);
int error(double *g, double *u, double *x, int N);

int main(void)
{
	int N, NC, i, NK, p, k, k1, j, k2, t, k3, c, itr;
	double *r, *u, *b, *x, *y, *z, *w, *g;
	clock_t t1, t2;
	printf("請輸入N(必須是2的次方數且大於4)：");         //輸入要切成幾等分 
	scanf("%d", &N);
	NC=N-1;
	NK=NC*NC;                                           //計算內點個數 
	
	p=dec_prime(N);                                     //將N質因數分解 
	
	k=N;
	k1=0;
	c=0;
	itr=0;                                              //將迭代次數設為0次 
	printf("未知的點數M為：%d\n", NK);
	printf("測試函數：f(x,y)=(-5.0)*(M_PI*M_PI)*(sin(M_PI*x))*(sin(2.0*M_PI*y))\n");
	for(i=0;i<p-1;i++)
	{
		k1=(k-1)*(k-1)+k1;
		k=k/2;
	} 
	k3=k1-(N-1)*(N-1);
	r = (double *) malloc(NK*sizeof(double));
	u = (double *) malloc(NK*sizeof(double));
	b = (double *) malloc(NK*sizeof(double));
	x = (double *) malloc(NK*sizeof(double));
	y = (double *) malloc(k1*sizeof(double));
	z = (double *) malloc(NK*sizeof(double));
	w = (double *) malloc(k3*sizeof(double));
	g = (double *) malloc(NK*sizeof(double));
	initial(NK, u);                                     //將要進行第一次迭代前設定初始值為0向量 
	Segmentation(N, b);                                 //將(1.0/(N*N))b計算出來 
	
	for(j=0;j<NK;j++)
	{
		z[j]=b[j];                                      //因為b向量在運算中會改變，因此用z向量暫存 
		g[j]=u[j];                                      //暫存第n次迭代得到的解u向量 
	}
	
	t1 = clock();                                       //計算迭代時間 
	
	while(c==0)
	{
		k2=0;
		t=0;
		
		MG_RES(u, b, r, y, N, x, k2, w, t);             //將細網格轉換成粗網格減少我們的維度 
		MG_INT(u, b, y, N, x, z, k1, w, k3);            //將粗網格轉換成細網格讓我得到一次迭代後的解 
		                                                //完整一次的MG_RES+MG_INT叫做V-cycle 
		c=error(g, u, x, N);                            //讓此迭代停止的條件函數 
		
		for(j=0;j<NK;j++)
		{
			b[j]=z[j];                                  //因為一次迭代後，若要再一次，要把原問題的b向量代回去 
			g[j]=u[j];                                  //暫存第n次迭代得到的解u向量 
		}
		itr=itr+1;                                      //計算迭代次數 
	}
	
	t2 = clock();                                       //計算迭代時間 
	
	printf("Fast Poisson Solver By Multigrid: %f secs", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	printf("\n收斂次數：%d\n", itr);
	
	free(u);                                            //釋放記憶體 
	free(r);
	free(b);
	free(x);
	free(y);
	free(z);
	free(w);
	free(g);
}

int dec_prime(int N)                                    //將N質因數分解的函數 
{
	int p;
	p = 0;
	
	while(N%2 == 0)
	{
		N = N/2;
		p += 1;
	}
	
	return p;
}

int initial(int NK, double *u)                          //將要進行第一次迭代前設定初始值為0向量的函數 
{
	int i;
	
	for(i=0;i<NK;i++)
	{
		u[i]=0.0;
	}
}

int Segmentation(int N, double *b)                      //將(1.0/(N*N))b計算出來的函數 
{
	int i, j, t;
	double k;
	t=0;
	k=(N-1.0)/N;
	
	for(j=0;j<N-1;j++)
	{
		for(i=0;i<N-1;i++)
		{
				b[t]=(1.0/(N*N))*(-5.0)*(M_PI*M_PI)*(sin(M_PI*k))*(sin(2.0*M_PI*(i+1.0)/N));
				t=t+1;
		}
		k=k-(1.0/N);
	}
}

int MG_RES(double *u, double *b, double *r,double *y, int N , double *x, int k, double *w, int t)   //將細網格轉換成粗網格減少我們的維度函數 
{
	int i, NK, k1;
	
	if (N==2)                                                                                       //停止條件 
	{
		u[0]=-0.25*b[0];                                                                            //最粗的網格直接解 
	}
	else
	{
		for(i=0;i<3;i++)
		{
			GSM_smooth(b, u, x, N-1);                                                               //進行smooth，可以降低的高頻率，這裡進行3次 
		}

		Scratch_1(y, u, N-1, k);                                                                    //將每一次得到的u向量都暫存起來 
		k=(N-1)*(N-1)+k;
		creat_residual(b, u, r, N-1);                                                               //計算(1.0/(N*N))b-Ax=r向量(為殘差) 
		N = (N/2);                                                                                  //降低維度 
    	NK = (N-1)*(N-1);                                                                           //並且計算此時的粗網格內點個數 
    	Newsol(NK, u);                                                                              //此時令粗網格初始值為0向量 
		Restrict(r, x, NK, 2*N);                                                                    //利用細網格的r向量估計粗網格的r向量 
    	for(i=0;i<(N-1)*(N-1);i++)
		{
			r[i]=x[i];
		}
		
		if (N != 2)                                                                                 //當還未最粗網格時的做法 
		{
			Scratch_2(w, r, t, NK);                                                                 //將每一次得到的r向量都暫存起來
			t = NK+t;

			for(i=0;i<(N-1)*(N-1);i++)
			{	
				b[i]=r[i];
			}
			MG_RES(u, b, r, y, N, x, k, w, t);                                                      //繼續讓網格更粗 
		}
		else                                                                                        //達最粗網格時的做法  
		{
			for(i=0;i<(N-1)*(N-1);i++)
			{
			
				b[i]=r[i];
			}
		
			MG_RES(u, b, r, y, N, x, k, w, t);                                                      //再進行最後一次的MG_RES直接解   
		}
	}
} 

int MG_INT(double *u, double *b, double *y, int N , double *x, double *z, int k1, double *w, int k3)//將粗網格轉換成細網格讓我得到一次迭代後的解 
{
	int N1,i;
	N1=N;
	N=2;

	while(N<N1)                                                                                     //停止條件 
	{
		Interpolate(u, x, N);                                                                       //從最粗網格得到的解，往回轉成細網格，至與原題目粗細大小一樣為止                                                                     
		N = N*2;
		
		for(i=0;i<(N-1)*(N-1);i++)
		{
			u[i]=x[i];
		}
		
		k1=k1-(N-1)*(N-1);
    	sum(u, y, N, k1);                                                                           //將估計的殘差加回MG_RES時得到的smooth解，修正解的值 
    	
		if (N==N1)                                                                                  //如果已達與原題目網格粗細大小一樣的做法 
    	{
		 	for(i=0;i<(N1-1)*(N1-1);i++)
		 	{
    			b[i]=z[i];
			}
			for(i=0;i<3;i++)
			{
				GSM_smooth(b, u, x, N1-1);                                                          //修正後，再smooth，三次 
			}
		} 
		else                                                                                        //若未達與原題目網格粗細大小一樣的做法 
		{
			k3=k3-(N-1)*(N-1);
			for(i=0;i<(N-1)*(N-1);i++)
			{
				b[i]=w[i+k3]; 
			}
			for(i=0;i<3;i++)
			{
				GSM_smooth(b, u, x, N-1);                                                           //修正後，再smooth，三次  
			}
		}
	}			
}

int GSM_smooth(double *b ,double *u, double *x, int NC)                                             //smooth的函數，這裡是Gauss Seidel Method 
{
	int i, j, k, NK;
	double w;
	NK=NC*NC;
	k=0;
	w=1.0;
	
	u[0]=-0.25*w*(b[0]-u[1]-u[NC]);
	for(i=1;i<NC-1;i++)
	{
		u[i]=-0.25*w*(b[i]-u[i-1]-u[i+1]-u[NC+i]);
	}
	u[NC-1]=-0.25*w*(b[NC-1]-u[NC-2]-u[2*NC-1]);
	
	while(k<=(NC-3))
	{
		u[NC*(k+1)]=-0.25*w*(b[NC*(k+1)]-u[NC*(k+1)-NC]-u[NC*(k+1)+1]-u[NC*(k+1)+NC]);
		for(j=1;j<NC-1;j++)
		{
			u[NC*(k+1)+j]=-0.25*w*(b[NC*(k+1)+j]-u[NC*(k+1)+j-1]-u[NC*(k+1)+j+1]-u[NC*(k+1)+j+NC]-u[NC*(k+1)+j-NC]);
		}
		u[NC*(k+1)+j]=-0.25*w*(b[NC*(k+1)+j]-u[NC*(k+1)+j-1]-u[NC*(k+1)+j-NC]-u[NC*(k+1)+j+NC]);
		k=k+1;
	}
	
	u[NK-NC]=-0.25*w*(b[NK-NC]-u[NK-NC+1]-u[NC]);
	for(i=NK-NC+1;i<NK-1;i++)
	{
		u[i]=-0.25*w*(b[i]-u[i-1]-u[i+1]-u[i-NC]);
	}
	u[NK-1]=-0.25*w*(b[NK-1]-u[NK-2]-u[NK-NC-1]);	
}

int creat_residual(double *b, double *u, double *r ,int NC)                                         //計算(1.0/(N*N))b-Ax=r向量(為殘差)函數 
{
	int i, j, k, t;
	k=NC;
	t=k*k;
	
	r[0]=b[0]+(4*u[0]-u[1]-u[k]);
	
	for(i=1;i<k-1;i++)
	{
		r[i]=b[i]+(-u[i-1]+4*u[i]-u[i+1]-u[i+k]);
	}
	
	r[i]=b[i]+(-u[i-1]+4*u[i]-u[i+k]);
	i=k;
	
	while (i<=(t-2*k))
	{
		r[i]=b[i]+(-u[i-k]+4*u[i]-u[i+1]-u[i+k]);
		i=i+k;
	}
	
	i=2*k-1;
	
	while (i<t-k)
	{
		r[i]=b[i]+(-u[i-k]-u[i-1]+4*u[i]-u[i+k]);
		i=i+k;
	}
	
	i=k+1;
	
	while(i<(t-k))
	{
		for(j=0;j<k-2;j++)
		{
			r[i+j]=b[i+j]+(-u[i-k+j]-u[i-1+j]+4*u[i+j]-u[i+1+j]-u[i+k+j]);
		}
		i=i+k;
	}
	
	r[t-1]=b[t-1]+(4*u[t-1]-u[t-2]-u[t-k-1]);
	r[t-k]=b[t-k]+(4*u[t-k]-u[t-k+1]-u[t-k*2]);
	
	for(i=t-k+1;i<t-1;i++)
	{
		r[i]=b[i]+(-u[i-k]-u[i-1]+4*u[i]-u[i+1]);
	}
	
}

int Newsol(int n, double *u)                                                                        //此時令粗網格初始值為0向量函數 
{
	int i;
	
	for(i=0;i<n;i++)
	{
		u[i]=0.0;
	} 
}

int Restrict(double *r, double *x, int NK, int N)                                                   //利用細網格的r向量估計粗網格的r向量函數 
{
	int i, j, k, m, N1, n;
	N1=(N/2)-1;
	k=0;
	n=0;
	
	for(m=0;m<N1;m++)
	{
		j=0;
		
		for(i=0;i<N1;i++)
		{
			x[n]=4*(0.25*0.25*r[j+k]+0.5*0.25*r[j+1+k]+0.25*0.25*r[j+2+k]+0.5*0.25*r[N+j-1+k]+0.25*r[N+j+k]+0.5*0.25*r[N+j+1+k]+0.25*0.25*r[2*N+j-2+k]+0.5*0.25*r[2*N+j-1+k]+0.25*0.25*r[2*N+j+k]);
			j=j+2;
			n=n+1;
		}
		
		k=k+N+j;
	}	
}

int Interpolate(double *u, double *x, int N)                                                        //從最粗網格得到的解，往回轉成估計細網格值的函數 
{
	int i, k, j, N1, N2, NK, g , h, N4, N3;
	NK=N-2;
	N1=N*2-1;
	N3=N1*N1;
	N2=N-1;
	N4=N2*N2;
	k=0;
	
	if(N>2)
	{
		x[0]=0.25*u[0];
		for(i=0;i<NK;i++)
		{
			x[k+1]=0.5*u[i];
			x[k+2]=0.25*u[i]+0.25*u[i+1];
			k=k+2;
		}
	
		x[k+1]=0.5*u[i];
		x[k+2]=0.25*u[i];

		k=0;
		x[N3-N1]=0.25*u[N4-N2];
	
		for(i=0;i<NK;i++)
		{
			x[N3-N1+k+1]=0.5*u[N4-N2+i];
			x[N3-N1+2+k]=0.25*u[N4-N2+i]+0.25*u[N4-N2+1+i];
			k=k+2;
		}
	
		x[N3-2]=0.5*u[N4-1];
		x[N3-1]=0.25*u[N4-1];
	
		j=1;
		k=0;
		h=0;
		g=0;
	
		while(j<N1)
		{
			x[j*N1]=0.5*u[N2*h];
		
			for(i=0;i<NK;i++)
			{
				x[j*N1+1+g]=u[h*N2+i];
				x[j*N1+2+g]=0.5*u[h*N2+i]+0.5*u[h*N2+1+i];
				g=g+2;
			}
		
			h=h+1;
			x[j*N1+g+1]=u[h*N2-1];
			x[j*N1+g+2]=0.5*u[h*N2-1];
			j=j+2;
			g=0;
		}
	
		j=2;
		k=0;
		h=0;
		g=0;
	
		while(j<N1-1)
		{
			x[j*N1]=0.25*u[N2*h]+0.25*u[N2*h+N2];
		
			for(i=0;i<NK;i++)
			{
				x[j*N1+1+g]=0.5*u[h*N2+i]+0.5*u[h*N2+i+N2];
				x[j*N1+2+g]=0.25*u[h*N2+i]+0.25*u[h*N2+1+i]+0.25*u[h*N2+i+N2]+0.25*u[h*N2+1+i+N2];
				g=g+2;
			}
		
			h=h+1;
			x[j*N1+g+1]=0.5*u[h*N2-1]+0.5*u[h*N2-1+N2];
			x[j*N1+g+2]=0.25*u[h*N2-1]+0.25*u[h*N2-1+N2];
			j=j+2;
			g=0;
		}
	}
	else
	{
		x[0]=0.25*u[0];
		x[1]=0.5*u[0];
		x[2]=0.25*u[0];
		x[3]=0.5*u[0];
		x[4]=u[0];
		x[5]=0.5*u[0];
		x[6]=0.25*u[0];
		x[7]=0.5*u[0];
		x[8]=0.25*u[0];
	}		
}

int sum(double *u, double *y, int N , int k1)                                                       //將估計的殘差加回MG_RES時得到的smooth解，修正解的值(函數) 
{
	int i, NK;
	NK=(N-1)*(N-1);
	
	for(i=0;i<NK;i++)
	{
		u[i]=u[i]+y[i+k1];
	}
}

int Scratch_1(double *y, double *u, int NC ,int k)                                                  //將每一次得到的u向量都暫存起來 的函數 
{
	int i, NK;
	NK=NC*NC;
	
	for(i=0;i<NK;i++)
	{
		y[i+k]=u[i];	
	}
}

int Scratch_2(double *w, double *r, int t, int NK)                                                  //將每一次得到的r向量都暫存起來 
{
	int i;
	
	for(i=0;i<NK;i++)
	{
		w[i+t]=r[i];	
	}
}

int error(double *g, double *u, double *x, int N)                                                   //讓此迭代停止的條件函數
{
	int i;
	
	for(i=0;i<(N-1)*(N-1);i++)
	{
		x[i]=fabs(g[i]-u[i]);
	}
	
	for(i=0;i<(N-1)*(N-1);i++)
	{
        if(x[0]<x[i])
		{
            x[0] = x[i];
        }       
    }
    
    for(i=0;i<(N-1)*(N-1);i++)
	{
        if(fabs(u[0])<fabs(u[i]))
		{
            x[1] = fabs(u[i]);
        } 
		else
		{
			x[1] = fabs(u[0]);
		}      
    }

    if(x[0]/x[1] < 0.0000001)
    {
    	return 1;
	}
	else
	{
		return 0;
	}
}
