/*
  Write a program to solve the follow PDE numerically
  ��u=f(x,y)
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
	printf("�п�JN(�����O2������ƥB�j��4)�G");         //��J�n�����X���� 
	scanf("%d", &N);
	NC=N-1;
	NK=NC*NC;                                           //�p�⤺�I�Ӽ� 
	
	p=dec_prime(N);                                     //�NN��]�Ƥ��� 
	
	k=N;
	k1=0;
	c=0;
	itr=0;                                              //�N���N���Ƴ]��0�� 
	printf("�������I��M���G%d\n", NK);
	printf("���ը�ơGf(x,y)=(-5.0)*(M_PI*M_PI)*(sin(M_PI*x))*(sin(2.0*M_PI*y))\n");
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
	initial(NK, u);                                     //�N�n�i��Ĥ@�����N�e�]�w��l�Ȭ�0�V�q 
	Segmentation(N, b);                                 //�N(1.0/(N*N))b�p��X�� 
	
	for(j=0;j<NK;j++)
	{
		z[j]=b[j];                                      //�]��b�V�q�b�B�⤤�|���ܡA�]����z�V�q�Ȧs 
		g[j]=u[j];                                      //�Ȧs��n�����N�o�쪺��u�V�q 
	}
	
	t1 = clock();                                       //�p�⭡�N�ɶ� 
	
	while(c==0)
	{
		k2=0;
		t=0;
		
		MG_RES(u, b, r, y, N, x, k2, w, t);             //�N�Ӻ����ഫ���ʺ����֧ڭ̪����� 
		MG_INT(u, b, y, N, x, z, k1, w, k3);            //�N�ʺ����ഫ���Ӻ������ڱo��@�����N�᪺�� 
		                                                //����@����MG_RES+MG_INT�s��V-cycle 
		c=error(g, u, x, N);                            //�������N��������� 
		
		for(j=0;j<NK;j++)
		{
			b[j]=z[j];                                  //�]���@�����N��A�Y�n�A�@���A�n�����D��b�V�q�N�^�h 
			g[j]=u[j];                                  //�Ȧs��n�����N�o�쪺��u�V�q 
		}
		itr=itr+1;                                      //�p�⭡�N���� 
	}
	
	t2 = clock();                                       //�p�⭡�N�ɶ� 
	
	printf("Fast Poisson Solver By Multigrid: %f secs", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	printf("\n���Ħ��ơG%d\n", itr);
	
	free(u);                                            //����O���� 
	free(r);
	free(b);
	free(x);
	free(y);
	free(z);
	free(w);
	free(g);
}

int dec_prime(int N)                                    //�NN��]�Ƥ��Ѫ���� 
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

int initial(int NK, double *u)                          //�N�n�i��Ĥ@�����N�e�]�w��l�Ȭ�0�V�q����� 
{
	int i;
	
	for(i=0;i<NK;i++)
	{
		u[i]=0.0;
	}
}

int Segmentation(int N, double *b)                      //�N(1.0/(N*N))b�p��X�Ӫ���� 
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

int MG_RES(double *u, double *b, double *r,double *y, int N , double *x, int k, double *w, int t)   //�N�Ӻ����ഫ���ʺ����֧ڭ̪����ר�� 
{
	int i, NK, k1;
	
	if (N==2)                                                                                       //������� 
	{
		u[0]=-0.25*b[0];                                                                            //�̲ʪ����檽���� 
	}
	else
	{
		for(i=0;i<3;i++)
		{
			GSM_smooth(b, u, x, N-1);                                                               //�i��smooth�A�i�H���C�����W�v�A�o�̶i��3�� 
		}

		Scratch_1(y, u, N-1, k);                                                                    //�N�C�@���o�쪺u�V�q���Ȧs�_�� 
		k=(N-1)*(N-1)+k;
		creat_residual(b, u, r, N-1);                                                               //�p��(1.0/(N*N))b-Ax=r�V�q(���ݮt) 
		N = (N/2);                                                                                  //���C���� 
    	NK = (N-1)*(N-1);                                                                           //�åB�p�⦹�ɪ��ʺ��椺�I�Ӽ� 
    	Newsol(NK, u);                                                                              //���ɥO�ʺ����l�Ȭ�0�V�q 
		Restrict(r, x, NK, 2*N);                                                                    //�Q�βӺ��檺r�V�q���p�ʺ��檺r�V�q 
    	for(i=0;i<(N-1)*(N-1);i++)
		{
			r[i]=x[i];
		}
		
		if (N != 2)                                                                                 //���٥��̲ʺ���ɪ����k 
		{
			Scratch_2(w, r, t, NK);                                                                 //�N�C�@���o�쪺r�V�q���Ȧs�_��
			t = NK+t;

			for(i=0;i<(N-1)*(N-1);i++)
			{	
				b[i]=r[i];
			}
			MG_RES(u, b, r, y, N, x, k, w, t);                                                      //�~���������� 
		}
		else                                                                                        //�F�̲ʺ���ɪ����k  
		{
			for(i=0;i<(N-1)*(N-1);i++)
			{
			
				b[i]=r[i];
			}
		
			MG_RES(u, b, r, y, N, x, k, w, t);                                                      //�A�i��̫�@����MG_RES������   
		}
	}
} 

int MG_INT(double *u, double *b, double *y, int N , double *x, double *z, int k1, double *w, int k3)//�N�ʺ����ഫ���Ӻ������ڱo��@�����N�᪺�� 
{
	int N1,i;
	N1=N;
	N=2;

	while(N<N1)                                                                                     //������� 
	{
		Interpolate(u, x, N);                                                                       //�q�̲ʺ���o�쪺�ѡA���^�ন�Ӻ���A�ܻP���D�زʲӤj�p�@�ˬ���                                                                     
		N = N*2;
		
		for(i=0;i<(N-1)*(N-1);i++)
		{
			u[i]=x[i];
		}
		
		k1=k1-(N-1)*(N-1);
    	sum(u, y, N, k1);                                                                           //�N���p���ݮt�[�^MG_RES�ɱo�쪺smooth�ѡA�ץ��Ѫ��� 
    	
		if (N==N1)                                                                                  //�p�G�w�F�P���D�غ���ʲӤj�p�@�˪����k 
    	{
		 	for(i=0;i<(N1-1)*(N1-1);i++)
		 	{
    			b[i]=z[i];
			}
			for(i=0;i<3;i++)
			{
				GSM_smooth(b, u, x, N1-1);                                                          //�ץ���A�Asmooth�A�T�� 
			}
		} 
		else                                                                                        //�Y���F�P���D�غ���ʲӤj�p�@�˪����k 
		{
			k3=k3-(N-1)*(N-1);
			for(i=0;i<(N-1)*(N-1);i++)
			{
				b[i]=w[i+k3]; 
			}
			for(i=0;i<3;i++)
			{
				GSM_smooth(b, u, x, N-1);                                                           //�ץ���A�Asmooth�A�T��  
			}
		}
	}			
}

int GSM_smooth(double *b ,double *u, double *x, int NC)                                             //smooth����ơA�o�̬OGauss Seidel Method 
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

int creat_residual(double *b, double *u, double *r ,int NC)                                         //�p��(1.0/(N*N))b-Ax=r�V�q(���ݮt)��� 
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

int Newsol(int n, double *u)                                                                        //���ɥO�ʺ����l�Ȭ�0�V�q��� 
{
	int i;
	
	for(i=0;i<n;i++)
	{
		u[i]=0.0;
	} 
}

int Restrict(double *r, double *x, int NK, int N)                                                   //�Q�βӺ��檺r�V�q���p�ʺ��檺r�V�q��� 
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

int Interpolate(double *u, double *x, int N)                                                        //�q�̲ʺ���o�쪺�ѡA���^�ন���p�Ӻ���Ȫ���� 
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

int sum(double *u, double *y, int N , int k1)                                                       //�N���p���ݮt�[�^MG_RES�ɱo�쪺smooth�ѡA�ץ��Ѫ���(���) 
{
	int i, NK;
	NK=(N-1)*(N-1);
	
	for(i=0;i<NK;i++)
	{
		u[i]=u[i]+y[i+k1];
	}
}

int Scratch_1(double *y, double *u, int NC ,int k)                                                  //�N�C�@���o�쪺u�V�q���Ȧs�_�� ����� 
{
	int i, NK;
	NK=NC*NC;
	
	for(i=0;i<NK;i++)
	{
		y[i+k]=u[i];	
	}
}

int Scratch_2(double *w, double *r, int t, int NK)                                                  //�N�C�@���o�쪺r�V�q���Ȧs�_�� 
{
	int i;
	
	for(i=0;i<NK;i++)
	{
		w[i+t]=r[i];	
	}
}

int error(double *g, double *u, double *x, int N)                                                   //�������N���������
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
