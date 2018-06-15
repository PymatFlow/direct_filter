#include "spatial.h"
//#include "ls.c"

double var(double y[],int n)
{
    double mean,s;
    int i;
    
    s=0;
    for(i=0;i<n;i++)
        s+=y[i];
    mean=s/n;
    s=0;
    for(i=0;i<n;i++)
        s+=(y[i]-mean)*(y[i]-mean);
    return(sqrt(s/n));
}

void direct_filter(double **orig,double **noise,double **pred,int numrows,int numcols)
{
    int i,j,k,l,sum;
    unsigned char **flag;
    double **C,e1,e2;
    double buf,tmp,w[K],r[K],y[MX];
    FILE *diafile;
    
    init(3,1);
    C=mem_allocate(K,MX);
    flag=mem_assign(H,W);
    e1=0;
    e2=0;
    sum=0;
    //directional filtering
    
    for(i=0;i<numrows;i++)
    {
        printf("%d\n",i);
        for(j=0;j<numcols;j++)
        {
            if(i>=B+T&&i<numrows-T-B&&j>=B+T&&j<numcols-T-B)
            {
                for(k=0;k<MX;k++)
                {y[k]=pred[i+nx[k]][j+ny[k]];
                 for(l=0;l<K;l++)
                     C[l][k]=pred[i+nx[k]+mx[l]][j+ny[k]+my[l]];
                }
                for(k=0;k<K;k++)
                    r[k]=pred[i+mx[k]][j+my[k]];
                if(var(r,K)<th)
                {
                    buf=0;
                    for(k=0;k<K;k++)
                        buf+=pred[i+mx[k]][j+my[k]]/K;
                    if(fabs(noise[i][j]-buf)>Th)
                    {noise[i][j]=buf;flag[i][j]=1;sum++;}
                }
                else
                {
                    least_sq(C,y,w,K,MX);
                    buf=0;
                    for(k=0;k<K;k++)
                        buf+=r[k]*w[k];
                    if(fabs(noise[i][j]-buf)>Th)
                    {noise[i][j]=buf;flag[i][j]=1;sum++;}
                }
            }
            else
                noise[i][j]=pred[i][j];
            e1+=(orig[i][j]-pred[i][j])*(orig[i][j]-pred[i][j]);
            e2+=(orig[i][j]-noise[i][j])*(orig[i][j]-noise[i][j]);
        }
    }
    //printf("original MSE=%lf\n",e1/(numrows*numcols));
    //printf("reduced MSE=%lf\n",e2/(numrows*numcols));
    printf("sum=%d, detected noise probability = %lf\n",sum,
            (double) sum/(numrows*numcols));
    diafile=fopen("dia","wb");
    for(i=0;i<numrows;i++)
        fwrite(flag[i],sizeof(unsigned char),numcols,diafile);
    fclose(diafile);
    mem_free(C,K,MX);
    mem_free(flag,H,W);
}
