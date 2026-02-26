#include <stdio.h>
#include <math.h>
#include <cblas.h>

#define EPS 1e-6

#define ASSERT_D(name, exp, got) \
    do { double _e=(double)(exp), _g=(double)(got); \
    if(fabs(_e-_g)>EPS){printf("[%s] FAIL: expected %.6f got %.6f\n",name,_e,_g);fails++;}\
    else printf("[%s] OK\n",name);} while(0)

#define ASSERT_I(name, exp, got) \
    do { long _e=(long)(exp), _g=(long)(got); \
    if(_e!=_g){printf("[%s] FAIL: expected %ld got %ld\n",name,_e,_g);fails++;}\
    else printf("[%s] OK\n",name);} while(0)

int main()
{
    int fails = 0;

    
    { double x[]={1,2,3}, y[]={4,5,6};
      ASSERT_D("ddot",32.0,cblas_ddot(3,x,1,y,1)); }

    { float x[]={1,2,3}, y[]={4,5,6};
      ASSERT_D("sdot",32.0f,cblas_sdot(3,x,1,y,1)); }

    
    {
        float x[]={1,0, 2,0};
        float y[]={3,0, 4,0};
        float res[2];
        cblas_cdotu_sub(2,x,1,y,1,res);
        ASSERT_D("cdotu",11.0,res[0]);
    }

    {
        double x[]={1,0, 2,0};
        double y[]={3,0, 4,0};
        double res[2];
        cblas_zdotu_sub(2,x,1,y,1,res);
        ASSERT_D("zdotu",11.0,res[0]);
    }

    
    { double x[]={1,2,3}, y[]={4,5,6};
      cblas_daxpy(3,2.0,x,1,y,1);
      ASSERT_D("daxpy",6.0,y[0]); }

    { float x[]={1,2,3}, y[]={4,5,6};
      cblas_saxpy(3,2.0f,x,1,y,1);
      ASSERT_D("saxpy",6.0f,y[0]); }

    
    { double x[]={1,2,3};
      cblas_dscal(3,3.0,x,1);
      ASSERT_D("dscal",3.0,x[0]); }

    { float x[]={1,2,3};
      cblas_sscal(3,3.0f,x,1);
      ASSERT_D("sscal",3.0f,x[0]); }

    
    { double x[]={1,2,3}, y[]={0};
      cblas_dcopy(3,x,1,y,1);
      ASSERT_D("dcopy",1.0,y[0]); }

    { float x[]={1,2,3}, y[]={0};
      cblas_scopy(3,x,1,y,1);
      ASSERT_D("scopy",1.0f,y[0]); }

    
    { double x[]={1,2}, y[]={3,4};
      cblas_dswap(2,x,1,y,1);
      ASSERT_D("dswap",3.0,x[0]); }

    { float x[]={1,2}, y[]={3,4};
      cblas_sswap(2,x,1,y,1);
      ASSERT_D("sswap",3.0f,x[0]); }

    
    { double x[]={3,4};
      ASSERT_D("dnrm2",5.0,cblas_dnrm2(2,x,1)); }

    { float x[]={3,4};
      ASSERT_D("snrm2",5.0f,cblas_snrm2(2,x,1)); }

    
    { double x[]={-1,2,-3};
      ASSERT_D("dasum",6.0,cblas_dasum(3,x,1)); }

    { float x[]={-1,2,-3};
      ASSERT_D("sasum",6.0f,cblas_sasum(3,x,1)); }

    
    { double x[]={1,5,3,9,2};
      ASSERT_I("idamax",3,cblas_idamax(5,x,1)); }

    { float x[]={1,5,3,9,2};
      ASSERT_I("isamax",3,cblas_isamax(5,x,1)); }

    
    { double x[]={5,-2,3,-1};
      ASSERT_I("idamin",3,cblas_idamin(4,x,1)); }

    { float x[]={5,-2,3,-1};
      ASSERT_I("isamin",3,cblas_isamin(4,x,1)); }

    
    { double x[]={1,0}, y[]={0,1};
      cblas_drot(2,x,1,y,1,0.0,1.0);
      ASSERT_D("drot",0.0,x[0]); }

    { float x[]={1,0}, y[]={0,1};
      cblas_srot(2,x,1,y,1,0.0f,1.0f);
      ASSERT_D("srot",0.0f,x[0]); }

    
    { double a=3,b=4,c,s;
      cblas_drotg(&a,&b,&c,&s);
      ASSERT_D("drotg",5.0,a); }

    { float a=3,b=4,c,s;
      cblas_srotg(&a,&b,&c,&s);
      ASSERT_D("srotg",5.0f,a); }

    
    {
        double d1=1,d2=1,x1=1,y1=1,param[5];
        cblas_drotmg(&d1,&d2,&x1,y1,param);
        int zero=1; for(int i=0;i<5;i++) if(param[i]!=0.0) zero=0;
        if(zero){ printf("[drotmg] FAIL\n"); fails++; }
        else printf("[drotmg] OK\n");
    }

    {
        float d1=1,d2=1,x1=1,y1=1,param[5];
        cblas_srotmg(&d1,&d2,&x1,y1,param);
        int zero=1; for(int i=0;i<5;i++) if(param[i]!=0.0f) zero=0;
        if(zero){ printf("[srotmg] FAIL\n"); fails++; }
        else printf("[srotmg] OK\n");
    }

    
    {
    double x[]={1,2}, y[]={3,4}, p[5]={-1,0,1,1,0};

    double x0=x[0], x1=x[1], y0=y[0], y1=y[1];

    cblas_drotm(2,x,1,y,1,p);

    int changed = (x[0]!=x0) || (x[1]!=x1) || (y[0]!=y0) || (y[1]!=y1);

    if(!changed){ printf("[drotm] FAIL\n"); fails++; }
    else printf("[drotm] OK\n");
    }   

    {
    float x[]={1,2}, y[]={3,4}, p[5]={-1,0,1,1,0};

    float x0=x[0], x1=x[1], y0=y[0], y1=y[1];

    cblas_srotm(2,x,1,y,1,p);

    int changed = (x[0]!=x0) || (x[1]!=x1) || (y[0]!=y0) || (y[1]!=y1);

    if(!changed){ printf("[srotm] FAIL\n"); fails++; }
    else printf("[srotm] OK\n");
}

    printf("\nTotal failed: %d\n",fails);
    return fails;
}
