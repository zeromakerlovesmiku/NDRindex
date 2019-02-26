#include<Rcpp.h>
#include<cstdio>
#include<cmath>
#include<cstring>
#include<algorithm>
#include<vector>
#include<ctime>
#include<string>
#include<string.h>
#include<time.h>
using namespace Rcpp;
using namespace std;
double dis(vector<double>A,vector<double>B,int m){
  double tem=0;
  for (int j=0;j<m;j++)
    tem+=(A[j]-B[j])*(A[j]-B[j]);
  return sqrt(tem);
}

// [[Rcpp::export]]
double NDRindex(NumericMatrix A,int n,int m){
    double x,y,o,oo=0,are,z,S,inf,ans,sum,sumo;int t,id;
    vector<double>tempdis;
    int tempdiscnt=0;
    vector<vector<double> >a,aa,q,icenter;
    vector<double>b,p;
    tempdis.clear();
    a.clear();
    aa.clear();
    q.clear();
    icenter.clear();
    b.clear();
    p.clear();
    for (int i=0;i<n;i++){
        b.clear();
        for (int j=0;j<m;j++)b.push_back(A(i,j));
        aa.push_back(b);
    }
    ans=0;
    tempdis.clear();
    for (int i=0;i<n-1;i++)
        for (int j=i+1;j<n;j++)
            tempdis.push_back(dis(aa[i],aa[j],m));
    sort(tempdis.begin(),tempdis.end());
    tempdiscnt=tempdis.size();
    if (tempdiscnt%4==0)
        S=tempdis[tempdiscnt/4];
    else
        S=(tempdis[tempdiscnt/4]+tempdis[(tempdiscnt+(4-tempdiscnt%4))/4])/2;
    S/=log10(n);
    for (int pp=1;pp<=min(n,100);pp++){
        a=aa;
        t=1;
        //srand(time(0));
        int sd;
        sd=pp%n;
        p=a[sd];q.clear();
        q.push_back(a[sd]);
        a.erase(a.begin()+sd);
        sum=1;sumo=0;
        inf=1000000000;
        while(1){
            o=inf;
            for (int i=0;i<a.size();i++)
                if((z=dis(p,a[i],m))<o)
                    o=z,id=i;
            if (o==inf)break;
            if (o<S){
                q.push_back(a[id]);t++;
                for (int j=0;j<m;j++)p[j]=0;

                for (int i=0;i<q.size();i++)
                    for (int j=0;j<m;j++)p[j]+=q[i][j];

                for (int j=0;j<m;j++)p[j]=p[j]/t;
                a.erase(a.begin()+id);
            }
            else {
                are=0;
                for (int i=0;i<t;i++)are+=dis(q[i],p,m);
                if (t>3)
                {
                    are=are/t;
                    sumo+=are;
                    sum++;
                    icenter.push_back(p);
                }
                q.clear();
                q.push_back(a[id]);
                t=1;
                p=a[id];
            }
        }
        oo+=(1.0-((1.0*sumo)/(S*sum)));
    }
    return oo/100;
}
