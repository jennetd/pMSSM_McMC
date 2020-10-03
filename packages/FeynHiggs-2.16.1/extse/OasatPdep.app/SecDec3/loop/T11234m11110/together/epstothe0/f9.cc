#include "intfile.hh"

dcmplx Pf9(const double x[], double es[], double esx[], double em[], double lambda, double lrs[], double bi) {
double x0=x[0];
double x1=x[1];
double x2=x[2];
dcmplx y[77];
dcmplx FOUT;
dcmplx MYI(0.,1.);
y[1]=1./bi;
y[2]=em[0];
y[3]=x1*x1;
y[4]=x2*x2;
y[5]=esx[0];
y[6]=x0*x0;
y[7]=2.*x0*y[1]*y[2];
y[8]=-x0;
y[9]=1.+y[8];
y[10]=-x1;
y[11]=1.+y[10];
y[12]=2.*x1*y[1]*y[2];
y[13]=2.*x2*y[1]*y[2];
y[14]=2.*y[1]*y[2];
y[15]=4.*x0*x1*y[1]*y[2];
y[16]=4.*x0*x2*y[1]*y[2];
y[17]=-(x2*y[1]*y[5]);
y[18]=lrs[0];
y[19]=lrs[1];
y[20]=2.*y[1]*y[2]*y[6];
y[21]=-(x0*y[1]*y[5]);
y[22]=-(y[1]*y[5]*y[6]);
y[23]=y[7]+y[20]+y[21]+y[22];
y[24]=y[1]*y[2]*y[3];
y[25]=2.*x0*y[1]*y[2]*y[3];
y[26]=2.*x1*x2*y[1]*y[2];
y[27]=4.*x0*x1*x2*y[1]*y[2];
y[28]=y[1]*y[2]*y[4];
y[29]=2.*x0*y[1]*y[2]*y[4];
y[30]=-(x1*x2*y[1]*y[5]);
y[31]=-2.*x0*x1*x2*y[1]*y[5];
y[32]=y[12]+y[13]+y[17]+y[24]+y[25]+y[26]+y[27]+y[28]+y[29]+y[30]+y[31];
y[33]=-(lambda*MYI*y[9]*y[18]*y[32]);
y[34]=-x2;
y[35]=1.+y[34];
y[36]=-(y[1]*y[5]);
y[37]=-(x1*y[1]*y[5]);
y[38]=-2.*x0*x1*y[1]*y[5];
y[39]=y[12]+y[13]+y[14]+y[15]+y[16]+y[36]+y[37]+y[38];
y[40]=lambda*lambda;
y[41]=-2.*x0*x2*y[1]*y[5];
y[42]=y[12]+y[13]+y[14]+y[15]+y[16]+y[17]+y[41];
y[43]=y[1]*y[2];
y[44]=2.*x0*x1*y[1]*y[2];
y[45]=2.*x1*y[1]*y[2]*y[6];
y[46]=2.*x0*x2*y[1]*y[2];
y[47]=2.*x2*y[1]*y[2]*y[6];
y[48]=-(x0*x2*y[1]*y[5]);
y[49]=-(x2*y[1]*y[5]*y[6]);
y[50]=y[7]+y[43]+y[44]+y[45]+y[46]+y[47]+y[48]+y[49];
y[51]=-(lambda*MYI*y[11]*y[19]*y[50]);
y[52]=lrs[2];
y[53]=2.*y[1]*y[2]*y[3];
y[54]=4.*x1*x2*y[1]*y[2];
y[55]=2.*y[1]*y[2]*y[4];
y[56]=-2.*x1*x2*y[1]*y[5];
y[57]=y[53]+y[54]+y[55]+y[56];
y[58]=-(lambda*MYI*x0*y[9]*y[18]*y[57]);
y[59]=lambda*MYI*x0*y[18]*y[32];
y[60]=1.+y[33]+y[58]+y[59];
y[61]=y[7]+y[20];
y[62]=-(lambda*MYI*x1*y[11]*y[19]*y[61]);
y[63]=lambda*MYI*x1*y[19]*y[50];
y[64]=1.+y[51]+y[62]+y[63];
y[65]=-(x0*x1*y[1]*y[5]);
y[66]=-(x1*y[1]*y[5]*y[6]);
y[67]=y[7]+y[21]+y[43]+y[44]+y[45]+y[46]+y[47]+y[65]+y[66];
y[68]=-(lambda*MYI*x1*y[11]*y[19]*y[50]);
y[69]=x1+y[68];
y[70]=-(lambda*MYI*x0*y[9]*y[18]*y[32]);
y[71]=x0+y[70];
y[72]=-(lambda*MYI*x2*y[35]*y[52]*y[67]);
y[73]=x2+y[72];
y[74]=pow(y[69],2);
y[75]=pow(y[71],2);
y[76]=pow(y[73],2);
FOUT=(x0*x1*pow(bi,-2)*(1.+y[33])*(1.+y[51])*(lambda*MYI*x2*y[23]*y[35]*y[52\
]*(x0*x1*y[9]*y[11]*y[18]*y[19]*y[39]*y[40]*y[42]-lambda*MYI*x1*y[11]*y[19]\
*y[23]*y[60])-lambda*MYI*x2*y[35]*y[39]*y[52]*(-(x0*x1*y[9]*y[11]*y[18]*y[1\
9]*y[23]*y[40]*y[42])+lambda*MYI*x0*y[9]*y[18]*y[39]*y[64])+(x0*x1*pow(y[42\
],2)*y[9]*y[11]*y[18]*y[19]*y[40]+y[60]*y[64])*(1.-lambda*MYI*x2*y[35]*y[52\
]*y[61]+lambda*MYI*x2*y[52]*y[67]-lambda*MYI*y[35]*y[52]*y[67])))/((y[1]+y[\
1]*y[69]+y[1]*y[69]*y[71]+y[1]*y[73]+y[1]*y[71]*y[73])*(y[43]+y[1]*y[2]*y[6\
9]+2.*y[1]*y[2]*y[69]*y[71]+y[1]*y[2]*y[73]+2.*y[1]*y[2]*y[71]*y[73]-y[1]*y\
[5]*y[71]*y[73]+2.*y[1]*y[2]*y[69]*y[71]*y[73]-y[1]*y[5]*y[69]*y[71]*y[73]+\
y[1]*y[2]*y[71]*y[74]+2.*y[1]*y[2]*y[69]*y[73]*y[75]-y[1]*y[5]*y[69]*y[73]*\
y[75]+y[1]*y[2]*y[74]*y[75]+y[1]*y[2]*y[71]*y[76]+y[1]*y[2]*y[75]*y[76]));
return (FOUT);
}