/*******************************
*    CalcHEP  3.8.7*
*******************************/
#include"num_out.h"
#include"num_in.h"
extern FNN F57_a32a32_a10a10p1;
static void C57(REAL*V,REAL * C)
{
REAL S[11];                                                                 
     
S[0]=V[8]*V[8];
S[1]=V[30]*V[30];
S[2]=V[29]*V[29];
S[3]=V[28]*V[28];
S[4]=V[8]*V[8]*V[8]*V[8];
C[20]=+S[0]*(V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*
 V[12])+2*V[30]*V[27]*V[11])+V[29]*V[27]*V[12]*V[11]-4*V[30]*V[28]*S[0])+
 S[1]*(V[12]*(2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*
 V[30]*V[29]*V[3])))+S[0]*(V[27]*(V[12]*(V[28]*(-S[1]-S[2]))-V[30]*V[29]*
 V[27]*V[11])-V[30]*V[29]*S[3]*V[11]))+2*V[30]*V[29]*V[28]*V[27]*S[4]);
S[5]=V[27]*V[27];
C[19]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*V[12])+2*
 V[30]*V[27]*V[11])+V[29]*V[27]*V[12]*V[11]-4*V[30]*V[28]*S[0])+S[1]*(V[12]*
 (2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*V[30]*V[29]*
 V[3])))+S[0]*(V[11]*(V[29]*(V[30]*(2*(-S[3]-S[5]))))))+2*V[30]*V[29]*V[28]*
 V[27]*S[4];
C[18]=+V[3]*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))-V[30]*V[29]*V[27]*V[11])-
 V[30]*V[29]*S[3]*V[11])-2*V[30]*V[29]*V[28]*V[27]*S[0];
C[17]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+2*V[29]*V[12])+2*
 V[30]*V[27]*V[11])+V[29]*V[27]*V[12]*V[11]-2*V[30]*V[28]*S[0])+S[1]*(V[12]*
 (2*V[28]*V[3]+V[27]*V[11])))+S[3]*(V[11]*(V[12]*(S[1]+S[2])+2*V[30]*V[29]*
 V[3])))+S[0]*(V[11]*(V[29]*(V[30]*(-S[3]-S[5])))))+2*V[30]*V[29]*V[28]*
 V[27]*S[4];
C[16]=+V[3]*(V[27]*(V[28]*(V[12]*(S[1]+S[2])-2*V[30]*V[29]*V[3])-2*V[30]*
 V[29]*V[27]*V[11])-2*V[30]*V[29]*S[3]*V[11]);
C[15]=+S[0]*(V[27]*(V[28]*(V[3]*(V[12]*(S[1]+S[2])+2*V[30]*V[29]*V[3])-2*
 V[30]*V[29]*S[0])));
C[14]=+V[3]*(V[3]*(V[11]*(V[12]*(S[5]*(S[1]+S[2])+S[3]*(S[1]+S[2]))+V[3]*(
 V[29]*(V[30]*(2*(S[3]+S[5])))))+V[3]*(V[27]*(V[28]*(V[12]*(2*(S[1]+S[2]))+
 4*V[30]*V[29]*V[3]))))+S[0]*(V[27]*(V[12]*(V[28]*(S[1]+S[2]))-V[30]*V[29]*
 V[27]*V[11])-V[30]*V[29]*S[3]*V[11]));
C[13]=+V[3]*(V[27]*(V[28]*(V[12]*(3*(S[1]+S[2]))+4*V[30]*V[29]*V[3])-V[30]*
 V[29]*V[27]*V[11])-V[30]*V[29]*S[3]*V[11])-4*V[30]*V[29]*V[28]*V[27]*S[0];
C[12]=+V[3]*(V[27]*(V[28]*(V[12]*(2*(S[1]+S[2]))+4*V[30]*V[29]*V[3])+V[30]*
 V[29]*V[27]*V[11])+V[30]*V[29]*S[3]*V[11])+V[11]*(V[12]*(S[1]*S[3]+S[2]*
 S[5]));
C[11]=+V[3]*(V[29]*(V[30]*(V[11]*(S[3]+S[5])+2*V[28]*V[27]*V[3])));
C[10]=+V[27]*(V[28]*(V[3]*(V[12]*(2*(S[1]+S[2]))+4*V[30]*V[29]*V[3])-2*
 V[30]*V[29]*S[0]));
C[9]=+V[3]*(V[3]*(V[27]*(V[29]*(V[3]*(V[28]*(-4*V[30]*V[3]-2*V[29]*V[12])-2*
 V[30]*V[27]*V[11])-4*V[30]*V[28]*S[0]-V[29]*V[27]*V[12]*V[11])+S[1]*(V[12]*
 (-2*V[28]*V[3]-V[27]*V[11])))+S[3]*(V[11]*(V[12]*(-S[1]-S[2])-2*V[30]*
 V[29]*V[3])))+S[0]*(V[27]*(V[12]*(V[28]*(-S[1]-S[2]))+V[30]*V[29]*V[27]*
 V[11])+V[30]*V[29]*S[3]*V[11]))+S[0]*(V[11]*(V[12]*(2*(S[1]*S[3]+S[2]*
 S[5])))+4*V[30]*V[29]*V[28]*V[27]*S[0]);
C[8]=+V[27]*(V[29]*(V[3]*(V[30]*V[27]*V[11]+V[29]*V[28]*V[12])+4*V[30]*
 V[28]*S[0]+2*V[29]*V[27]*V[12]*V[11])+S[1]*V[28]*V[12]*V[3])+S[3]*(V[11]*(
 V[30]*(2*V[30]*V[12]+V[29]*V[3])));
S[6]=V[3]*V[3];
C[7]=+V[12]*(V[27]*(S[2]*(V[28]*V[3]+2*V[27]*V[11])+S[1]*V[28]*V[3])+2*S[1]*
 S[3]*V[11])-2*V[30]*V[29]*V[28]*V[27]*S[6];
C[6]=+V[27]*(V[29]*(V[3]*(V[28]*(4*V[30]*V[3]+V[29]*V[12])+2*V[30]*V[27]*
 V[11])+4*V[30]*V[28]*S[0]+V[29]*V[27]*V[12]*V[11])+S[1]*V[28]*V[12]*V[3])+
 S[3]*(V[11]*(V[30]*(V[30]*V[12]+2*V[29]*V[3])));
C[5]=+V[3]*(V[27]*(V[28]*(V[12]*(2*(S[1]+S[2]))+8*V[30]*V[29]*V[3])+2*V[30]*
 V[29]*V[27]*V[11])+2*V[30]*V[29]*S[3]*V[11]);
C[4]=+V[29]*(V[30]*(V[3]*(V[11]*(2*(S[3]+S[5]))+4*V[28]*V[27]*V[3])+6*V[28]*
 V[27]*S[0]));
C[3]=+2*V[30]*V[29]*V[28]*V[27];
C[2]=+4*V[30]*V[29]*V[28]*V[27];
S[7]=V[7]*V[7]*V[7]*V[7];
S[8]=V[6]*V[6]*V[6]*V[6];
S[9]=V[1]*V[1]*V[1]*V[1];
C[1]=+S[7]*S[8]*S[9];
S[10]=V[2]*V[2]*V[2]*V[2]*V[2]*V[2];
C[0]=+32*S[10];
}
REAL F57_a32a32_a10a10p1(double GG, REAL*V, REAL*DP,REAL*Q0,COMPLEX*Q1,REAL*
 Q2,REAL*cb_coeff,int gsw,int gtw)
{
REAL N,D,R; COMPLEX Prop;
static REAL C[21];                                                          
     
if(!DP){C57(V,C); return 0;} 
  REAL N_p1p2_=1/DP[0];
N=-C[0];
D=+C[1];
R=+DP[0]*(C[3]*(DP[0]*(DP[0]-DP[1]-DP[2])+DP[3]*(DP[1]+DP[3])+DP[4]*(-DP[2]-
 DP[4]))+DP[3]*(C[2]*(DP[2]-DP[0]+DP[4])+C[13])+C[16]*DP[1]-C[19]-C[18]*
 DP[0]-C[11]*DP[2]+C[8]*DP[4])+DP[3]*(DP[4]*(C[3]*(DP[2]-DP[1])+C[2]*(DP[4]-
 DP[3])+C[5])+DP[2]*(C[11]-C[3]*DP[3])+C[14]-C[12]*DP[1]-C[10]*DP[3])+DP[4]*
 (DP[1]*(C[3]*DP[4]-C[7])+C[9]-C[6]*DP[2]-C[4]*DP[4])+C[17]*DP[1]-C[20]-
 C[15]*DP[2];
R*=(N/D);
Prop=1*(gtw ? creal(Q1[1]):Q1[1])*(gsw ? creal(Q1[4]):Q1[4])*(gtw ? creal(
 Q1[20]):conj(Q1[20]))*(gtw ? creal(Q1[10]):conj(Q1[10]));
R*=creal(Prop);
 return R;
}