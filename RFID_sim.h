#ifndef RFID
#define RFID

#include<math.h>
#include<stdlib.h>
#include<fstream.h>
#include<iostream.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#define pi acos(-1)

#define TAGNUM 100	//max number of tags
#define REANUM 12   //available readers
#define MAXX 60
#define MAXY 10     //simulation square
#define Gr 6.7		//reader antenna gain
#define Gt 3.7		//tag antenna gain
#define re 0.8		//reflection efficiency of tag
#define wl 0.328    //wavelength (frequency = 915 MHz)
#define Tr -80      //threshold value of reader-to-tag communication [dBm]
#define Tt -14      //threshold value of tag-to-reader communication [dBm]
#define LOSS 2      //attenuation factor including cable loss, polarization loss etc [dB]
#define nSpa 2.2    //varying for different environment 
int tagnumber;      //number of tags in this instance


struct TAG
{
	double posx;
	double posy;	//position of the tag
}; 
TAG tags[TAGNUM];

double tmp=wl*wl/4.0/4.0/pi/pi;

struct READER
{
	double posx;
	double posy;	//position of the reader
	double P1;      //the power of the reader send;
};
READER *readers=new READER[1];

int uncoverage;
int unuplink;
double inference;

inline double randval(double low,double high)			//产生随机数
{
	return (double(rand())/RAND_MAX)*(high-low)+low;
}

inline double dBm_to_mW(double dBm)
{
	return pow(10, double(dBm)/10);
}

inline double mW_to_dBm(double mW)
{
	//if(mW==0)return 0;
	return 10*log10(mW);
}

double cal_radius(double power)
{
	double m, radius;
	m=Tt-Gr-Gt-power+LOSS;
	radius=dBm_to_mW(m);
	radius=pow(double(tmp/radius),1.0/nSpa);
	return radius;	
}

int cal_cover(double x, double y, double r)
{
	int cover=0;
	double d_2;
	r*=r;
	for(int i=0; i<tagnumber; i++){
		d_2 = (tags[i].posx-x)*(tags[i].posx-x)+(tags[i].posy-y)*(tags[i].posy-y);
		if(d_2<=r)
			cover++;
	}
	return cover;	
}

void solution_to_reader(int readernum, int readeron[], double solu[])
{
	int k=0, j;
	delete readers;
	readers = new READER[readernum];
	for(int i=0; i<readernum; i++){
		while (readeron[k] == 0){
			k++;
		}
		j = 3*k;
		readers[i].posx = solu[j];
		readers[i].posy = solu[j+1];
		readers[i].P1 = solu[j+2];		
		k++;
	}	
}

void propagation(int readernum, int readeron[], double solu[])
{
	solution_to_reader(readernum, readeron, solu);
	
	int i, j;
	double d_2,d_n;
	double rP, recPower[TAGNUM], bP;
	double mind_2;
	
	inference = 0;
	uncoverage=0;
	unuplink=0;
	for(i=0; i<tagnumber; i++){
		mind_2 = 20000;
		recPower[i]=Tt;
		for(j=0; j<readernum; j++){
			d_2 = (tags[i].posx-readers[j].posx)*(tags[i].posx-readers[j].posx)+
				  (tags[i].posy-readers[j].posy)*(tags[i].posy-readers[j].posy);
			d_n = pow(d_2, double(nSpa)/2);
			rP = Gr+Gt+readers[j].P1+mW_to_dBm(tmp/d_n)-LOSS;
			if(rP>Tt){
				inference+=dBm_to_mW(rP);
				if(rP>recPower[i]){
					recPower[i]=rP;
				}
			}
			if(d_2<mind_2)mind_2=d_2;
		}
		if(recPower[i] != Tt)
			inference-=dBm_to_mW(recPower[i]);
		else
			uncoverage++;
		bP=Gr+Gt+recPower[i]+mW_to_dBm(re*re*tmp/mind_2);
		if(bP<Tr)unuplink++;
	}
}

void close_zero(int &readernum, int readeron[], double solu[])
{
	int i,j,k,index;
	int cover[REANUM];
	double rP,maxrP;
	double d_2, d_n;
	memset(cover,0,sizeof(cover));
	for(i=0; i<tagnumber; i++){
		maxrP=Tt;
		for(j=0; j<REANUM; j++){
			if(readeron[j] == 1){
				k=3*j;
				d_2 = (tags[i].posx-solu[k])*(tags[i].posx-solu[k])+
					  (tags[i].posy-solu[k+1])*(tags[i].posy-solu[k+1]);
				d_n = pow(d_2, double(nSpa)/2);
				rP = Gr+Gt+solu[k+2]+mW_to_dBm(tmp/d_n)-LOSS;
				if(rP > maxrP){
					index=j;
					maxrP=rP;
				}
			}
		}
		if(maxrP!=Tt)cover[index]++;	
	}
	
	for(i=0; i<REANUM; i++){
		if(readeron[i]==1 && cover[i]==0){
			readeron[i]=0;
			readernum--;
		}
	}	
}

int leastcover(int readeron[], double solu[], int lastclose)
{
	int i,j,k,index;
	int cover[REANUM];
	double rP,maxrP;
	double d_2,d_n;
	memset(cover,0,sizeof(cover));
	for(i=0; i<tagnumber; i++){
		maxrP=Tt;
		for(j=0; j<REANUM; j++){
			if(readeron[j] == 1){
				k=3*j;
				d_2 = (tags[i].posx-solu[k])*(tags[i].posx-solu[k])+
					  (tags[i].posy-solu[k+1])*(tags[i].posy-solu[k+1]);
				d_n = pow(d_2, double(nSpa)/2);
				rP = Gr+Gt+solu[k+2]+mW_to_dBm(tmp/d_n)-LOSS;
				if(rP > maxrP){
					index=j;
					maxrP=rP;
				}
			}
		}
		if(maxrP!=Tt)cover[index]++;	
	}

	int minc[REANUM], minnum=0, min_cover=tagnumber;
	for(i=0; i<REANUM; i++){
		if(readeron[i]==1 && i!=lastclose){
			if(cover[i]<min_cover){
				min_cover=cover[i];
				minnum=0;
				minc[minnum++]=i;
			}
			else if(cover[i]==min_cover){
				minc[minnum++]=i;
			}
		}
	}
	return minc[rand()%minnum];
}

#endif