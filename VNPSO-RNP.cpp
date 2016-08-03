
#include"RFID_sim.h"

#define PARSIZE 20
#define ROW 4
#define COL 5
#define MAXGENS 20000
#define FES 400000
#define dimensions 3*REANUM     
#define RUNTIME 10
#define RECOVERGEN 500
#define reducerangeP 2
#define deltaP 0.5
 
int READERON[REANUM];
int readerNum=REANUM;
int lastclose=-1;  
bool can_reduce;

struct FitnessType
{
	int uncover;
	int unuplink;
	int readers;
	double infere;
	double powersum;
};

struct ParticleType
{
	double velocity[dimensions];
	double position[dimensions];
	FitnessType fitness;	
	double pbest[dimensions];
	FitnessType pbestval;
 //	bool   inrange;
	int neighbor[4];
};
ParticleType particle[PARSIZE];
int bestpar;
FitnessType gbestval;

double lbound[dimensions];
double ubound[dimensions];
double vmax[dimensions];
const double c1=2.0;
const double c2=2.0;
double w=0.9; 

char filename_tag[100];
char filename_reader[100];
char filename_result[100];
 
int generation;
int fes;

#define betterthan <
#define worsethan >

void Initialize();
void Evaluate(int i);
void Update();
void Process();


void caleval(int num)		
{
	fes++;
	double powersum;
	propagation(readerNum, READERON, particle[num].position);
	particle[num].fitness.uncover=uncoverage;
	particle[num].fitness.unuplink=unuplink;
	particle[num].fitness.readers=readerNum;
	particle[num].fitness.infere=inference;
	powersum=0;
	for(int i=2; i<dimensions; i+=3){
		if(READERON[i/3] == 1)
			powersum+=dBm_to_mW(particle[num].position[i]);
	}
	particle[num].fitness.powersum=mW_to_dBm(powersum);
}

double cal_std(double array[],double avg,int times)
{
	double std=0;
	for(int i=0;i<times;i++)
		std+=(array[i]-avg)*(array[i]-avg);
		std=sqrt(std/(times-1));
	return std;
}

bool dominate(FitnessType a, FitnessType b)
{
	if(a.uncover betterthan b.uncover)
		return true;
	else if((a.uncover == b.uncover) && (a.unuplink betterthan b.unuplink))
		return true;
	else if((a.uncover == b.uncover) && (a.unuplink == b.unuplink) && (a.readers betterthan b.readers))
		return true;
	else if((a.uncover == b.uncover) && (a.unuplink == b.unuplink) && 
		(a.readers == b.readers) && (a.infere betterthan b.infere))
		return true;
	else if((a.uncover == b.uncover) && (a.unuplink == b.unuplink) && 
		(a.readers == b.readers) && (a.infere == b.infere)  && (a.powersum betterthan b.powersum))
		return true;
	return false;
}

void von_Neumann()
{
	int i;
	int row, col;
	for(i=0; i<PARSIZE; i++){
		row=i/ROW;
		col=i%ROW;
		particle[i].neighbor[0]=COL*row+(col+COL-1)%COL;	//left
		particle[i].neighbor[1]=COL*row+(col+1)%COL;		//right
		particle[i].neighbor[2]=COL*((row+ROW-1)%ROW)+col;	//up
		particle[i].neighbor[3]=COL*((row+1)%ROW)+col;      //down
	}	
}

int bestneigbor(int num)
{
	int index, bestn;
	bestn=index=particle[num].neighbor[0];
	FitnessType ff=particle[index].pbestval;
	for(int i=1; i<4; i++){
		index=particle[num].neighbor[i];
		if(dominate(particle[index].pbestval, ff))
			bestn=index;
	}
	return bestn;	
}

void read_tagfile()
{
	FILE *tagDistribution = fopen(filename_tag,"r");
	for(int i=0; i<TAGNUM; i++){
		fscanf(tagDistribution,"%lf",&tags[i].posx);
		fscanf(tagDistribution,"%lf",&tags[i].posy);
	}
	fclose(tagDistribution);	
}


void Initialize()
{
	int i,j;

	for(i=0;i<dimensions;i++){
		if(i%3==0){
			lbound[i]=0;ubound[i]=MAXX;
		}
		else if(i%3==1){
			lbound[i]=0;ubound[i]=MAXY;
		}
		else{
			lbound[i]=20;ubound[i]=33;  //0.1~2 W
		}
		vmax[i] = 0.2*(ubound[i]-lbound[i]);
	}
	for(i=0; i<REANUM; i++) READERON[i]=1;

	for(i=0;i<PARSIZE;i++){
		for(j=0;j<dimensions;j++){
			particle[i].position[j]=randval(lbound[j],ubound[j]);
			particle[i].velocity[j]=randval(-vmax[j],vmax[j]);
			particle[i].pbest[j]=particle[i].position[j];			
		}
		caleval(i);
		particle[i].pbestval=particle[i].fitness;
// 		particle[i].inrange=true;
	}
	gbestval=particle[0].pbestval;
	bestpar=0;
	for(i=1;i<PARSIZE;i++){
		if(dominate(particle[i].pbestval, gbestval)){
			gbestval = particle[i].pbestval;
			bestpar = i;
		}
	}
	fes=0;
	generation=0;	
	can_reduce=false;
	readerNum=REANUM;
	von_Neumann();
}

void Evaluate(int i)
{
	int j;	
	
	caleval(i);
	if(particle[i].fitness.uncover==0 /*&& reduce_or_recover ==0*/)
		can_reduce=true;
		//no_need_to_recover=1;

	if(dominate(particle[i].fitness, particle[i].pbestval)){
		for(j=0;j<dimensions;j++)
			particle[i].pbest[j]=particle[i].position[j];
		particle[i].pbestval=particle[i].fitness;
		if(dominate(particle[i].pbestval, gbestval)){
			gbestval=particle[i].pbestval;
			bestpar=i;	
		}
	}
}

void Update()
{
	int i,j;
	double rand1,rand2;
	int bestn;
	for(i=0;i<PARSIZE;i++){
		w= 0.9-0.5*fes/(FES-1);
// 		particle[i].inrange=true;
		bestn=bestneigbor(i);
		for(j=0;j<dimensions;j++){		
			if(READERON[j/3] == 1){
				rand1 = randval(0,1);
				rand2 = randval(0,1);
				particle[i].velocity[j]=particle[i].velocity[j]*w+c1*rand1*(particle[i].pbest[j]
					-particle[i].position[j])+c2*rand2*(particle[bestn].pbest[j]-particle[i].position[j]);
				if(particle[i].velocity[j]<-vmax[j])particle[i].velocity[j]=-vmax[j];
				else if(particle[i].velocity[j]>vmax[j])particle[i].velocity[j]=vmax[j];			
				particle[i].position[j]+=particle[i].velocity[j];			          
				//if(particle[i].position[j]<lbound[j]||particle[i].position[j]>ubound[j])particle[i].inrange=false;
				if(particle[i].position[j]<lbound[j])particle[i].position[j]=lbound[j];
				else if(particle[i].position[j]>ubound[j])particle[i].position[j]=ubound[j];
			}
		}
		//if(particle[i].inrange)
			Evaluate(i);
	}
}

void reduceReader()
{
	lastclose = leastcover(READERON, particle[bestpar].pbest, lastclose);
	readerNum--;
	READERON[lastclose]=0;	
}

void recoverReader()
{
	readerNum++;
	READERON[lastclose]=1;	
}

void reducePower(int num)
{
	int r;
	do{
		r=rand()%REANUM;
	}while(READERON[r]==0);
	r=r*3+2;
	double p=randval(0,reducerangeP);
	particle[num].position[r]-=p;
	if(particle[num].position[r]<lbound[r])particle[num].position[r]=lbound[r];
}

void final_reducePower()
{
	int i, k;
	double tt=lbound[2]+deltaP;
	for(i=0; i<REANUM; i++){
		if(READERON[i]){
			k=3*i+2;
			do{
				particle[bestpar].pbest[k]-=deltaP;
				propagation(readerNum, READERON, particle[bestpar].pbest);
			} while(uncoverage==0 && unuplink==0 && particle[bestpar].pbest[k]>tt);
			if(uncoverage!=0 || unuplink!=0) particle[bestpar].pbest[k]+=deltaP;
		}
	}
	double powersum;
	propagation(readerNum, READERON, particle[bestpar].pbest);
	particle[bestpar].pbestval.uncover=uncoverage;
	particle[bestpar].pbestval.unuplink=unuplink;
	particle[bestpar].pbestval.readers=readerNum;
	particle[bestpar].pbestval.infere=inference;
	powersum=0;
	for(i=2; i<dimensions; i+=3){
		if(READERON[i/3] == 1)
			powersum+=dBm_to_mW(particle[bestpar].pbest[i]);
	}
	particle[bestpar].pbestval.powersum=mW_to_dBm(powersum);
}

void Mutation(int num)
{
	int r;
//	do{
		r=rand()%dimensions;
//	}while(READERON[r/3]==0);
	double maxdelta=0.2*(ubound[r]-lbound[r]);
	double x=randval(-maxdelta,maxdelta);
	particle[num].position[r]+=x;
	if(particle[num].position[r]<lbound[r])particle[num].position[r]=lbound[r];
	if(particle[num].position[r]>ubound[r])particle[num].position[r]=ubound[r];
	
}

void Process()
{
	int num;
	int recoverGen=0;
//	FILE *pro=fopen("process.txt","w");
//	FILE *pro2=fopen("process2.txt","w");
	Initialize();	
	while(generation<MAXGENS){
		Update();
//		num=rand()%PARSIZE;
//		reducePower(num);
		num=rand()%PARSIZE;
		Mutation(num);

		if(can_reduce && recoverGen==0){
			recoverGen=RECOVERGEN;
			reduceReader();
			can_reduce=false;				
		}
		if(recoverGen>1)
			recoverGen--;
		if(recoverGen==1){
			recoverGen--;
			if(!can_reduce)
				recoverReader();				
		}
		generation++;
//		fprintf(pro,"%d\t%d\t%d\t%d\t%f\t%f\n",generation,particle[bestpar].pbestval.readers, particle[bestpar].pbestval.uncover, 
//				particle[bestpar].pbestval.unuplink,particle[bestpar].pbestval.infere,particle[bestpar].pbestval.powersum);
//		fprintf(pro2,"%d\t%d\t%d\t%d\t%f\n",generation,particle[bestpar].fitness.readers, particle[bestpar].fitness.uncover, 
//				particle[bestpar].fitness.unuplink,particle[bestpar].fitness.infere);
	}
	if(recoverGen!=0 && !can_reduce){
		recoverReader();
		//can_reduce=true;
	}
	final_reducePower();
//	fclose(pro);
//	fclose(pro2);
}


void main()
{
	int k;
	double x,y,Power,Radius;
	int cover;
	FILE *readerDistribution, *finalResult; 
	FILE *totalResult=fopen("totalresult.txt","a");
	fprintf(totalResult,"instance\tavg_reaNum\tavg_unCov\tavg_unUpl\tavg_inf\tavg_powSum\t");
	fprintf(totalResult,"bst_reaNum\tbst_unCov\tbst_unUpl\tbst_inf\tbst_powSum\n");
	srand((unsigned)time(NULL));
	for(int instance=6; instance<=6; instance++){
		switch(instance){
		case 0:
			tagnumber=30;
			strcpy(filename_tag, "C_50x50_30.txt");
			strcpy(filename_reader, "reader_distribution_C_50x50_30.txt");
			strcpy(filename_result,"final_result_C_50x50_30.txt");
			break;
		case 1:
			tagnumber=50;
			strcpy(filename_tag, "C_50x50_50.txt");
			strcpy(filename_reader, "reader_distribution_C_50x50_50.txt");
			strcpy(filename_result,"final_result_C_50x50_50.txt");
			break;		
		case 2:
			tagnumber=100;
			strcpy(filename_tag, "C_50x50_100.txt");
			strcpy(filename_reader, "reader_distribution_C_50x50_100.txt");
			strcpy(filename_result,"final_result_C_50x50_100.txt");
			break;
		case 3:
			tagnumber=30;
			strcpy(filename_tag, "R_50x50_30.txt");
			strcpy(filename_reader, "reader_distribution_R_50x50_30.txt");
			strcpy(filename_result,"final_result_R_50x50_30.txt");
			break;
		case 4:
			tagnumber=50;
			strcpy(filename_tag, "R_50x50_50.txt");
			strcpy(filename_reader, "reader_distribution_R_50x50_50.txt");
			strcpy(filename_result,"final_result_R_50x50_50.txt");
			break;
		case 5:
			tagnumber=100;
			strcpy(filename_tag, "R_50x50_100.txt");
			strcpy(filename_reader, "reader_distribution_R_50x50_100.txt");
			strcpy(filename_result,"final_result_R_50x50_100.txt");
			break;
		case 6:
			tagnumber=18;
			strcpy(filename_tag, "realapp.txt");
			strcpy(filename_reader, "reader_distribution_realapp.txt");
			strcpy(filename_result,"final_result_realapp.txt");
			break;
		default:
			printf("Oops, no such instance");
		}

		readerDistribution = fopen(filename_reader,"w");
		finalResult = fopen(filename_result, "w");

		read_tagfile();
		printf("Program is running......\n");
		printf("runTime\treaderNum\tunCoverage\tunUplink\tInference\tpowerSum\n");
		fprintf(finalResult, "runTime\treaderNum\tunCoverage\tunUplink\tInference\tpowerSum\n");

		double avg_readerNum=0, avg_unCoverage=0, avg_unUplink=0, avg_Inference=0, avg_powSum=0;
		FitnessType bestSolu;
		bestSolu.uncover=TAGNUM;bestSolu.unuplink=TAGNUM;bestSolu.readers=REANUM;
		bestSolu.infere=100;bestSolu.powersum=10000;
		for(int run=1; run<=RUNTIME; run++){
			Process();			
			fprintf(readerDistribution,"%d\n**********\n",run);
			for(int i=0; i<REANUM; i++){
				if(READERON[i]){
					k=3*i;
					x=particle[bestpar].pbest[k];
					y=particle[bestpar].pbest[k+1];
					Power=particle[bestpar].pbest[k+2];
					Radius=cal_radius(Power);
					cover=cal_cover(x,y,Radius);
				/*	if(cover==0){
						readerNum--;
						particle[bestpar].pbestval.powersum-=Power;
					}
					else*/
						fprintf(readerDistribution,"%f\t%f\t%f\t%f\t%d\n",x,y,Power,Radius,cover);
				}
			}
			fprintf(readerDistribution,"\n");
			printf("%d\t%d\t%d\t%d\t%f\t%f\n",run, readerNum, particle[bestpar].pbestval.uncover,
					particle[bestpar].pbestval.unuplink,particle[bestpar].pbestval.infere,particle[bestpar].pbestval.powersum);		
			fprintf(finalResult, "%d\t%d\t%d\t%d\t%f\t%f\n",run, readerNum, particle[bestpar].pbestval.uncover,
					particle[bestpar].pbestval.unuplink,particle[bestpar].pbestval.infere, particle[bestpar].pbestval.powersum);
			particle[bestpar].pbestval.readers=readerNum;
			avg_readerNum+=readerNum;
			avg_unCoverage+=particle[bestpar].pbestval.uncover;
			avg_unUplink+=particle[bestpar].pbestval.unuplink;
			avg_Inference+=particle[bestpar].pbestval.infere;
			avg_powSum+=particle[bestpar].pbestval.powersum;	
			if(dominate(particle[bestpar].pbestval,bestSolu))
				bestSolu=particle[bestpar].pbestval;
		}
		avg_readerNum/=(double)RUNTIME;
		avg_unCoverage/=(double)RUNTIME;
		avg_unUplink/=(double)RUNTIME;
		avg_Inference/=(double)RUNTIME;
		avg_powSum/=(double)RUNTIME;
		printf("AVG\t%f\t%f\t%f\t%f\t%f\n",avg_readerNum, avg_unCoverage, avg_unUplink, avg_Inference,avg_powSum);		
//		fprintf(finalResult, "AVG\t%f\t%f\t%f\t%f\t%f\n",avg_readerNum, avg_unCoverage, avg_unUplink, avg_Inference, avg_powSum);
		printf("BST\t%d\t%d\t%d\t%f\t%f\n",bestSolu.readers, bestSolu.uncover, bestSolu.unuplink, bestSolu.infere, bestSolu.powersum);		
//		fprintf(finalResult, "BST\t%d\t%d\t%d\t%f\t%f\n",bestSolu.readers, bestSolu.uncover, bestSolu.unuplink, bestSolu.infere, bestSolu.powersum);
		printf("Successful!\n");
		fclose(readerDistribution);
		fclose(finalResult);

		fprintf(totalResult,"%d\t%f\t%f\t%f\t%f\t%f\t",instance,avg_readerNum, avg_unCoverage, avg_unUplink, 
			avg_Inference,avg_powSum);
		fprintf(totalResult,"%d\t%d\t%d\t%f\t%f\n",bestSolu.readers, 
			bestSolu.uncover, bestSolu.unuplink, bestSolu.infere, bestSolu.powersum);
	}
	fclose(totalResult);
}	