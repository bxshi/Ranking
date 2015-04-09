#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<time.h>
#include<algorithm>
#include<iostream>
#include<math.h>
#include<string>
#include<string.h>

using namespace std;


vector< vector<int> > edge;
vector< vector<int> > invedge;
double *maph,*mapa,*olda;
int newid;

double delta(int u, int v)
{
  if (u==v)
    return 1;
  else
    return 0;
}

double abs1(double x)
{
  if (x>=0)
    return x;
  else
    return -x;
}

void calcSalsa(int start,int times, double df){
	double tempv;
	int i,j,ttx,ttv,k;

  for (i=0; i<newid; i++) {
    olda[i]=1.0/double(newid);
    mapa[i]=1.0/double(newid);
    maph[i]=0;
  }

  for (i=0; i<2*times; i++) {
    if (i%2==0) {
      for (j=0; j<newid; j++) {
        tempv=0;
        for (k=1; k<edge[j].size(); k++) {
          ttx=edge[j][k];
          tempv=tempv+mapa[ttx]/double(invedge[ttx].size()-1);
        }
        maph[j]=df*delta(j,start)+(1-df)*tempv;
      }
    } 
		else {
      for (j=0; j<newid; j++) {
        tempv=0;
        for (k=1; k<invedge[j].size(); k++) {
          ttv=invedge[j][k];
          tempv=tempv+maph[ttv]/double(edge[ttv].size()-1);
        }
        mapa[j]=tempv;
      }
      tempv=0;
      for (j=0; j<newid; j++) {
        tempv=tempv+abs1(mapa[j]-olda[j]);
        olda[j]=mapa[j];
      }
      if (tempv<1e-15)
        break;
      tempv=0;
      for (j=0; j<newid; j++)
        tempv=tempv+mapa[j];
    }
  }
}

int main(int argc, char **argv)
{
  FILE *fin,*fout;
  char line[300];
  char file[100],ofile[100],qfile[100];
  int times,len,topK,tempid,i,j,k,a,b,numedge=0,maxvv=0,start,sid;
  int *trans;
  int *invtrans;
  int *dist,*outputid;
  double tempv,df;
  srand(time(NULL));
  if (argc!=8) {
    printf("Usage: Edge_file Output_file StartId T topK Dumping Query_file\n");
    return 2;
  }

  strcpy(file,argv[1]);
  strcpy(ofile,argv[2]);
  sscanf(argv[3],"%d",&start);
  sscanf(argv[4],"%d",&times);
  sscanf(argv[5],"%d",&topK);
  sscanf(argv[6],"%lf",&df);
	strcpy(qfile,argv[7]);

  fin=fopen(file,"r");
  if (!fin) {
    printf("File doesn't exist\n");
    return 1;
  }
  while (fgets(line,100,fin)!=NULL) {
    sscanf(line,"%d %d",&a,&b);
    numedge++;
    if (a>maxvv)
      maxvv=a;
    if (b>maxvv)
      maxvv=b;
    //printf("%d %d\n",a,b);
  }
  maxvv++;
  trans=(int *) malloc(maxvv*sizeof(int));
  for (i=0; i<maxvv; i++)
    trans[i]=-1;
  fclose(fin);
  a=0;
  b=0;
  fin=fopen(file,"r");
	newid=0;
  while (fgets(line,100,fin)!=NULL) {
    sscanf(line,"%d %d",&a,&b);
    if (trans[a]==-1) {
      trans[a]=newid;
      newid++;
    }
    if (trans[b]==-1) {
      trans[b]=newid;
      newid++;
    }
  }

  invtrans=(int *)malloc(newid*sizeof(int));
  mapa=(double *)malloc(newid*sizeof(double));
  maph=(double *)malloc(newid*sizeof(double));
  olda=(double *)malloc(newid*sizeof(double));
  for (i=0; i<maxvv; i++)
    if (trans[i]!=-1)
      invtrans[trans[i]]=i;

  for (i=0; i<newid; i++) {
    vector<int> temprow (1,-1);
    edge.push_back(temprow);
    invedge.push_back(temprow);
  }

  fclose(fin);
  fin=fopen(file,"r");
  while (fgets(line,100,fin)!=NULL) {
    sscanf(line,"%d %d",&a,&b);
    if (find(edge[trans[a]].begin(),edge[trans[a]].end(),trans[b])==edge[trans[a]].end()) {
      edge[trans[a]].push_back(trans[b]);
      invedge[trans[b]].push_back(trans[a]);
    }
  }
  fclose(fin);
  
  int numquery,qstart,qend,arank,hrank,starank,sthrank;
  double ascore,hscore,stascore,sthscore;
  fin=fopen(qfile,"r");
  numquery=20;
  fout=fopen(ofile,"w+");
  bool edge_exist;
  for (i=0; i<numquery; i++) {
    fscanf(fin,"%d %d",&qstart,&qend);    
    if (find(edge[trans[qstart]].begin(),edge[trans[qstart]].end(),trans[qend])==edge[trans[qstart]].end())
        edge_exist=0;
    else
        edge_exist=1;
    if (!edge_exist){
       edge[trans[qstart]].push_back(trans[qend]);
       invedge[trans[qend]].push_back(trans[qstart]);
       printf("Insert Edge %d -> %d\n",qstart,qend); 
    }
		calcSalsa(trans[qstart],times,df);
    ascore=mapa[trans[qend]];
    arank=0;
    for (j=0; j<newid; j++)
      if (mapa[j]>ascore)
        arank++;

    hscore=maph[trans[qend]];
    hrank=0;
    for (j=0; j<newid; j++)
      if (maph[j]>hscore)
        hrank++;

    stascore=mapa[trans[qstart]];
    starank=0;
    for (j=0; j<newid; j++)
      if (mapa[j]>stascore)
        starank++;

    sthscore=maph[trans[qstart]];
    sthrank=0;
    for (j=0; j<newid; j++)
      if (maph[j]>sthscore)
        sthrank++;

    if (!edge_exist){
        printf("Delete Edge %d -> %d\n",invtrans[invedge[trans[qend]].back()],invtrans[edge[trans[qstart]].back()]);
        edge[trans[qstart]].pop_back();
        invedge[trans[qend]].pop_back();
    }
    printf("%d %d %f %d %f %d %f %d %f %d\n",qstart,qend,stascore,starank,sthscore,sthrank,ascore,arank,hscore,hrank);
    fprintf(fout,"%d %d %f %d %f %d %f %d %f %d\n",qstart,qend,stascore,starank,sthscore,sthrank,ascore,arank,hscore,hrank);
  }
  fclose(fout);
  fclose(fin);
}

