#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;

} Point;

//Leer las dimensiones de la imagen a x b pixeles  
void readImageSize(FILE *ifp,int* K,int* a,int* b){
	fscanf(ifp,"%d\n",K);
	printf("%d\n",*K);

	fscanf(ifp,"%d\n",a);
	printf("%d\n",*a);

	fscanf(ifp,"%d\n",b);
	printf("%d\n",*b);
}

//Leer archivo ifp y lo almacena en el struct
void readPoints(FILE* ifp,Point *points, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
	{
	fscanf(ifp,"%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
	//printf("%lf,%lf,%lf,%lf,%lf\n", points[i]._r, points[i]._g, points[i]._b, points[i]._m, points[i]._n);
	}
}

//Inicializar puntos aleatorios como los medios (k numero de clusters) 
void initialize(Point* mean,int K, int num_points, Point* points)
{
	int i, a, p=2;
	srand(time(NULL));
	for(i=0;i<K;i++){
		a = num_points/p;
		//printf("\n num_points: %d\n", num_points/p);
		mean[i]._r = points[a]._r;
		mean[i]._g = points[a]._g;
		mean[i]._b = points[a]._b;
		mean[i]._m = points[a]._m;
		mean[i]._n = points[a]._n;
		/*mean[i]._r=((double)(rand()%1000))/1000;
		mean[i]._g=((double)(2*rand()%1000))/1000;
		mean[i]._b=((double)(3*rand()%1000))/1000;
		mean[i]._m=((double)(4*rand()%1000))/1000;
		mean[i]._n=((double)(5*rand()%1000))/1000;*/
		//printf("%lf,%lf,%lf,%lf,%lf\n",mean[i]._r,mean[i]._g,mean[i]._b,mean[i]._m,mean[i]._n);
		p++;
		/*mean[i]._r=((double)(rand()%1000))/1000;
		mean[i]._g=((double)(2*rand()%1000))/1000;
		mean[i]._b=((double)(3*rand()%1000))/1000;
		mean[i]._m=((double)(4*rand()%1000))/1000;
		mean[i]._n=((double)(5*rand()%1000))/1000;*/
	}
}

//Todos los puntos sin clusters
int IntClusterMem(int *cluster, int num_points)
{
	int i;
	for(i=0;i<num_points;i++){
		cluster[i]=-1;
	}
}

//Distancia
double calculateDistance(Point point1,Point point2){
	return sqrt((pow((point1._r-point2._r),2)+pow((point1._g-point2._g),2)+pow((point1._b-point2._b),2)+pow((point1._m-point2._m),2)+pow((point1._n-point2._n),2)));
}

//Para calcular a que cluster pertenece el punto. (k numero de clusters)
int pointsCluster(Point point,Point* mean,int K)
{
	int parent=0;
	double dist=0;
	double minDist=calculateDistance(point,mean[0]);
	int i;
	for(i=1;i<K;i++)
	{
		dist=calculateDistance(point,mean[i]);
		if(minDist>=dist){
			parent=i;
			minDist=dist;
		}
	}
	return parent;
}

//Calcular nueva media
void calcNewMean(Point* points,int* cluster,Point* mean,int K,int num_points)
{
	//Nuevo vector de puntos de las nuevas medias
	Point* newMean=malloc(sizeof(Point)*K);
	
	int* members=malloc(sizeof(int)*K);
	int i;
	//Recorremos cada cluster
	for(i=0;i<K;i++)
	{
		members[i]=0;
		newMean[i]._r=0;
		newMean[i]._g=0;
		newMean[i]._b=0;
		newMean[i]._m=0;
		newMean[i]._n=0;
	}
	
	for(i=0;i<num_points;i++)
	{
		members[cluster[i]]++;
		newMean[cluster[i]]._r+=points[i]._r;
		newMean[cluster[i]]._g+=points[i]._g;
		newMean[cluster[i]]._b+=points[i]._b;
		newMean[cluster[i]]._m+=points[i]._m;
		newMean[cluster[i]]._n+=points[i]._n;
	}
	for(i=0;i<K;i++){
		if(members[i]!=0.0){
			newMean[i]._r/=members[i];
			newMean[i]._g/=members[i];
			newMean[i]._b/=members[i];
			newMean[i]._m/=members[i];
			newMean[i]._n/=members[i];
		}
		else{
		newMean[i]._r=0;
		newMean[i]._g=0;
		newMean[i]._b=0;
		newMean[i]._m=0;
		newMean[i]._n=0;
		}
	}
	//Actualizamos la media de cada cluster
	for(i=0;i<K;i++){
		mean[i]._r=newMean[i]._r;
		mean[i]._g=newMean[i]._g;
		mean[i]._b=newMean[i]._b;
		mean[i]._m=newMean[i]._m;
		mean[i]._n=newMean[i]._n;
	}
}

//Comprobamos la convergencia
//Comprueba que cada clúster de puntos permanece igual
int chkConvrg(int *before_clusters,int *after_cluster,int num_points, float tol){
	int i;
	tol = num_points*tol;
	for(i=0;i<num_points;i++)
		if(abs(before_clusters[i]-after_cluster[i])>tol)
			return -1;
	return 0;
}

int main(int argc, char* argv[]){
	int K;
	int num_points;
	int i;
	int job_done=0;
	int x,y;
	float tol;
	double tic, toc;
	int iter=0;
	double tspantemp = 0.0, tspan = 0.0;

	Point* mean;
	Point* points;
	Point* get_points;
	int * formed_clusters;
	int * before_clusters;
	int * after_cluster;
	int thread_no;
	//clock_t tic, toc;

	tol = atof(argv[3]);
	printf("Tolerancia:  %f \n",tol);
	thread_no = atoi(argv[4]);
	//printf("Max threads: %d\n",thread_no);
	//Leer archivo

	FILE *ifp;
	ifp=fopen(argv[1],"r");
	readImageSize(ifp,&K,&x,&y);
	num_points = x*y;
	points=(Point*)malloc(sizeof(Point)*num_points);
	readPoints(ifp,points,num_points);
	fclose(ifp);

	before_clusters=(int*)malloc(sizeof(int)*num_points);
	after_cluster=(int*)malloc(sizeof(int)*num_points);
	mean=malloc(sizeof(Point)*K);

	//Inicializando a valores por defecto
	initialize(mean,K,num_points,points);
	IntClusterMem(before_clusters,num_points);
	IntClusterMem(after_cluster,num_points);
	int iter_max=0;
	char *fcentroid;
	if(K==10){
		iter_max=4;
		fcentroid="Centroid_10.txt";
	}
	else if(K==30){
		iter_max=14;
		fcentroid="Centroid_30.txt";
	}
	else if(K==50){
		iter_max=63;
		fcentroid="Centroid_50.txt";
	}
	else{
		iter_max=87;
		fcentroid="Centroid_100.txt";
	}
	//tic = omp_get_wtime();
	while(iter<iter_max)
	{
		
		//printf("at 211: tic is %.3f\n\n", myclocktemp);
		tic = omp_get_wtime();
		omp_set_num_threads(thread_no);
		#pragma omp parallel default(shared)
		//tic = clock();

		#pragma omp for

		for(i=0;i<num_points;i++)
		{

			after_cluster[i]=pointsCluster(points[i],mean,K);
		}

		calcNewMean(points,after_cluster,mean,K,num_points);
		// printf("Calculando nuevos centroides!\n");
		//toc = clock();
		toc = omp_get_wtime();

		if(chkConvrg(after_cluster,before_clusters,num_points,tol)==0)
		{		
			
			//printf("Algoritmo K-means converge!\n");
			//printf("Iteraciones! : %d \n",cont);
			job_done=1;
		}
		else{
			//printf("No convergio!\n");
			for(i=0;i<num_points;i++)
				before_clusters[i]=after_cluster[i];
			
		}

		/*if(job_done==1)
			break;*/
		iter++;
		
	}
	//toc = omp_get_wtime();
	tspantemp = (double)(toc-tic);
	FILE* ofp=fopen(argv[2],"w");
	FILE* ofpc=fopen(fcentroid,"w");
	fprintf(ofp,"%d\n",K);
	fprintf(ofp,"%d\n",x);
	fprintf(ofp,"%d\n",y);
	for(i=0;i<K;i++)
		fprintf(ofpc,"%d,%d,%d,%d,%d\n",(int)mean[i]._r,(int)mean[i]._g,(int)mean[i]._b,(int)mean[i]._m,(int)mean[i]._n);
	for(i=0;i<num_points;i++)
		fprintf(ofp,"%d,%d,%d,%d,%d,%d\n",(int)points[i]._r,(int)points[i]._g,(int)points[i]._b,(int)points[i]._m,(int)points[i]._n,after_cluster[i]+1);
	fclose(ofp);
	fclose(ofpc);
	printf("Total Iteraciones = %d\n",iter);
	printf("OPENMP:Tiempo total transcurrido en la ejecucion. k= %d con tolerancia= %f en clusters : %f sec\n", K,tol,tspantemp);
	return 0;
}
