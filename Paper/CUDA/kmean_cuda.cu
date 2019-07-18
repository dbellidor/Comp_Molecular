#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cuda.h>

typedef struct {
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;
} Point;

#define CUDA_CALL(x) {if((x) != cudaSuccess){ \
	printf("CUDA error at %s:%d\n",__FILE__,__LINE__); \
	printf("  %s\n", cudaGetErrorString(cudaGetLastError())); \
	exit(EXIT_FAILURE);}}

//Leer las dimensiones de la imagen a x b pixeles  
void readImageSize(FILE *ifp,int* K,int* a,int* b) {
	fscanf(ifp,"%d\n",K);
	printf("%d\n",*K);

	fscanf(ifp,"%d\n",a);
	printf("%d\n",*a);

	fscanf(ifp,"%d\n",b);
	printf("%d\n",*b);
}

//Leer archivo ifp y lo almacena en el struct
void readPoints(FILE* ifp,Point *points, int num_points) {
	int i;
	for(i=0;i<num_points;i++) {
		fscanf(ifp,"%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
		
	}
}

//Inicializar puntos aleatorios como los medios (k numero de clusters) 
void initialize(Point* mean,int K, int num_points, Point* points) {
	int i, a, p=2;
	srand(time(NULL));
	for(i=0;i<K;i++) {
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
void IntClusterMem(int *cluster, int num_points) {
	int i;
	for(i=0;i < num_points; i ++) {
		cluster[i]=-1;
	}
}


//Para calcular a que cluster pertenece el punto. (k numero de clusters)
__global__ void pointsCluster(int* after_cluster_d, Point* point_d,Point* Dmean,int K, int x, int y) {
	

	int j, k, i;
	j = blockIdx.x*blockDim.x+threadIdx.x;
	k = blockIdx.y*blockDim.y+threadIdx.y;
	
	int parent=0;
	double dist=0;
	int t = (k*(x)+j);
	
	
	double minDist= sqrt((pow((point_d[t]._r-Dmean[0]._r),2)+pow((point_d[t]._g-Dmean[0]._g),2)+pow((point_d[t]._b-Dmean[0]._b),2)+pow((point_d[t]._m-Dmean[0]._m),2)+pow((point_d[t]._n-Dmean[0]._n),2)));
	
	for(i=1;i<K;i++) {
		dist = sqrt((pow((point_d[t]._r-Dmean[i]._r),2)+pow((point_d[t]._g-Dmean[i]._g),2)+pow((point_d[t]._b-Dmean[i]._b),2)+pow((point_d[t]._m-Dmean[i]._m),2)+pow((point_d[t]._n-Dmean[i]._n),2)));
			if(minDist>=dist) {
				parent=i;
				minDist=dist;
			}
	}
	after_cluster_d[t] = parent;
}


//Calcular nueva media
void calcNewMean(Point* points,int* cluster,Point* mean,int K,int num_points) {
	Point* newMean=(Point*)malloc(sizeof(Point)*K);
	int* members=(int*)malloc(sizeof(int)*(K));
	int i;
	for(i=0;i<K;i++) {
		members[i]=0;
		newMean[i]._r=0;
		newMean[i]._g=0;
		newMean[i]._b=0;
		newMean[i]._m=0;
		newMean[i]._n=0;
	}
	for(i=0;i<num_points;i++) {
		members[cluster[i]]++;
		newMean[cluster[i]]._r+=points[i]._r;
		newMean[cluster[i]]._g+=points[i]._g;
		newMean[cluster[i]]._b+=points[i]._b;
		newMean[cluster[i]]._m+=points[i]._m;
		newMean[cluster[i]]._n+=points[i]._n;
	}

	for(i=0;i<K;i++) {
		if(members[i]!=0.0) {
			newMean[i]._r/=members[i];
			newMean[i]._g/=members[i];
			newMean[i]._b/=members[i];
			newMean[i]._m/=members[i];
			newMean[i]._n/=members[i];
		}
		else {
			newMean[i]._r=0;
			newMean[i]._g=0;
			newMean[i]._b=0;
			newMean[i]._m=0;
			newMean[i]._n=0;
		}
	}

	for(i=0;i<K;i++) {
		mean[i]._r=newMean[i]._r;
		mean[i]._g=newMean[i]._g;
		mean[i]._b=newMean[i]._b;
		mean[i]._m=newMean[i]._m;
		mean[i]._n=newMean[i]._n;
	}
}

//Comprobamos la convergencia
//Comprueba que cada clúster de puntos permanece igual
int chkConvrg(int *before_clusters,int *after_cluster,int num_points, float tol) {
	int i;
	tol = num_points*tol;
	for(i=0;i<num_points;i++) {
	if(abs(before_clusters[i]-after_cluster[i])>tol) {
		//check = abs(before_clusters[i]-after_cluster[i]);
		//printf("check = %d, after_cluster[%d]=%d, before_clusters[%d]=%d\n",check,i,after_cluster[i],i,before_clusters[i]);
		return -1;
		}
	}

return 0;
}

int main(int argc, char* argv[]) {
	
	//Variables CPU
	int K;
	int num_points;
	int * before_clusters;
	int i;
	int job_done=0;
	int x,y,iter=0;

	Point* mean;
	Point* points;

	int * after_cluster;
	float tol;

	//Variables GPU
	Point* points_d;
	Point* mean_d;
	int * after_cluster_d;
	int * before_cluster_d;

	cudaEvent_t startinit, endinit, startmean, endmean, startcal, endcal, startindex, endindex;
	cudaEvent_t start1, end1;
	float timeinit, timemean, timecal, timeindex;
	float time1;

	//float totTime = 0;
	tol = atof(argv[3]);
	//iterations = atof(argv[3]);
	//printf("Ingrese tolerancia:  ");
	//scanf("%f",&tol);
	printf("Tolerancia = %f \n",tol);

	cudaEventCreate(&start1);
	cudaEventCreate(&end1);
	cudaEventRecord(start1, 0);

	//Leyendo archivo
	FILE *ifp;
	ifp=fopen(argv[1],"r");
	readImageSize(ifp,&K,&x,&y);
	K = atoi(argv[6]);
	num_points = x*y;
	int blockX=atoi(argv[4]);
	int blockY=atoi(argv[5]);

	//asignar memoria a la CPU
	points=(Point*)malloc(sizeof(Point)*num_points);
	readPoints(ifp,points,num_points);
	fclose(ifp);

	//printf("Entrada leida con exito \n");
	before_clusters=(int*)malloc(sizeof(int)*num_points);
	after_cluster=(int*)malloc(sizeof(int)*num_points);
	mean=(Point*)malloc(sizeof(Point)*K);

	//inicializando valores por defecto
	initialize(mean,K, num_points, points);
	IntClusterMem(before_clusters,num_points);
	IntClusterMem(after_cluster,num_points);


	CUDA_CALL(cudaMalloc((void**) &after_cluster_d, sizeof(int)*num_points));
	CUDA_CALL(cudaMalloc((void**) &before_cluster_d, sizeof(int)*num_points));
	CUDA_CALL(cudaMalloc((void**) &points_d, sizeof(Point)*num_points));
	CUDA_CALL(cudaMalloc((void**) &mean_d, sizeof(Point)*K));

	cudaEventCreate(&startinit);
	cudaEventCreate(&endinit);
	cudaEventRecord(startinit, 0);

	//copiar los puntos al device
	CUDA_CALL(cudaMemcpy(points_d, points, sizeof(Point)*num_points, cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(after_cluster_d, after_cluster, sizeof(int)*num_points, cudaMemcpyHostToDevice));

	cudaEventRecord(endinit, 0);
	cudaEventSynchronize(endinit);
	cudaEventElapsedTime(&timeinit, startinit, endinit);

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
	while(1) {

		
		cudaEventCreate(&startmean);
		cudaEventCreate(&endmean);
		cudaEventRecord(startmean, 0);
		
		//copiar los centroides iniciales al device
		CUDA_CALL(cudaMemcpy(mean_d, mean, sizeof(Point)*K, cudaMemcpyHostToDevice));
		cudaEventRecord(endmean, 0);
		cudaEventSynchronize(endmean);
		cudaEventElapsedTime(&timemean, startmean, endmean);
		
		//copia de memoria cuda
		//CUDA_CALL(cudaMemcpy(after_cluster_d, after_cluster, sizeof(int)*num_points, cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(before_cluster_d, before_clusters, sizeof(int)*num_points, cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(x_d, &x, sizeof(int), cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(y_d, &y, sizeof(int), cudaMemcpyHostToDevice));
		//CUDA_CALL(cudaMemcpy(K_d, &K, sizeof(int), cudaMemcpyHostToDevice));
		cudaEventCreate(&startcal);
		cudaEventCreate(&endcal);
		cudaEventRecord(startcal, 0);

		dim3 block(blockX, blockY);
		dim3 grid((x+blockX-1)/blockX, (y+blockY-1)/blockY);

		pointsCluster<<<grid,block>>>(after_cluster_d, points_d,mean_d,K,x,y);

		
		cudaDeviceSynchronize();
		cudaEventRecord(endcal, 0);
		cudaEventSynchronize(endcal);
		cudaEventElapsedTime(&timecal, startcal, endcal);

		cudaEventCreate(&startindex);
		cudaEventCreate(&endindex);
		cudaEventRecord(startindex, 0);
		CUDA_CALL(cudaMemcpy(after_cluster, after_cluster_d, sizeof(int)*num_points, cudaMemcpyDeviceToHost));
		cudaEventRecord(endindex, 0);
		cudaEventSynchronize(endindex);
		cudaEventElapsedTime(&timeindex, startindex, endindex);
		calcNewMean(points,after_cluster,mean,K,num_points);
		//printf("Nuevos centroides son calculados!\n");

		if(iter>iter_max) {
			// printf("El algoritmo kmeans converge con = %d!\n",iter);
			job_done=1;
		}

		else {
			//printf("No converge!\n");
			for(i=0;i<num_points;i++) {
				//printf("1 after_cluster[%d]=%d, before_clusters[%d]=%d\n",i,after_cluster[i],i,before_clusters[i]);

				before_clusters[i]=after_cluster[i];

				//printf("after_cluster[%d]=%d, before_clusters[%d]=%d\n",i,after_cluster[i],i,before_clusters[i]);
			}
		}

		if(job_done==1)
			break;
		++iter;
	}

	//Salida en archivos
	FILE* ofp=fopen(argv[2],"w");
	FILE* ofpc=fopen(fcentroid,"w");
	fprintf(ofp,"%d\n",K);
	fprintf(ofp,"%d\n",x);
	fprintf(ofp,"%d\n",y);
	for(i=0;i<K;i++)
		fprintf(ofpc,"%.0f,%.0f,%.0f,%.0f,%.0f\n",mean[i]._r,mean[i]._g,mean[i]._b,mean[i]._m,mean[i]._n);
	for(i=0;i<num_points;i++)
		fprintf(ofp,"%.0f,%.0f,%.0f,%.0f,%.0f,%d\n",points[i]._r,points[i]._g,points[i]._b,points[i]._m,points[i]._n,after_cluster[i]+1);

	fclose(ofp);
	fclose(ofpc);
	cudaEventRecord(end1, 0);
	cudaEventSynchronize(end1);
	cudaEventElapsedTime(&time1, start1, end1);

	printf("Total Iteraciones = %d\n",iter-1);
	printf("CUDA:Tiempo total transcurrido en la ejecucion. k= %d con tolerancia= %f en clusters : %f sec\n", K,tol,time1/1000);
	//printf("Tiempo Total : %f\t sec\n",time1/1000);

	CUDA_CALL(cudaFree(after_cluster_d));
	CUDA_CALL(cudaFree(mean_d));
	CUDA_CALL(cudaFree(points_d));

	free(before_clusters);
	free(mean);
	free(points);
	free(after_cluster);

	return 0;
}
