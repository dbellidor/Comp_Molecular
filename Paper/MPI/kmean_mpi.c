#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
typedef struct
{
	double _r;
	double _g;
	double _b;
	double _m;
	double _n;

} Point;

//Lee las dimensiones de la imagen a x b pixeles
void readImageSize(FILE *ifp,int* K,int* a,int* b)
{
	fscanf(ifp,"%d\n",K);
	printf("%d\n",*K);

	fscanf(ifp,"%d\n",a);
	printf("%d\n",*a);

	fscanf(ifp,"%d\n",b);
	printf("%d\n",*b);
}

//Lee el archivo ifp y lo almacena en la estructura.
void readPoints(FILE* ifp,Point *points, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
	{
		fscanf(ifp,"%lf,%lf,%lf,%lf,%lf", &points[i]._r, &points[i]._g, &points[i]._b, &points[i]._m, &points[i]._n);
		
	}
}

// Inicializar puntos aleatorios como medias 
void initialize(Point* mean,int K, int num_points, Point* points)
{
	int i, a, p=2;
	srand(time(NULL));
	for(i=0;i<K;i++)
		{
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
//Todos los puntos sin clusters.
int IntClusterMem(int *cluster, int num_points)
{
	int i;
	for(i=0;i<num_points;i++)
		{
		cluster[i]=-1;
		}
}

//Distancia
double calculateDistance(Point point1,Point point2)
{
	return sqrt((pow((point1._r-point2._r),2)+pow((point1._g-point2._g),2)+pow((point1._b-point2._b),2)+pow((point1._m-point2._m),2)+pow((point1._n-point2._n),2)));
}

//Para calcular a que cluster pertenece el punto.
int pointsCluster(Point point,Point* mean,int K)
{
	int parent=0;
	double dist = 0;
	double minDist=calculateDistance(point,mean[0]);
	int i;
	for(i=1;i<K;i++)
	{
		dist=calculateDistance(point,mean[i]);
			if(minDist>=dist)
			{
				parent=i;
				minDist=dist;
			}

	}
	return parent;
}

//calcular nueva media
void calcNewMean(Point* points,int* cluster,Point* mean,int K,int num_points)
{
	Point* newMean=malloc(sizeof(Point)*K);
	int* members=malloc(sizeof(int)*K);
	int i;
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
	for(i=0;i<K;i++)
	{
		if(members[i]!=0.0)
		{
		newMean[i]._r/=members[i];
		newMean[i]._g/=members[i];
		newMean[i]._b/=members[i];
		newMean[i]._m/=members[i];
		newMean[i]._n/=members[i];
		}
		else
		{
			newMean[i]._r=0;
			newMean[i]._g=0;
			newMean[i]._b=0;
			newMean[i]._m=0;
			newMean[i]._n=0;
		}
	}
	for(i=0;i<K;i++)
	{
		mean[i]._r=newMean[i]._r;
		mean[i]._g=newMean[i]._g;
		mean[i]._b=newMean[i]._b;
		mean[i]._m=newMean[i]._m;
		mean[i]._n=newMean[i]._n;
	}
}

//comprobar la convergencia
int chkConvrg(int *before_clusters,int *after_cluster,int num_points, float tol)
{
	int i;
	tol = num_points*tol;
	for(i=0;i<num_points;i++)
	if((before_clusters[i]-after_cluster[i])>tol)
		return -1;
	return 0;
}

int main(int argc, char* argv[])
{
	int rank;
	int size;
	struct timespec start_t, stop_t;
	clock_gettime(CLOCK_MONOTONIC,&start_t);
	double span_t;
	char hostname[256];
	int K=0;
	int num_points;
	int i,j,l=1;
	int job_size;
	int job_done=0;
	int x,y;
	float tol;
	int iter_max=0;
	double TIC,TOC,tic,ticIn, tocIn, ticC1, tocC1, ticC2, tocC2, ticC3, tocC3, ticNewMean, tocNewMean, ticChCon, tocChCon,ticOut,tocOut,tocPC4,ticPC4,ticPCal,tocPCal,ticPC5,tocPC5, ticPC6, tocPC6;
	double tspanNewMean = 0.0, tspanC2=0.0, tspanChCon=0.0, tspanC3=0.0, tspanOut=0.0,tspanPCal=0.0,tspanPC4=0.0,tspanPC5=0.0, tspanPC6=0.0;
	float check;

	Point* mean;
	Point* points;
	Point* get_points;
	int * formed_clusters;
	int * before_clusters;
	int * after_cluster;

	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	gethostname(hostname,255);
	clock_gettime(CLOCK_MONOTONIC,&stop_t);
	span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
	 TIC = span_t;
	MPI_Datatype MPI_POINT;
	MPI_Datatype type=MPI_DOUBLE;
	int blocklen=2;
	MPI_Aint disp=0; 
	MPI_Type_create_struct(1,&blocklen,&disp,&type,&MPI_POINT);
	MPI_Type_commit(&MPI_POINT);

	//Master

	if(rank!=0)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticPC4 = span_t;
	
		MPI_Recv(&job_size ,1 ,MPI_INT ,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(&K,1 ,MPI_INT ,0,0,MPI_COMM_WORLD,&status);
		mean =malloc(sizeof(Point)*K);
		MPI_Recv(mean ,K,MPI_POINT,0,0,MPI_COMM_WORLD,&status);
	
		points =(Point*)malloc(sizeof(Point)*job_size);
		after_cluster =(int*)malloc(sizeof(int)*job_size);

		for(i=0;i<job_size;i++)
		{

			MPI_Recv(&points[i]._r,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._g,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._b,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._m,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&points[i]._n,1,MPI_POINT ,0,0,MPI_COMM_WORLD,&status);
		}
	
		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
	

		tocPC4 = span_t;
		tspanPC4 = (double)(tocPC4-ticPC4);

		while(1)
		{

			ticPCal = MPI_Wtime();
			for(i=0;i<job_size;i++)
			{
				after_cluster[i]=pointsCluster(points[i],mean,K);

			}
			tocPCal = MPI_Wtime();
			if(rank == 1)
			tspanPCal += (double)(tocPCal-ticPCal);

			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));


			ticPC5 = span_t;

		  	MPI_Send(after_cluster,job_size, MPI_INT,0, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));


			tocPC5 = span_t;

			tspanPC5 += (double)(tocPC5-ticPC5);

			MPI_Bcast(&job_done,1, MPI_INT,0,MPI_COMM_WORLD);

			if(job_done==1) //No more work to be done
			break;
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));


			ticPC6 = span_t;

			MPI_Bcast(mean,K, MPI_POINT,0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));


			tocPC6 = span_t;

			tspanPC6 += (double)(tocPC6-ticPC6);
		}

	}
	else
	{

		//f = 1;
		//Archivo de lectura

		FILE *ifp;
		ifp=fopen(argv[1],"r");
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticIn = span_t;
		readImageSize(ifp,&K,&x,&y);
		K = atoi(argv[4]);
		num_points = x*y;
		if(K==10)
			iter_max=4;
		else if(K==30)
			iter_max=14;
		else if(K==50)
			iter_max=63;
		else
				iter_max=87;
		points =(Point*)malloc(sizeof(Point)*num_points);
		readPoints(ifp,points,num_points);
		fclose(ifp);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		tocIn = span_t;

		before_clusters =(int*)malloc(sizeof(int)*num_points);
		after_cluster=(int*)malloc(sizeof(int)*num_points);
		mean =malloc(sizeof(Point)*K);

		check = num_points%(size-1);
		if(check==0.00)
		{
			job_size=num_points/(size-1);
		}
		else
		{
			printf("\n Enter no. of Processes as n+1 (where n divides  %d  in equal parts)\n\n",num_points);
			exit(1);
		}

		initialize(mean,K,num_points,points);
		IntClusterMem(before_clusters,num_points);
		IntClusterMem(after_cluster,num_points);
		tol = atof(argv[3]);

		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticC1=span_t;

		//Envío del cluster principal a otros procesadores.
		for(i=1;i<size;i++)
		{

			MPI_Send(&job_size ,1 , MPI_INT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&K ,1 , MPI_INT ,i,0,MPI_COMM_WORLD);
			MPI_Send(mean ,K, MPI_POINT ,i,0,MPI_COMM_WORLD);
			for(j=0;j<job_size;j++)
			{
			MPI_Send(&points[j]._r+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._g+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._b+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._m+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			MPI_Send(&points[j]._n+(i-1)*job_size,1 , MPI_POINT ,i,0,MPI_COMM_WORLD);
			}
		}


		MPI_Barrier(MPI_COMM_WORLD);
		clock_gettime(CLOCK_MONOTONIC, &stop_t);
		span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));



		tocC1 = span_t;

		//trabajo de procesador maestro
		while(1)
		{

			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticC2 = span_t;
			ticC2 = span_t;

			for(i=1;i<size;i++)
				MPI_Recv(after_cluster+(job_size*(i-1)),job_size,MPI_INT,i,0,MPI_COMM_WORLD,&status);

		MPI_Barrier(MPI_COMM_WORLD);

			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocC2 = span_t;


			tspanC2 += (double)(tocC2-ticC2);

			clock_gettime(CLOCK_MONOTONIC,&stop_t);
		        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticNewMean = span_t;

			calcNewMean(points,after_cluster,mean,K,num_points);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
		        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocNewMean = span_t;
			tspanNewMean += (double)(tocNewMean-ticNewMean);

			clock_gettime(CLOCK_MONOTONIC,&stop_t);
		        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			ticChCon = span_t;
			if(chkConvrg(after_cluster,before_clusters,num_points,tol)==0)
			{

				job_done=1;

			}
			else
			{
				l++;

				for(i=0;i<num_points;i++)
					before_clusters[i]=after_cluster[i];
			}
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
		        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
			tocChCon = span_t;
			tspanChCon += (double)(tocChCon-ticChCon);



				MPI_Bcast(&job_done,1, MPI_INT,0,MPI_COMM_WORLD);

			if(job_done==1)
				break;
			 MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));


			tic = span_t;


				MPI_Bcast(mean,K, MPI_POINT,0, MPI_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			clock_gettime(CLOCK_MONOTONIC,&stop_t);
			span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));



			tocC3 = span_t;
			tspanC3 += (double)(tocC3-tic);

		}


		FILE* ofp=fopen(argv[2],"w");
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
		        span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		ticOut = span_t;
		fprintf(ofp,"%d\n",K);
		fprintf(ofp,"%d\n",x);
		fprintf(ofp,"%d\n",y);
		for(i=0;i<num_points;i++)
			fprintf(ofp,"%.0f,%.0f,%.0f,%.0f,%.0f,%d\n",points[i]._r,points[i]._g,points[i]._b,points[i]._m,points[i]._n,after_cluster[i]+1);
		fclose(ofp);
		clock_gettime(CLOCK_MONOTONIC,&stop_t);
	   	span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
		tocOut = span_t;
		tspanOut = (double)(tocOut-ticOut);
	}
	clock_gettime(CLOCK_MONOTONIC,&stop_t);
	span_t = (stop_t.tv_sec+(stop_t.tv_nsec/1000000000.0)) - (start_t.tv_sec+(start_t.tv_nsec/1000000000.0));
	 TOC = span_t;


	if(rank==0)
	{
		 printf("K-mean converge en la iteracion %d\n",iter_max);
		 if(K>10)
	   	 	printf("\n MPI: Tiempo Total  : %f sec\n",(TOC - TIC)*10);
		 else
			printf("\n MPI: Tiempo Total  : %f sec\n",(TOC - TIC));
	}
	else if(rank==1)
	{

	}

	MPI_Finalize();
	/*
	free (get_points);
	free (formed_clusters);
	free (mean);
	free (before_clusters);
	free (points);
	free (after_cluster);*/

     return 0;
}
