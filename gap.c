/******************************************************************************
  A template program for developing a GAP solver. Subroutines to read instance
  data and compute the cost of a given solution are included.

  This program can also be used to compute the cost and check the feasibility
  of a solution given from a file. The format of a file is:
  for each job j from 1 to n in this order, the index of the agent (the value
  should be given as values from [1, m]) to which j is assigned. For example,
  if n=4 and m=3, and jobs 1, 2, 3 and 4 are assigned to agents 2, 1, 3 and 1,
  respectively, then the data in the file should be as follows:  2 1 3 1.

  NOTE: Index i of agents ranges from 0 to m-1, and
        index j of jobs   ranges from 0 to n-1 in the program,
	while in the solution file,
	index i of agents ranges from 1 to m, and
        index j of jobs   ranges from 1 to n in the program.
	Sorry for the confusion.

  If you would like to use various parameters, it might be useful to modify
  the definition of struct "Param" and mimic the way the default value of
  "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "cpu_time.c"

/***** default values of parameters ******************************************/
#define	TIMELIM	300	/* the time limit for the algorithm in seconds */
#define	GIVESOL	0	/* 1: input a solution; 0: do not give a solution */

typedef struct {
  int		timelim;	/* the time limit for the algorithm in secs. */
  int		givesol;	/* give a solution (1) or not (0) */
  /* Never modify the above two lines.  */
  /* You can add more components below. */
} Param;			/* parameters */

typedef struct {
  int	n;	/* number of jobs */
  int	m;	/* number of agents */
  int	**c;	/* cost matrix c_{ij} */
  int	**a;	/* resource requirement matrix a_{ij} */
  int	*b;	/* available amount b_i of resource for each agent i */
} GAPdata;	/* data of the generalized assignment problem */

typedef struct {
  double	timebrid;	/* the time before reading the instance data */
  double	starttime;	/* the time the search started */
  double	endtime;	/* the time the search ended */
  int		*bestsol;	/* the best solution found so far */
  /* Never modify the above four lines. */
  /* You can add more components below. */
} Vdata;		/* various data often necessary during the search */

typedef struct {
  int b;
  float bres; /* rest amount of resouce */
  int *aagt; /* array of resource requirement a_{ij} for each agent i */
  int *aind; /* array of index sorted for resource requirement a_{ij} for each agent i */
  int current; /* current position */
} Agent;

/*************************** functions ***************************************/
void copy_parameters(int argc, char *arcv[], Param *param);
void read_instance(GAPdata *gapdata);
void prepare_memory(Vdata *vdata, GAPdata *gapdata);
void copy_instance(GAPdata *gapdata_o, GAPdata *gapdata_c);
void free_memory(Vdata *vdata, GAPdata *gapdata_o, GAPdata *gapdata_c);
void read_sol(Vdata *vdata, GAPdata *gapdata);
void recompute_cost(Vdata *vdata, GAPdata *gapdata);
void *malloc_e(size_t size);
void my_algorithm(Vdata *vdata, GAPdata *gapdata, Param *param);
int minIndex(float nums[], int n);

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_cost(Vdata *vdata, GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */

  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  printf("recomputed cost = %d\n", cost);
  if(penal>0){
    printf("INFEASIBLE!!\n");
    printf(" resource left:");
    for(i=0; i<gapdata->m; i++){printf(" %3d", rest_b[i]);}
    printf("\n");
  }
  printf("time for the search:       %7.2f seconds\n",
	 vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
	 vdata->starttime - vdata->timebrid);

  free((void *) rest_b);
}

/***** read a solution from STDIN ********************************************/
void read_sol(Vdata *vdata, GAPdata *gapdata)
{
  int	j;		/* index of jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */

  for(j=0; j<gapdata->n; j++){
    fscanf(fp, "%d", &value_read);
    /* change the range of agents from [1, m] to [0, m-1] */
    vdata->bestsol[j] = value_read - 1;
  }
}

/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory(Vdata *vdata, GAPdata *gapdata)
{
  int j;

  vdata->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  /* the next line is just to avoid confusion */
  for(j=0; j<gapdata->n; j++){vdata->bestsol[j] = 0;}
}

/***** free memory space *****************************************************/
/***** Feel free to modify this subroutine. **********************************/
void free_memory(Vdata *vdata, GAPdata *gapdata_o, GAPdata *gapdata_c)
{
  free((void *) vdata->bestsol);
  free((void *) gapdata_o->c[0]);
  free((void *) gapdata_o->c);
  free((void *) gapdata_o->a[0]);
  free((void *) gapdata_o->a);
  free((void *) gapdata_o->b);
  free((void *) gapdata_c->c[0]);
  free((void *) gapdata_c->c);
  free((void *) gapdata_c->a[0]);
  free((void *) gapdata_c->a);
  free((void *) gapdata_c->b);
}


/***** copy the instance data ************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void copy_instance(GAPdata *gapdata_o, GAPdata *gapdata_c) {
  int	i, j;		/* indices of agents and jobs */

  /* copy the the numbers of agents and jobs */
  gapdata_c->m = gapdata_o->m;
  gapdata_c->n = gapdata_o->n;
  
  /* initialize memory */
  gapdata_c->c    = (int **) malloc_e(gapdata_c->m * sizeof(int *));
  gapdata_c->c[0] = (int *)  malloc_e(gapdata_c->m * gapdata_c->n * sizeof(int));
  for(i=1; i<gapdata_c->m; i++){gapdata_c->c[i] = gapdata_c->c[i-1] + gapdata_c->n;}
  gapdata_c->a    = (int **) malloc_e(gapdata_c->m * sizeof(int *));
  gapdata_c->a[0] = (int *)  malloc_e(gapdata_c->m * gapdata_c->n * sizeof(int));
  for(i=1; i<gapdata_c->m; i++){gapdata_c->a[i] = gapdata_c->a[i-1] + gapdata_c->n;}
  gapdata_c->b    = (int *)  malloc_e(gapdata_c->m * sizeof(int));

  /* copy the cost coefficients */
  for(i=0; i<gapdata_o->m; i++){    
    for(j=0; j<gapdata_o->n; j++){
      gapdata_c->c[i][j] = gapdata_o->c[i][j];
    }
  }

  /* copy the resource consumption */
  for(i=0; i<gapdata_o->m; i++){
    for(j=0; j<gapdata_o->n; j++){
      gapdata_c->a[i][j] = gapdata_o->a[i][j];
    }
  }

  /* copy the resource capacity */
  for(i=0; i<gapdata_o->m; i++){    
    gapdata_c->b[i] = gapdata_o->b[i];
  }

}

/***** read the instance data ************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_instance(GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */

  /* read the number of agents and jobs */
  fscanf(fp, "%d", &value_read);	/* number of agents */
  gapdata->m = value_read;
  fscanf(fp,"%d",&value_read);		/* number of jobs */
  gapdata->n = value_read;

  /* initialize memory */
  gapdata->c    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->c[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->c[i] = gapdata->c[i-1] + gapdata->n;}
  gapdata->a    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->a[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->a[i] = gapdata->a[i-1] + gapdata->n;}
  gapdata->b    = (int *)  malloc_e(gapdata->m * sizeof(int));

  /* read the cost coefficients */   
  for(i=0; i<gapdata->m; i++){    
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->c[i][j] = value_read;
    }
  }

  /* read the resource consumption */
  for(i=0; i<gapdata->m; i++){
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->a[i][j] = value_read;
    }
  }

  /* read the resource capacity */
  for(i=0; i<gapdata->m; i++){    
    fscanf(fp,"%d", &value_read);
    gapdata->b[i] = value_read;
  }
}

/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters(int argc, char *argv[], Param *param)
{
  int i;

  /**** copy the parameters ****/
  param->timelim = TIMELIM;
  param->givesol = GIVESOL;
  /**** read the parameters ****/
  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./gap [param_name, param_value] [name, value]...\n");
    exit(EXIT_FAILURE);}
  else{
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0) param->timelim = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0) param->givesol = atoi(argv[i+1]);
    }
  }
}

/***** malloc with error check ***********************************************/
void *malloc_e( size_t size ) {
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : Not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}

/***** detect index with minimize value **************************************/
int minIndex(float nums[], int n) {
    float min_value; 
    int min_index; 
    int i;

    min_value = nums[0];
    min_index = 0;

    for (i = 0; i < n; i++) {
        if (nums[i] < min_value) {
            /* 最小値よりもnums[i]の方が小さければ最小値を更新 */
            min_value = nums[i];
            min_index = i;
        }
    }
    return min_index;
}

/***** qsort用比較関数 *********************************************************/
int compare(const void *n1, const void *n2) {
    if (*(int *)n1 > *(int *)n2)
	{
		return 1;
	}
	else if (*(int *)n1 < *(int *)n2)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}

/***** write your algorithm here *********************************************/
void my_algorithm(Vdata *vdata, GAPdata *gapdata, Param *param) {

  /* Here is an example of repeating a loop until the timelimit is reached
     or the maximum iteration number is reached. Uncomment below and try it.*/

  // int itr, itr_max, j;
  // itr_max = 500000000; 
  // for(itr = 0; itr < itr_max && cpu_time() - vdata->starttime <= param->timelim; itr++) {
  //   for(j = 0; j < gapdata->n; j++) {
  //     vdata->bestsol[j] = j % gapdata->m;
  //   }
  // }
  // printf("the iteration count after the loop ends: %d\n", itr);


  /*
    Write your program here. Of course you can add your subroutines
    outside of my_algorithm(). The instance data is stored in "gapdata".
	gapdata->n		number of jobs n
	gapdata->m		number of agents m
	gapdata->c[i][j]	cost c_{ij} 
	gapdata->a[i][j]	resource requirement a_{ij} 
	gapdata->b[i]		available amount b_i of resource at agent i
    Note that i ranges from 0 to m-1, and j ranges from 0 to n-1.
    Store your best solution in vdata->bestsol, then "recompute_cost" will
    compute its cost and its feasibility. The format of vdata->bestsol is:
    For each job j from 0 to n-1 in this order, the index of the agent 
    (the value should be given as values from [0, m-1]) to which j is
    assigned. For example, if n=4 and m=3, and jobs 0, 1, 2 and 3 are
    assigned to agents 1, 0, 2 and 0, respectively, then vdata->bestsol
    should be as follows:  
	vdata->bestsol[0] = 1
	vdata->bestsol[1] = 0
	vdata->bestsol[2] = 2
	vdata->bestsol[3] = 0.
  */

  //Flagの作成
  int* isAssigned = (int*) malloc_e(gapdata->n * sizeof(int));
  for (int i = 0; i < gapdata->n; ++i){
    isAssigned[i] = 0;
  }

  //構造体agentの作成
  Agent *agent = (Agent*) malloc_e(gapdata->m * sizeof(Agent));
  //agentの初期化
  for (int i = 0; i < gapdata->m; ++i){
    agent[i].aagt = (int*) malloc_e(gapdata->n * sizeof(int));
    agent[i].aind = (int*) malloc_e(gapdata->n * sizeof(int));
    agent[i].b = gapdata->b[i];
    agent[i].bres = gapdata->b[i];
    agent[i].current = 0;
    for (int j = 0; j < gapdata->n; ++j){
      agent[i].aagt[j] = gapdata->a[i][j];
      agent[i].aind[j] = j;
    }
  }
  //agentごとの要求資源量をソート
  int tmp;
  for (int i = 0; i < gapdata->m; ++i){
    for (int j = 0; j < gapdata->n; j++){
      for (int k = j + 1; k < gapdata->n; k++){
        if (agent[i].aagt[agent[i].aind[j]] > agent[i].aagt[agent[i].aind[k]])
        {
          tmp = agent[i].aind[j];
          agent[i].aind[j] = agent[i].aind[k];
          agent[i].aind[k] = tmp;
        }
      }
	  }
  }

  for(int x = 0; x < gapdata->n; ++x){
    //各エージェントに対して最小になりうるa/bresの配列を計算
    float *possiblemin = (float*)malloc_e(gapdata->m * sizeof(float));
    for (int i = 0; i < gapdata->m; ++i){
      possiblemin[i] = (float)agent[i].aagt[agent[i].aind[agent[i].current]]/agent[i].bres;
    }
    // for(int i = 0; i < gapdata->m; ++i){printf("possibleMin[] =  %f ",possiblemin[i]);}

    //a/bresが最小になるエージェントを特定
    int minAgent = minIndex( possiblemin, gapdata->m );

    //currentを進める
    agent[minAgent].current += 1;
    while(isAssigned[agent[minAgent].aind[agent[minAgent].current]] == 1){
      agent[minAgent].current += 1;
    }

    //bestsolに割当
    int minJob = agent[minAgent].aind[agent[minAgent].current];
    // printf("\nminAgent = %d minJob = %d\n", minAgent, minJob);
    vdata->bestsol[minJob] = minAgent;

    //bresの再計算
    agent[minAgent].bres = agent[minAgent].bres - agent[minAgent].aagt[minJob];
    if (agent[minAgent].bres <= 0){
      agent[minAgent].bres = 0.000000001;
      printf("over!!!");
    }
    // printf("bres = %f\n", agent[minAgent].bres);

    //Flagを立てる 
    isAssigned[minJob] = 1;
  }

  printf("Flag");
  for (int i = 0; i < gapdata->n; ++i ){ printf(" %d", isAssigned[i]);}
  printf("\n");
  // printf("bestsol\n");
  // for (int i = 0; i < gapdata->n; ++i){
  //   printf("job[%d] %d\n", i,vdata->bestsol[i]);
  // }

  free(isAssigned);


}

/***** main ******************************************************************/
/***** NEVER MODIFY THIS MAIN ROUTINE! ***************************************/
int main(int argc, char *argv[])
{
  Param		param;		/* parameters */
  GAPdata	gapdata_org;	/* the original GAP instance data */
  GAPdata	gapdata_cpd;	/* copied GAP instance data */
  Vdata		vdata;		/* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_instance(&gapdata_org);
  copy_instance(&gapdata_org, &gapdata_cpd);
  prepare_memory(&vdata, &gapdata_cpd);
  if(param.givesol==1){read_sol(&vdata, &gapdata_cpd);}
  vdata.starttime = cpu_time();
  my_algorithm(&vdata, &gapdata_cpd, &param);
  vdata.endtime = cpu_time();
  recompute_cost(&vdata, &gapdata_org);
  free_memory(&vdata, &gapdata_org, &gapdata_cpd);

  return EXIT_SUCCESS;
}
