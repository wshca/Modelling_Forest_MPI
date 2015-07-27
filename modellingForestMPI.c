/*

Wolfgang Huettinger
StudentID 0685854
whuettin@uoguelph.ca

*/

#include <stdio.h>
#include <stdlib.h>
#include "pilot.h"

//new stuff by myself
#include <string.h> // testing of the strings


// no. of Workers
#define W 10

PI_PROCESS *Master;                    // master node
PI_PROCESS *Worker[W];		            // array of process pointers
PI_CHANNEL *MessageOut[W];             // send the needed data out to the neighbor
PI_CHANNEL *MessageIn[W];              // receive the needed data from the neighbor
PI_CHANNEL *stopSlaveLoop[W];
PI_CHANNEL *startCalculation[W];
PI_CHANNEL *firstRunSlave[W];
PI_CHANNEL *readDataIn[W];
PI_CHANNEL *readDataOnceSlave[W];
PI_CHANNEL *writeDataOut[W];
PI_CHANNEL *nextCycle[W];             // used to tell the slaves that the next cycle is to be calculated
PI_CHANNEL *dataExchangeOut[W];
PI_CHANNEL *dataExchangeIn[W];
PI_CHANNEL *finishCalc[W];
PI_CHANNEL *patchChange[W];

// number of interation the program should run
#define numberInteration 100

// grid size
#define gridX 100// if you change this value you have to change the %200d values, too
#define gridY 100


// life or death decision module
int modificationRule(int upperNbr, int lowerNbr, int leftNbr, int rightNbr)
{


  return 0;
}
// -----------------------------------------------------------------------------
// random number generator basic function
double randomgenerator(int index)
{
  double rnd;

  srand((unsigned)time(NULL)*(index+2));
  rnd = (double)rand() / (double)RAND_MAX;
  //printf("float %f \n",rnd);
  return rnd;
}

int randomInt(int index)
{
  int rnd;

 // srand((unsigned)time(NULL)*(index+2));
  rnd = random();//(int)((double)rand() / ((double)RAND_MAX+1 * index));
  return rnd;
}

// -----------------------------------------------------------------------------
//
double randomnormal(double variance,double mean, int index)
//double randomnormal(double mean, double variance)
{
  double x = randomgenerator(index);
  const double pi = 3.14159265358979323846264338327950288;
  double sigma = 0.0;
  double rnd = 0.0;
  double exponent = 0.0;
  double mu;

  mu = mean;
  sigma = sqrt(variance);

  // We write a function so we simulate the random-normal function in netlogo
  rnd = 1 / (sigma * sqrt(2 * pi));
  exponent =  -1 * pow( (x - mu),2) / (2 * pow(sigma,2));
  rnd = rnd * exp(exponent);

return rnd;
}

// -----------------------------------------------------------------------------
//
double computeDistance(double height, double sddP, double alpha, int index)
{
  int i = 0;
  double rnd = 0;
  double dispdist;
  // Original function out of nlogo
  // ifelse random-float 1 < sddP
  // [set dispdist ln( random-float 1) * (- alpha)]
  // [set dispdist ln( random-float 1) * (- [height] of myself * 2.5 ^ 1.5)]
  rnd = randomgenerator(index);
  //ifelese function
  if (rnd < sddP)
    {
    // set dispdist ln( random-float 1) * (- alpha)
    dispdist = log(rnd) * (-alpha);
    }
    else
    {
    // set dispdist ln( random-float 1) * (- [height] of myself * 2.5 ^ 1.5)
    // note that 'pow' is the C function for power to
    dispdist = log(rnd) * (- height * pow(2.5,1.5));
    };
  return dispdist;
}

int introSpecies(int species)
{


  return 0;
}

double treesDefineResource(int countTreesHere, double compRed, double T_Gd)
{
  int i = 0;
  int j = 0;
  double increment = 0.0;
  //double countTreesHere = 0.0;
  double exponent = 0.0;
  double euler = 2.71828182845904523536028747135266249775724709369995;
  double T_supComp;
  double temp;

  T_supComp = countTreesHere * compRed;
  exponent =  T_supComp * T_Gd;
  temp = pow(euler, -exponent);

  return temp;
}

double treesGrowheight(double T_G0, double T_height, double T_resource, double T_maxH)
{
  double increment;
  increment = T_G0 * log (T_maxH / T_height) * T_height * T_resource;
  return increment;
}
/*
// random number generator for the initialization of the field
int randomgenerator(int reseed)
{
  double rnd;
  double i;
  double j;
  int sendback;
  //Initialize the random number generator with the seed
  srand((unsigned)time(NULL)*reseed);

  // Call the C routine rand() for a random number
  // hand out the random number to the calling routine

  i = rand();
  j = RAND_MAX;
  rnd = i/j;

  if (rnd < 0.5)
  {
      sendback = 0;
  }
  else
  {
      sendback = 1;
  };

  //printf("rand fnct is %f \n",i);
  //printf("RAND_Max is %f \n",j);
  //printf("rnd is %d \n",sendback);
  return sendback;
}
*/

// this is the node the slaves build of - doing the test - distributed
int slaveNode(int index, void* arg2 )
{
    //local variables
    int gridSize = gridX*gridY/W+2;
    // working grid were the edge are predefined or filled by the neighboring slave
    int partedGrid[gridSize];

    //char string[sizeString];
    //int palinFound = 0;
    int i;
    int j;
    int k;
    int l;
    int stopSlave = 0;       //preinitialization is important and used
    //int stopCalculation = 0; //preinitialization is important and used
    int startCalc = 0;
    int readDataInOnce = 0;
    int nextCycleCalc = 0;
    int initGrid = 0;
    int dataExchangeOutInput = 0;
    int dataExchangeInInput = 0;
    int partedGridIn[2*gridX] = {0};
    int partedGridOut[2*gridX] = {0};


double totresources = 0; //initial number of resoucses is 0
double returninterval = 0;
int gaps = 0;
double deathRes = 0;
double deathInv = 0;
int ddInv = 0;
int ddRes = 0;
int lsRes = 0;
int lsInv = 0;
int resN = 500;
int invN = 500;

double resVar = 0.10; //was 0.10 so times 100 to fit integer
double resMd = 0.90; //was 0.90 so times 100 to fit integer
double resM0 = 0.90; //was 0.90 so times 100 to fit integer
double resGd = 1.30; //was 1.30 so times 100 to fit integer
double resG0 = 0.25; //was 0.25 so times 100 to fit integer
int resF = 10;
double resIntrinsicM = 0.010; //was 0.010 so times 1000 to fit integer
double resMaxH = 2500.0;
double resReproH = 0;

double invVar = 0.10; //was 0.10
double invMd = 0.50; //was 0.50
double invM0 = 0.41; //was 0.41
double invGd = 0.99; //was 0.99
double invG0 = 0.10; //was 0.10
int invF = 5;
double invIntrinsicM = 0.005; //was 0.0050
double invMaxH = 2300.0;
double invReproH = 300.0;
double resIntrinsic = 0;
double reproH = 0;

double compRed = 0.20; //was 0.20
double hcompRatio = 1.00; //was 1.0

//int tree = 0;
double firstRes = 0;
double firstInv = 0;

double introT = 20;
double introP = 10;

int ticks = 0; //0;
int ticksmax = 500; //500;

double pDisturb = 0.010; //was 0.010
double distRadius = 1.5; //was 1.5
double sddP = 0.99; //was 0.99
double alpha = 6; //was 6

double resNlocal;
double invNlocal;
int rndcounter = 1;

char storeChar[100];

long int treemax = 10000000;
int T_xcoord[treemax];
int T_ycoord[treemax];
int T_age[treemax];
double T_height[treemax];
int T_species[treemax];
double T_Md[treemax];
double T_M0[treemax];
double T_Gd[treemax];
double T_G0[treemax];
double T_increment[treemax];
double T_resource[treemax];
double T_maxH[treemax];
double T_reproH[treemax];
double T_supComp[treemax];
double dispDistance[treemax];
int T_F[treemax];

int patchX = gridX;//2000;
int patchY = gridY/W;//2000;
int localPatch[gridX*gridY];
//int i,j,k;
int tree;
int treeTemp;
int treeTemp1;
int seedIndex;
int tempInt;
int treeDeath = 0;
int countTreesHere[patchX][patchY];
for(j=0;j<patchX;j++){
  for(k=0;k<patchY;k++){
    countTreesHere[j][k] = 0;
  }
}

// preinitialize the trees with -1 for not existing
for(i=0;i<treemax;i++){
    T_xcoord[i]= -1;
    T_height[i]= 0;
}
  tree = 0;



    // start the forever loop
    while (stopSlave == 0)
    {
      if( PI_ChannelHasData(nextCycle[index]) == 1)
      {
        // get signal to do the next calculation step
        PI_Read(nextCycle[index],"%d",&nextCycleCalc);
      };

      // initialize the local grid with random data once
      if(initGrid == 0)
      {
        for(i=0;i<gridSize;i++)
        {
          partedGrid[i] = randomgenerator(i);
        };
        // reset the value so we only do it once
        initGrid = 1;
      };


      if(dataExchangeOutInput == 1)
      {
        // generate partedGridOut
        for (i=0;i<gridX;i++)
        {
          partedGridOut[i] = partedGrid[i]; //upper values
          partedGridOut[i+gridX] = partedGrid[i+gridSize-gridX]; //lower values
        };

        // write out the data to the neighboring slaves
        PI_Write(MessageOut[index], "%200d", partedGridOut);

        // reset the variable so we only do it once per cycle time
        dataExchangeOutInput = 0;
      };

      if(dataExchangeInInput == 1)
      {
      if( PI_ChannelHasData(MessageIn[index]) == 1){
        // read in the data from the neighboring slaves
        PI_Read(MessageIn[index], "%200d", partedGridIn);
        }

      if(PI_ChannelHasData(patchChange[index]) == 1){
        //PI_Read(patchChange[index], "%10000d", localPatch);
        };

        // reset the variable so we only do it once per cycle time
        dataExchangeInInput = 0;

       };

        // overwrite the grid with input data
        for (i=0;i<gridX;i++)
        {
          partedGrid[i] = partedGridIn[i];
          partedGrid[i+gridSize-gridX] = partedGridIn[i+gridX];
        };
      //};

      // do actuall calculations in this loop once for every time step
      if(nextCycleCalc == 1)
      {

      // so we calcultate the cells
      //for(i=gridX;i<gridSize-gridX;i++)
      //{
      //  partedGrid[i] = modificationRule(partedGrid[i-gridX],partedGrid[i+gridX],partedGrid[i-1],partedGrid[i+1]);
      //};

      //PI_Write(finishCalc[index], "%d", 1);


//for(ticks=0;ticks<ticksmax;ticks++)
//{

  // let it run
  // introduce species
  if (ticks > introT-1 && ticks < (introT+introP)) {
      if (T_xcoord[tree] != -1){  //a way to do recycling of the memory because dead trees are -1
          tree++;
      }
      else
      {

  //introSpecies(0);
    for (i=0;i<resN;i++) {

      //create trees resident
      T_xcoord[tree] = (randomInt(index) % patchX);
      T_ycoord[tree] = (randomInt(index) % patchY);
      T_species[tree] = 0;
      T_height[tree] = 1;
      T_age[tree] = 0;
      T_Gd[tree] = resGd;
      T_G0[tree] = resG0;
      T_maxH[tree] = resMaxH;
      T_reproH[tree] = resReproH;
      T_M0[tree] = resM0;
      T_Md[tree] = resMd;
      T_F[tree] = resF;
      tree++;

    }
    for (i=0;i<invN;i++) {
      //create trees invader
      T_xcoord[tree] = (randomInt(index) % patchX);
      T_ycoord[tree] = (randomInt(index) % patchY);
      T_species[tree] = 0;
      T_height[tree] = 1;
      T_age[tree] = 0;
      T_G0[tree] = invG0;
      T_Gd[tree] = invGd;
      T_M0[tree] = invM0;
      T_Md[tree] = invMd;
      T_maxH[tree] = invMaxH;
      T_reproH[tree] = invReproH;
      T_F[tree] = invF;
      tree++;
    }
      }
  };

  //now do the tree modifications of the trees
  treeTemp = tree;
  treeTemp1 = tree;
  for(i=0;i<treeTemp;i++)
  {

     if(T_xcoord[i] != -1){

      //get the number of trees in the grid for the next step (resources calculations)
      for(j=0;j<patchX;j++){
        for(k=0;k<patchY;k++){
          if (T_xcoord[i] == j && T_ycoord[i] == k){
            if((T_height[j*k] > hcompRatio*T_height[i]) && (j*k != i)) countTreesHere[j][k]++;
          }
        }
      }

      // call routine for the resources
      for(j=0;j<patchX;j++){
        for(k=0;k<patchY;k++){
          T_resource[i] = T_resource[i] + treesDefineResource(countTreesHere[j][k], compRed, T_Gd[i]);
        }
      }

      // call "routine" for the aging process
      if(T_height[i] > 1){ //height is bigger than 1cm
        T_age[i]++;
      }

      // do growing of the trees
      T_increment[i] = treesGrowheight(T_G0[i], T_height[i], T_resource[i], T_maxH[i]);
      T_height[i] = T_height[i] + T_increment[i];

      // call routine for the dispersal process
      if(T_height[i] > T_reproH[i]){
        // compute disposal radius
        dispDistance[i] = computeDistance(T_height[i], sddP, alpha, index);
//printf("dispDistance %d \n",((int)(dispDistance[i]*T_xcoord[i]) % patchX));
        //// seed the new trees
       // printf("%d %d \n",randomInt(index),T_F[i]);
        for(j=0;j<(randomInt(index) % T_F[i]);j++){
            tree++;
            //determine the next free storage index
              k=0;
              /*
              while(k<tree+1){
                if(T_xcoord[k] == -1){
                  seedIndex = k;//treeTemp1+j;
                  k = tree+1;
                }
                k++;
              }
*/
seedIndex = treeTemp1+j;
            T_xcoord[seedIndex] = ((int)(dispDistance[i]*T_xcoord[i]) % patchX);
            T_ycoord[seedIndex] = ((int)(dispDistance[i]*T_ycoord[i]) % patchY); //modify here for better program
            T_species[seedIndex] = T_species[i];
            T_height[seedIndex] = 1;
            T_age[seedIndex] = 0;
            T_G0[seedIndex] = T_G0[i];
            T_Gd[seedIndex] = T_Gd[i];
            T_M0[seedIndex] = T_M0[i];
            T_Md[seedIndex] = T_Md[i];
            T_maxH[seedIndex] = T_maxH[i];
            T_reproH[seedIndex] = T_reproH[i];
            T_F[seedIndex] = T_F[i];
        }//maybe some changes so we put the new tree in the slot of a dead one

      }

      // kill the trees
      if(randomgenerator(index) < T_M0[i] * exp ( T_increment[i] * -1 * T_Md[i])){
        if(T_species[i] == 0) ddRes++;
        if(T_species[i] == 1) ddInv++;
        countTreesHere[T_xcoord[i]][T_ycoord[i]]--; //reduce the count in the cell
        T_xcoord[i] = -1;
        T_ycoord[i] = -1; //just to be nice setting xccord to -1 would be sufficient
        treeDeath++;
      }

      // kill by disturbance
      if(randomgenerator(index)<pDisturb){
        for(k=0;k<tree;k++){
          if(T_species[k] == i){
            T_species[k] = -1;
          }
        }
      }


    };//end of the xcoord != -1 loop


/*


  disturbPatch();
  changePatch();
*/
//tree++;
  }//end of the for loop going over all trees
treeTemp = tree;
//}//end of the ticks loop - should be removed later
ticks++;

if (ticks == numberInteration-1){
    tempInt = 0;
    for(k=0;k<tree;k++){
      if(T_species[k] != -1 && T_height[k]>200) tempInt++;
    }
    printf("slave %d: trees %d %d\n",index,tempInt,tree);
}







    //  };
      //};



      nextCycleCalc = 2; // make sure we only do it once to avoid problems
      };

      // for the data exchange
      if( PI_ChannelHasData(dataExchangeOut[index]) == 1)
      {
        PI_Read(dataExchangeOut[index],"%d",&dataExchangeOutInput);
      };
      if( PI_ChannelHasData(dataExchangeIn[index]) == 1)
      {
        PI_Read(dataExchangeIn[index],"%d",&dataExchangeInInput);
      };

      // final end of the slave sub program is determined by this channel
      if( PI_ChannelHasData(stopSlaveLoop[index]) == 1)
      {
        PI_Read(stopSlaveLoop[index],"%d",&stopSlave);
      };

    };//end while exterior

   return 0;
}


// this is the master node which controls the process
// it reads the file and hands out the workload to the slaves
int masterNode(int index, void* arg2)
{
   //char string[sizeString] = {0};  /* initialized to zeroes */
   int i;
   int j;
   int k;
   int rc;
   int t;
   int iterationIndex;
   int lineNumber = 0;
   int lineNumberInput = 0; /* initialized to zeroes */
   int endOfFile = 0;
   int finishCalcIn[W] = {0};
   int finishWriteData[W] = {0};
   //FILE *fileInput = fopen("n1000000_p75.inp", "rb");       // open file to read the palindrome for the test
   FILE *fileOutput = fopen("output.txt", "w"); // open file to write out the results

   //int wholeGrid[gridX*gridY]; (unneeded can be avoided to save memory
   int partedGrid[gridX*gridY/W+2];
   int partedGridOut[W][2*gridX];

   int patches[gridX*gridY];
   double pDisturb = 0.010; //was 0.010
   double distRadius = 1.5; //was 1.5
   int disturb;

   // datastrcutre is a matrix which is modelled as a list of integers 0 and 1 (can be modified
   // to a bigger grid and more dimensions as more states
   // it works as follows:
   //
   // index:     0 ... 1*(Worker_i+1)*gridY/W*gridX-1
   //   for the indeces for the Worker 0
   // index:     (Worker_i+1)*gridX ... 2*(Worker_i+1)*gridY/W*gridX-1
   //   for the indeces for the Worker 1
   // and so on up to
   // index:     (Worker_i+1)*gridX ... NumberOfWorker*(Worker_i+1)*gridY/W*gridX-1
   //
   // to my knowledge this is the most efficient way to store random data - there can be modification done
   // by change int to bool
   iterationIndex=0;
   i = 0;
   while (iterationIndex<numberInteration)
   {
     // make the gaps in the patches every iteration
     for(j=0;j<randomInt(i) % (gridX*gridY);j++){
       if(pDisturb>randomgenerator(i)){
         patches[j-1]=0;
         patches[j]=0;
         if(j>gridX && j<(gridX*gridY)-gridX){
           patches[j-gridX];
           patches[j+gridX];
         }
       }
     }


     for (k=0;k<W;k++)
     {
       //PI_Write( patchChange[k], "%10000d", patches); //note here is the patch size in use!
       PI_Write( dataExchangeOut[k], "%d", 1);
       if( PI_ChannelHasData(MessageOut[k]) == 1)
       {
         PI_Read(MessageOut[k], "%200d", &partedGridOut[k]); //%200d -> the flexible size seems not to be working in pilot
         finishWriteData[k] = 1;
         i = i++;

         if (i == W-1)
         {
           iterationIndex++;
           i = 0;
         };

       };
     };
     t = 0;
     k = 0;

     while (k<W)
     {

       if (finishWriteData[k] == 1)
       {
           t = t++;
       };


       PI_Write( MessageIn[k], "%200d", partedGridOut[W-k]);
       PI_Write( dataExchangeIn[k], "%d", 1);
       PI_Write( nextCycle[k], "%d", 1);

       if( PI_ChannelHasData(finishCalc[k]) == 1)
       {
         PI_Read(finishCalc[k],"%d",&finishCalcIn);

         if (finishCalcIn[k] == 1)
         {
           i = i + 1;
           //finishCalcIn = 0;
         };

         if (i == W-1)
         {

             i = 0;
         };

//printf("i is %d \n", i);
       };
       k++;
     };




     //printf("now we are in year %d \n",iterationIndex);





   };




   // stop the slaves
   for (k=0;k<W;k++)
   {
     // send signal to slave to finish calculations so the program can finish
     PI_Write( stopSlaveLoop[k], "%d", 1);
   };

   return 0;
}


int main( int argc, char *argv[] )
{
	int i;
    long int timestart;
    long int timeend;
    timestart = time(NULL);
	// Pilot configuration phase; return no. of processes available
	PI_Configure( &argc, &argv );

    // create Master process to be more clean
    Master = PI_CreateProcess( masterNode, 0, NULL);

	// create Worker processes (with arg i) and to/from channels
    for ( i=0; i<W; i++ )
    {
	  Worker[i] = PI_CreateProcess( slaveNode, i, NULL );
      stopSlaveLoop[i] = PI_CreateChannel( Master, Worker[i]);
      nextCycle[i] = PI_CreateChannel( Master, Worker[i]);
      startCalculation[i] = PI_CreateChannel( Master, Worker[i]);
      firstRunSlave[i] = PI_CreateChannel( Master, Worker[i]);
      readDataIn[i] = PI_CreateChannel( Master, Worker[i]);
      readDataOnceSlave[i] = PI_CreateChannel( Master, Worker[i]);
      writeDataOut[i] = PI_CreateChannel( Worker[i], Master);
      finishCalc[i] = PI_CreateChannel( Worker[i], Master);
      dataExchangeOut[i] = PI_CreateChannel( Master, Worker[i]);
      dataExchangeIn[i] = PI_CreateChannel( Master, Worker[i]);
      patchChange[i] = PI_CreateChannel( Master, Worker[i]);
    };

    // create the connection channels - note the number is W-1 instead of W!
    for ( i=0; i<W; i++ )
    {
      MessageOut[i] = PI_CreateChannel( Worker[i], Master );
      MessageIn[i] = PI_CreateChannel( Master, Worker[i] );
    };

    // start execution (slaveNode gets control on its processor)
    PI_StartAll();

    // stop the execution
    PI_StopMain(0);

    // calculate the time used in this program run
    timeend = time(NULL);
    timestart = timeend - timestart;

    // send a small message that the program is finished to be nice
    printf("Elapsed time is %ld .\n", timestart);
    printf("Program executed and finished without errors.");

    return 0;
}
