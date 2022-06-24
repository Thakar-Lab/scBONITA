#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#define STEP 100 //simsteps ROW
#define NODE 50000 //nodes COL
#define CELL 15000 

int updateBool(int currentNode, int *oldValue,int *nodeIndividual, int andNodes[7][3], int andNodeInvertList[7][3], int nodeStart, int nodeEnd)
{
// we update node by updating shadow and nodes then combining them to update or nodes.
    //update nodes with more than one input
    // first deal with case of simple logic without need of linear regression
    int counter =0;
    int orset[300];
    int andindex;
    int indindex;
    int addnode;
    int newval;

    // go through list of possible shadow and nodes to see which ones actually contribute
    for( indindex=nodeStart; indindex< nodeEnd; indindex++)
    {
        andindex=indindex-nodeStart;
        if ( nodeIndividual[indindex])
        {
            // if a shadow and contributes, compute its value using its upstream nodes
            // calculate value of first then use and to append rest in list of predecessors
            newval=(oldValue[andNodes[andindex][0]]!=andNodeInvertList[andindex][0]);

            for( addnode=1; addnode < 3; addnode++)
            {
                if(andNodes[andindex][addnode]>(-1)){newval=(newval && (oldValue[andNodes[andindex][addnode]]!=andNodeInvertList[andindex][addnode]));}
            }

            orset[counter]=newval;
            counter++;

        }
    }
    //combine the shadow and nodes with or operations
    newval=orset[0];
    int q;
    for(q = 1; q<counter; q++)
    {
        newval= (newval || orset[q]);
    }
    return newval;
}

void scSyncBool_old(int simData[NODE][STEP][CELL], int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[NODE][7][3], int andNodeInvertList[NODE][7][3], int simSteps, int *knockouts, int *knockins, int lenSamples, int binMat[NODE][CELL], int *nodePositions){

    // do simulation. individual specifies the particular logic rules on the model. params is a generic holder for simulation parameters.
    // set up data storage for simulation, add step 0
    // sim data is nodes * sim steps. Max sim steps hard coded to 200
    // iterate over number of steps necessary

    int step;
    int i;
    int nodeEnd;
    int temp;
    int nodeStart;
    int oldValue[NODE];
    int newValue[NODE];
    int cell;
    int nodePos;

    //for(i=0; i<nodeNum; i++){
    //  if(knockouts[i]){printf("Node %d: ", i);printf("Knocked out\n");}
    //  else if(knockins[i]){printf("Node %d: ", i);printf("Knocked in\n");}
    //}


    for(cell=0; cell<lenSamples; cell++){

        //printf("Cell: %d", cell);
        //printf("\nOriginal value:\n");

        for(i=0; i<nodeNum; i++){
            nodePos = nodePositions[i];
            //printf("%d ", nodePos);
            newValue[i]=binMat[nodePos][cell]; //[nodePos];
            //printf("%d ", newValue[i]);
            simData[i][0][cell]=temp;
        }
        //for(i=0; i<NODE; i++){
        //  newValue[i]=binMat[i][cell];
            //printf("%d ", binMat[i][cell]);
        //}
        //printf("\n");
        for(step=1; step < simSteps; ++step){
            for(i=0; i< nodeNum; i++){oldValue[i]=newValue[i];} //printf("%d", oldValue[i]);
            for(i=0; i< nodeNum; i++){
                //printf("\nNode: %d\n", i);
                //printf("\nNode: %d, andLenList: %d", i, andLenList[i]);
                //find start and finish for each node to update from the individualParse list
                if(knockouts[i]){temp=0;}
                else if(knockins[i]){temp=1;}
                else if (andLenList[i]==1){temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);} //;//printf("\n3catch ")
                else if (andLenList[i]==0){temp=oldValue[andNodes[i][0][0]];} //printf("3catch");
                else{
                    if (i==(nodeNum-1)){nodeEnd= indLen;} //if last node
                    else{nodeEnd=individualParse[i+1];}
                    nodeStart=individualParse[i];
                    //printf("\nNodestart: %d\n",nodeStart);
                    //printf("\nNodeend: %d\n",nodeEnd);
                    temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                }
                //newValue[i]=temp;
                simData[i][step][cell]=temp;
                //printf("%d", simData[i][step][cell]);
            }
            //printf("\n");
        }

        //printf("##########\n");
    }

}

void scSyncBool_python(int simData[NODE][STEP][CELL], int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[NODE][7][3], int andNodeInvertList[NODE][7][3], int simSteps, int *knockouts, int *knockins, int lenSamples, int binMat[NODE][CELL], int *nodePositions){

    int step;
    int i;
    int nodeEnd;
    int temp;
    int nodeStart;
    int oldValue[NODE];
    int newValue[NODE];
    int cell;
    int nodePos;

    //for(i=0; i<nodeNum; i++){
    //  printf("Node %d: ", i);
    //  printf("Knocked out: %d\n", knockouts[i]);
        //if(knockouts[i]==1){printf("Node %d: ", i);printf("Knocked out\n");}
        //else if(knockins[i]==1){printf("Node %d: ", i);printf("Knocked in\n");}
    //}

    for(cell=0; cell<lenSamples; cell++){
        for(i=0; i<nodeNum; i++){
            nodePos = nodePositions[i];
            newValue[i]=binMat[nodePos][cell];
            simData[i][0][cell]=newValue[i];
        }

        for(step=1; step < simSteps; step++){
            for(i=0; i< nodeNum; i++){oldValue[i]=newValue[i];} //simData[i][step-1][cell];}
            for(i=0; i< nodeNum; i++){
                if(knockouts[i]==1){temp=0;}
                else if(knockins[i]==1){temp=1;}
                if (andLenList[i]==1){temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);}
                else if (andLenList[i]==0){temp=oldValue[andNodes[i][0][0]];}

                else{
                    if (i==(nodeNum-1)){nodeEnd=indLen;}
                    else{nodeEnd=individualParse[i+1];}
                    nodeStart=individualParse[i];
                    temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                }
                newValue[i]=temp;
                simData[i][step][cell]=temp;
            }

        }
    }

}

//Given a binary matrix of M X N of integers, you need to return only unique rows of binary array

// A Trie node
typedef struct Node
{
    bool isEndOfCol;
    struct Node *child[2]; // Only two children needed for 0 and 1
} Node;


// A utility function to allocate memory for a new Trie node
Node* newNode()
{
    Node* temp = (Node *)malloc( sizeof( Node ) );
    temp->isEndOfCol = 0;
    temp->child[0] = temp->child[1] = NULL;
    return temp;
}

// Inserts a new matrix row to Trie.  If row is already
// present, then returns 0, otherwise inserts the row and
// returns 1

bool insert( Node** root, int (*M)[NODE], int row, int col )
{
    // base case
    if ( *root == NULL )
        *root = newNode();

    // Recur if there are more entries in this row
    if ( col < NODE )
        return insert ( &( (*root)->child[ M[row][col] ] ), M, row, col+1 );

    else // If all entries of this row are processed
    {
        // unique row found, return 1
        if ( !( (*root)->isEndOfCol ) )
            return (*root)->isEndOfCol = 1;

        // duplicate row found, return 0
        return 0;
    }
}


// A utility function to print a row
void printRow( int (*M)[NODE], int row )
{
    int i;
    for( i = 0; i < NODE; ++i )
        printf( "%d ", M[row][i] );
    printf("\n");
}

// The main function that prints all unique rows in a
// given matrix.
void findUniqueRows( int (*M)[NODE] )
{
    Node* root = NULL; // create an empty Trie
    int i;

    // Iterate through all rows
    for ( i = 0; i < STEP; ++i )
        // insert row to TRIE
        if ( insert(&root, M, i, 0) )
            // unique row found, print it
            printRow( M, i );
}

void findDuplicatedRows( int (*M)[NODE] )
{
    Node* root = NULL; // create an empty Trie
    int i;

    // Iterate through all rows
    for ( i = 0; i < NODE; ++i )
        // insert row to TRIE
        if (! insert(&root, M, i, 0))
            // duplicate row found, print it
            printRow( M, i );

}

struct Tuple {
    int first;
    int second;
};


struct Tuple getAttractCycle(int (*M)[NODE]){
    Node* root = NULL;
    int i;
    int k;
    int isCycle;
    int temp;
    for (i=0; i<STEP; ++i){
        if (! insert(&root, M, i, 0)){
            //find second occurrence of row
            //printf("Second occurrence: %d\n", i);
            //find first occurrence of row
            Node* root = NULL; // create an empty Trie
            for ( k = i; k > 0; --k ){
            //for ( k = 0; k < i; ++k ){
                if (! insert(&root, M, k, 0)){
                    //printf("First occurrence: %d\n", k);
                    //printRow( M, i );
                    //printRow( M, k );
                    //isCycle = i - k;

                    //for (temp = k; temp <=i; temp++){
                    //  printRow(M, temp);
                    //}
                    //return 0;
                    struct Tuple res = {k, i};
                    return res;
                }
            }
        }
    }
}


//void transpose(bool vals[][STEP], bool valsT[][NODE])
void transpose(int vals[NODE][STEP], int valsT[STEP][NODE])
{
    int i;
    int j;
    for (i = 0; i < STEP; i++) {
        for (j = 0; j < NODE; j++) {
            valsT[i][j] = vals[j][i];
        }
    }

}

struct Tuple getAttractCycle2(int valsT[STEP][NODE], int *nodePositions, int nodeNum){
    struct Tuple res = {0,0};
    int i;
    int j;
    int k;
    //Traverse through the matrix
    for(i=STEP-1; i >= 0; i--) {
        int flag=0;
        for(j=i-1; j>=0; j--) {
            flag=1;
            //for(k=0; k<=NODE;k++) {
            for(k=0; k < nodeNum; k++) {
                //if(valsT[i][k]!=valsT[j][k]){
                if(valsT[i][nodePositions[k]]!=valsT[j][nodePositions[k]]){
                    flag=0;
                }
            }
            if (flag == 1) {
                res.first = j;
                res.second = i;
                return res;
            }
        }
    }
    return res;
}

// Driver program to test above functions
void getAttract(int vals[STEP][NODE], int resSubmit[2], int *nodePositions, int nodeNum)
{
    struct Tuple res = getAttractCycle2(vals, nodePositions, nodeNum); //getAttractCycle(vals);
    resSubmit[0] = res.first;
    resSubmit[1] = res.second;
    int i;
    int j;
    //for(i=resSubmit[0]; i<=resSubmit[1]; i++){
    //  for (j = 0; j< nodeNum; j++){
    //      printf("%d\t", vals[i][j]);
    //  }
    //  printf("\n");
    //}
}


void scSyncBool(int simData[STEP][NODE], int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[NODE][7][3], int andNodeInvertList[NODE][7][3], int simSteps, int *knockouts, int *knockins, int lenSamples, int binMat[NODE][CELL], int *nodePositions, int errors[CELL], int localSearch, int importanceScores){

    int step;
    int i;
    int nodeEnd;
    int temp;
    int nodeStart;
    int oldValue[NODE];
    int newValue[NODE];
    int cell;
    int nodePos;
    int error[CELL];
    int totalError;
    int minError;

    //iteration over samples
    for(cell=0; cell<lenSamples; cell++){
        for(i=0; i<nodeNum; i++){
            nodePos = nodePositions[i];
            newValue[i]=binMat[nodePos][cell];
            simData[0][nodePos]=newValue[i];
        }

        for(step=1; step < simSteps; step++){

            for(i=0; i < nodeNum; i++){
                oldValue[i]=newValue[i];
            }
            for(i=0; i< nodeNum; i++){

                if(knockouts[i]==1){
                    temp=0;
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    continue;} else if(knockins[i]==1){
                    temp=1;
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    continue;
                }
                if (andLenList[i]==1){
                    temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    continue;}
                else if (andLenList[i]==0){
                    temp=oldValue[andNodes[i][0][0]];
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    continue;}
                else{
                    if (i==(nodeNum-1)){nodeEnd=indLen;}
                    else{nodeEnd=individualParse[i+1];}
                    nodeStart=individualParse[i];
                    temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    continue;
                }
                //if (simData[step][nodePos] > 1){
                //  printf("%d %d\n",simData[step][nodePos], binMat[nodePos][cell]);
                //}
            }
            importanceScores=importanceScores+temp;
        }

        //add error calculation function
        int resSubmit[2]; //positions of attractor in trajectory
        getAttract(simData, resSubmit, nodePositions, nodeNum); //get positions of attractor in trajectory

        minError=100000;
        //printf("%d %d \n",resSubmit[0], resSubmit[1]);
        //iterate over attractors
        for(temp=resSubmit[0]; temp<resSubmit[1]; temp++){
            //printf("Temp: %d\n", temp);
            //prepare attractor
            for(i=0; i<nodeNum; i++){
                nodePos = nodePositions[i];
                //if (simData[step][nodePos] > 0){
                //  printf("%d %d %d %d %d\n",temp, nodePos, nodeNum, simData[step][nodePos], binMat[nodePos][cell]);
                //}
                newValue[i] = simData[temp][nodePos]; //this is the attractor
                //printf("%d ", newValue[i]);
            }
            //printf("\n");

            //calculate error
            totalError=0;
            for(i=0; i<nodeNum; i++){
                nodePos = nodePositions[i]; //prepare cellInitValue from input data
                //error[i] = abs(newValue[i] - binMat[nodePos][cell]); //why is this always zero?
                //printf("Error: %f, Node: %d, newValue: %d, binMat: %d\n", error[i], i, newValue[i], binMat[nodePos][cell]);
                error[i] = abs(newValue[i] - binMat[nodePos][cell]);
                totalError=totalError+error[i];
            }
            //if (totalError > 0){
                //printf("Total error: %d\n", totalError);
            //}

            if (totalError < minError){
                minError = totalError;
                if (totalError <= 0.001*nodeNum && minError < 10000){ //this condition should change.
                    minError=0;
                    if (localSearch){
                        for(i=0; i<nodeNum; i++){
                            errors[i] = errors[i] + error[i];
                        }
                    } else {
                        errors[cell] = totalError;
                    }
                continue; //I want to break out of the loop that iterates over repeating states and skip to evaluating for the next sample
                }
            }  else {
                minError = minError;
            }

            }
        //finally
        if (localSearch){
            for(i=0; i<nodeNum; i++){
                errors[i] = errors[i] + error[i];
            }
        } else {
            errors[cell] = totalError;
        }
    }
}

void cluster(int simData[STEP][NODE], int resSubmit[2], int sampleIndex, int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[NODE][7][3], int andNodeInvertList[NODE][7][3], int simSteps, int *knockouts, int *knockins, int binMat[NODE][CELL], int *nodePositions){

    int step;
    int i;
    int nodeEnd;
    int temp;
    int nodeStart;
    int oldValue[NODE];
    int newValue[NODE];
    int nodePos;

    //simulate for only the specified sample
    for(i=0; i<nodeNum; i++){

        nodePos = nodePositions[i];
        newValue[i]=binMat[nodePos][sampleIndex];
        simData[0][nodePos]=newValue[i];
    }

    for(step=1; step < simSteps; step++){


        for(i=0; i < nodeNum; i++){
            oldValue[i]=newValue[i];
        }
        for(i=0; i< nodeNum; i++){

            if(knockouts[i]==1){
                temp=0;
                newValue[i]=temp;
                nodePos = nodePositions[i];
                simData[step][i]=temp;
                continue;} else if(knockins[i]==1){
                temp=1;
                newValue[i]=temp;
                nodePos = nodePositions[i];
                simData[step][nodePos]=temp;
                continue;
            }
            else if (andLenList[i]==1){
                temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);
                newValue[i]=temp;
                nodePos = nodePositions[i];
                simData[step][nodePos]=temp;
                continue;}
            else if (andLenList[i]==0){
                temp=oldValue[andNodes[i][0][0]];
                newValue[i]=temp;
                nodePos = nodePositions[i];
                simData[step][nodePos]=temp;
                continue;}
            else{
                if (i==(nodeNum-1)){nodeEnd=indLen;}
                else{nodeEnd=individualParse[i+1];}
                nodeStart=individualParse[i];
                temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                newValue[i]=temp;
                nodePos = nodePositions[i];
                simData[step][nodePos]=temp;
                continue;
            }

        }
    }
    getAttract(simData, resSubmit, nodePositions, nodeNum); //get positions of attractor in trajectory
}


void importanceScore(int simData[STEP][NODE], int *individual,int indLen, int nodeNum, int *andLenList, int *individualParse, int andNodes[NODE][7][3], int andNodeInvertList[NODE][7][3], int simSteps, int *knockouts, int *knockins, int lenSamples, int binMat[NODE][CELL], int *nodePositions, double *importanceScores){

    int step;
    int i;
    int nodeEnd;
    int temp;
    int nodeStart;
    int oldValue[NODE];
    int newValue[NODE];
    int cell;
    int nodePos;
    int error[CELL];
    int totalError;
    int minError;
    double attractorAverage_ko[nodeNum];
    double attractorAverage_ki[nodeNum];
    importanceScores[0] = 0.0;

    // print out ki and ko node
    //for(i=0; i< nodeNum; i++){
    //    attractorAverage_ko[i] = 0.0;
    //    attractorAverage_ki[i] = 0.0;
    //    if(knockouts[i]==1){
    //        printf("KO node: %d, nodePosition: %d\n", i, nodePositions[i]);
    //    }
    //    if(knockins[i]==1){
    //      printf("KI node: %i, nodePosition: %i\n", i, nodePositions[i]);
    //    }
    //}

    //iteration over samples
    for(cell=0; cell<lenSamples; cell++){
        //KNOCK-OUT
        //if(knockouts[0]){printf("Original data:\n");}
        for(i=0; i<nodeNum; i++){
            nodePos = nodePositions[i];
            newValue[i]=binMat[nodePos][cell];
            simData[0][nodePos]=newValue[i];
            //if(knockouts[0]){printf("%d ", binMat[nodePos][cell]);}
        }
        //if(knockouts[0]){printf("\n");}

        for(step=1; step < simSteps; step++){
            //if(knockouts[0]){printf("Step: %d\n", step);}
            for(i=0; i < nodeNum; i++){
                oldValue[i]=newValue[i];
            }
            for(i=0; i< nodeNum; i++){

                if(knockouts[i]){
                    temp=0;
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    }
                    //continue;}
                else if (andLenList[i]){
                    temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;}
                    //continue;}
                else if (andLenList[i]==0){
                    temp=oldValue[andNodes[i][0][0]];
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;}
                    //continue;}
                else{
                    if (i==(nodeNum-1)){nodeEnd=indLen;}
                    else{nodeEnd=individualParse[i+1];}
                    nodeStart=individualParse[i];
                    temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    //continue;
                }
                //if(knockouts[0]){printf("%d ", simData[step][nodePos]);}
            }
        }
        //if(knockouts[0]){printf("\n");}
        int ko_resSubmit[2]; //positions of attractor in trajectory
        //printf("KNOCKOUT\n");
        getAttract(simData, ko_resSubmit, nodePositions, nodeNum); //get positions of attractor in trajectory
        // if (knockouts[0]==1){printf("KO attractor:\n");printf("KO attractor loc: %i, %i\n", ko_resSubmit[0], ko_resSubmit[1]);}
        for(temp=ko_resSubmit[0]; temp<=ko_resSubmit[1]; temp++){
            //prepare attractor
            //if (knockouts[0]==1){printf("Temp: %d\n", temp);}
            for(i=0; i<nodeNum; i++){
                nodePos = nodePositions[i];
                //if (knockouts[0]==1){printf("%d ", simData[temp][nodePos]);}
                attractorAverage_ko[i] = attractorAverage_ko[i] + (double) (simData[temp][nodePos]/(ko_resSubmit[1] - ko_resSubmit[0])); //this is the attractor
            }
            //if (knockouts[0]==1){printf("\n");}
        }

    }
    //KNOCK - INS

    //iteration over samples
    for(cell=0; cell<lenSamples; cell++){
        //KNOCK-OUT
        //if(knockins[0]){printf("Original data:\n");}
        for(i=0; i<nodeNum; i++){
            nodePos = nodePositions[i];
            newValue[i]=binMat[nodePos][cell];
            simData[0][nodePos]=newValue[i];
            //if(knockins[0]){printf("%d ", binMat[nodePos][cell]);}
        }
        //if(knockins[0]){printf("\n");}

        for(step=1; step < simSteps; step++){
            //if(knockins[0]){printf("Step: %d\n", step);}
            for(i=0; i < nodeNum; i++){
                oldValue[i]=newValue[i];
            }
            for(i=0; i< nodeNum; i++){

                if(knockins[i]){
                    temp=0;
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    }
                    //continue;}
                else if (andLenList[i]){
                    temp=(oldValue[andNodes[i][0][0]]!=andNodeInvertList[i][0][0]);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;}
                    //continue;}
                else if (andLenList[i]==0){
                    temp=oldValue[andNodes[i][0][0]];
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;}
                    //continue;}
                else{
                    if (i==(nodeNum-1)){nodeEnd=indLen;}
                    else{nodeEnd=individualParse[i+1];}
                    nodeStart=individualParse[i];
                    temp=updateBool(i, oldValue, individual, andNodes[i], andNodeInvertList[i], nodeStart, nodeEnd);
                    newValue[i]=temp;
                    nodePos = nodePositions[i];
                    simData[step][nodePos]=temp;
                    //continue;
                }
                //if(knockins[0]){printf("%d ", simData[step][nodePos]);}
            }
        }
        //if(knockins[0]){printf("\n");}
        int ko_resSubmit[2]; //positions of attractor in trajectory
        //printf("KNOCKOUT\n");
        getAttract(simData, ko_resSubmit, nodePositions, nodeNum); //get positions of attractor in trajectory
        //if (knockins[0]==1){printf("\nKI attractor:\n");printf("KI attractor loc: %i, %i\n", ko_resSubmit[0], ko_resSubmit[1]);}
        for(temp=ko_resSubmit[0]; temp<=ko_resSubmit[1]; temp++){
            //prepare attractor
            //if (knockins[0]==1){printf("Temp: %d\n", temp);}
            for(i=0; i<nodeNum; i++){
                nodePos = nodePositions[i];
                //if (knockins[0]==1){printf("%d ", simData[temp][nodePos]);}
                attractorAverage_ko[i] = attractorAverage_ko[i] + (double) (simData[temp][nodePos]/(ko_resSubmit[1] - ko_resSubmit[0])); //this is the attractor
            }
            //if (knockins[0]==1){printf("\n");}
        }


        //average attractors - knock in
        //for(i=0; i <nodeNum; i++){
        //  temp = ki_resSubmit[1]- ki_resSubmit[0];
        //  //if (temp > 1){
        //  //  attractorAverage_ki[i] = attractorAverage_ki[i]/temp;
        //  //} else {attractorAverage_ki[i] = attractorAverage_ki[i];}
        //  attractorAverage_ki[i] = attractorAverage_ki[i]/temp;
        //}
    }

    //average attractors - knock in
    for(i=0; i <nodeNum; i++){
        attractorAverage_ki[i] = attractorAverage_ki[i]/lenSamples;
    }

    //average attractors - knock in
    for(i=0; i <nodeNum; i++){
        attractorAverage_ko[i] = attractorAverage_ko[i]/lenSamples;
    }

    //Calculate importance score - difference between attractorAverage_ko and attractorAverage_ki
    for (i=0; i <nodeNum; i++){
        //if (fabs(attractorAverage_ki[i] - attractorAverage_ko[i]) < lenSamples){
        //  importanceScores[0] = importanceScores[0] + fabs(attractorAverage_ki[i] - attractorAverage_ko[i]);
        //} else {
        //  importanceScores[0] = importanceScores[0] + 0.0;
        //}
        importanceScores[0] = importanceScores[0] + (fabs(attractorAverage_ki[i] - attractorAverage_ko[i]));
        //printf("KO avg: %f, KI avg: %f, Difference: %f\n", attractorAverage_ko[i], attractorAverage_ki[i], (fabs(attractorAverage_ki[i] - attractorAverage_ko[i])));
    }
    //average over number of cells
    //importanceScores[0] = importanceScores[0]/lenSamples;
    printf("IS: %f\n########\n\n", importanceScores[0]);
}