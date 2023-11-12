/*
1)Consider how social dynamics change over time. For example, as individuals age, the influence of "parent" cells might decrease.

!Findings: NSE - There is a critical point where neighbour influenced transitions swings from almost all up to almost all down

2) Birth and Death: Implement birth and death processes, where new individuals
 are born with the average status of their neighbors and individuals die based on age or other factors, leaving vacancies in the grid.
*/
double SA_injury = 0.1;
double GER_injury = 0.005;
double SA_schol =0.01;
double GER_schol = 0.1;
#define NSE_glo  1.5 
#define INJURY SA_injury
#define SCHOLARSHIP SA_schol
int countInjury =0;
int countScholar = 0;

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <dirent.h>
#include <string.h>

// Define constants for the grid size, number of states, etc.
#define GRID_WIDTH 50
#define GRID_HEIGHT 50
#define NUM_EDUCATION_STATES 5
#define NUM_SOCIOECONOMIC_STATES 5
#define NUM_EVENTS 3

#define TIME_STEPS 5000

#define INTERVAL_OF_GRID_PLOTS 100

//#define POVERTY_THRESHOLD 1058
//#define MINIMUM_WAGE 4067


#define NEIGHBORHOOD_SIZE 8

//MC constants
// Define probabilities
//#define PROB_SCHOLARSHIP 0.01 // 1% chance
//#define PROB_INJURY 0.005     // 0.5% chance

// Define the impact of events
#define SCHOLARSHIP_IMPACT 1  // E.g., increases education level by 1
#define INJURY_IMPACT -1      // E.g., decreases socioeconomic status by 1




// Define a structure for an individual in the grid
typedef struct {
    int education_level;
    int socioeconomic_status;
    int age;
    // Add more attributes as needed
} Individual;


// Define a struct to hold all metrics for a single timestep
typedef struct {
    double giniCoefficient;
    double averageEducationLevel;
    double averageSocioeconomicStatus;
    double povertyRate;
    // ... other metrics ...
} Metrics;


typedef struct {
    float avg_education;  // Average education level of the neighbors
    float avg_economic;   // Average economic status of the neighbors
} NeighborAverages;



// Define a grid of individuals.
typedef Individual Grid[GRID_HEIGHT][GRID_WIDTH];

void equilibrateSystem(Grid grid, double nse_value, int equilibration_steps);

//
void simulateLifeEvents(Individual *individual,double prob_scholarship,double prob_injury);

// Function prototypes
void initializeGrid(Grid grid);
void evolveCellularAutomata(Grid grid,double NSE);
void performMarkovChain(Grid grid);
void performMonteCarloSampling(Grid grid,double prob_scholarship,double prob_injury);
void recordData(Grid grid);
void analyzeData();

// Function prototypes for metric calculations
Metrics calculateMetrics(Grid grid);
void writeMetricsToFile(Metrics* metricsArray, int numTimesteps, FILE* outputFile);


void updateState(Individual *individual);
double calculateGiniCoefficient(Grid grid);
double calculateAverageEducationLevel(Grid grid);
double calculateAverageSocioeconomicStatus(Grid grid);
double calculatePovertyRate(Grid grid);
// ... any other function prototypes ...

//migrating functions
NeighborAverages calculate_neighbor_averages(Grid grid, int x, int y);
// Swaps 2 individuals in the grid
void swap_individuals(Grid grid, int x1, int y1, int x2, int y2);
// Method which orchestrates the migration of an individual
int migrate(Grid grid, Individual *individual, int x, int y);


void education_makes_wealth(Grid grid,int x,int y);



//just for the counts
typedef struct {
    int countEduUp;
    int countEduDown;
    int countEcoUp;
    int countEcoDown;
} UpdateCounts;

// Function prototypes
void updateCellBasedOnNeighbors(Grid grid, int x, int y,double NSE);
void new_blood(Grid grid);


UpdateCounts* totalGlobal;

//double education_effct = 0.0;

// Transition matrix where transitionMatrix[i][j] is the probability of moving from state i to state j
double education_transitionMatrix[NUM_EDUCATION_STATES][NUM_EDUCATION_STATES] = {
    // From primary school
    {0.7, 0.22, 0.07, 0.01, 0.0},
    // From high school
    {0.05, 0.74, 0.15, 0.05, 0.01},
    // From diploma
    {0.01, 0.05, 0.73, 0.15, 0.01},
    // From undergraduate
    {0.02, 0.03, 0.05, 0.80, 0.1},
    // From postgraduate
    {0.0, 0.00, 0.02, 0.08, 0.9}
};

double socioeconomic_transitionMatrix[NUM_SOCIOECONOMIC_STATES][NUM_SOCIOECONOMIC_STATES] = {
    // From poverty
    {0.7, 0.2, 0.07, 0.02, 0.01},
    // From low income
    {0.07, 0.63, 0.2, 0.07, 0.03},
    // From medium-low income
    {0.02, 0.1, 0.63, 0.2, 0.05},
    // From medium-high income
    {0.02, 0.05, 0.15, 0.6, 0.18},
    // From high income
    {0.01, 0.03, 0.06, 0.2, 0.7}
    //     // From primary school
    // {0.7, 0.22, 0.07, 0.01, 0.0},
    // // From high school
    // {0.05, 0.74, 0.15, 0.05, 0.01},
    // // From diploma
    // {0.01, 0.05, 0.73, 0.15, 0.01},
    // // From undergraduate
    // {0.02, 0.03, 0.05, 0.80, 0.1},
    // // From postgraduate
    // {0.0, 0.00, 0.02, 0.08, 0.9}
};




// Define the offsets for a Moore neighborhood
int neighbor_offsets[NEIGHBORHOOD_SIZE][2] = {
    {-1, -1}, {0, -1}, {1, -1}, // Above
    {-1,  0},          {1,  0}, // Sides
    {-1,  1}, {0,  1}, {1,  1}  // Below
};

// Define the influence weights for the Moore neighborhood
double influenceWeights[NEIGHBORHOOD_SIZE] = {
    0.3, // Parent above
    0.2, // Parent to the left
    0.2, // Parent to the right
    0.3, // Parent below
    0.1, // Diagonal neighbor
    0.1, // Diagonal neighbor
    0.1, // Diagonal neighbor
    0.1  // Diagonal neighbor
};



//Functions
// Function to initialize the grid with individuals
void initializeGrid(Grid grid) {
    for (int i = 0; i < GRID_HEIGHT; ++i) {
        for (int j = 0; j < GRID_WIDTH; ++j) {
            grid[i][j].education_level = (int)(drand48() *  NUM_EDUCATION_STATES);
            grid[i][j].socioeconomic_status = (int)(drand48() * NUM_SOCIOECONOMIC_STATES);
            
            //grid[i][j].age = (int)(drand48() * 100);
            // Initialize other attributes...
        }
    }
}

//**CELLULAR AUTOMATA**//

// Function to evolve the grid using Cellular Automata rules
void evolveCellularAutomata(Grid grid,double NSE) {
    
    
    //updating based off of neighbours
    for (int x = 0; x < GRID_HEIGHT; ++x) {
        for (int y = 0; y < GRID_WIDTH; ++y) {
            new_blood(grid);
            updateCellBasedOnNeighbors(grid, x, y,NSE);
            education_makes_wealth(grid,x,y);
        }
    }
}

void new_blood(Grid grid){


}

// This function 


 //neighbour_strength_effect. Higher means weaker effect from neighbours
void updateCellBasedOnNeighbors(Grid grid, int x, int y,double NSE) {

    Individual *individual = &grid[x][y];
    double total_EduInfluence = 0.0;
    double total_EcoInfluence = 0.0;

    double avgNeighborEdu = 0.0;
    double avgNeighborEco = 0.0;
    int neighbour_count =0;

    int rich_parents = 0;
    // Loop over all neighbors
    for (int i = 0; i < NEIGHBORHOOD_SIZE; ++i) {
        int nx = x + neighbor_offsets[i][0];
        int ny = y + neighbor_offsets[i][1];

        // Check if neighbor coordinates are within the grid bounds
        if (nx >= 0 && nx < GRID_HEIGHT && ny >= 0 && ny < GRID_WIDTH) {
            Individual *neighbor = &grid[nx][ny];

            // Calculate influence based on the neighbor's state and the weight
            double ecoInfluence = neighbor->socioeconomic_status * influenceWeights[i];
            double eduInfluence = neighbor->education_level * influenceWeights[i];
            total_EcoInfluence += ecoInfluence;
            total_EduInfluence += eduInfluence;

            //neighbour migration criterion
            avgNeighborEco +=neighbor->socioeconomic_status;
            avgNeighborEdu +=neighbor->education_level;
            neighbour_count++;

             if ((nx == x || ny == y) && neighbor->socioeconomic_status >= 3) {
                rich_parents++;
            }
        }

    }
    //**MIGRATION**//
    
    avgNeighborEco /= neighbour_count;
    avgNeighborEdu /= neighbour_count;

    int swap = 0;

    //here we enforce the rule of migrating if neighbours are too poor or too dumb
    if ( individual->education_level > avgNeighborEdu + 2 && individual->socioeconomic_status > avgNeighborEco + 2  ) 
    {
        swap = migrate(grid,individual,x,y);

        if (swap != 0){return;}
        
    }

    //RICH PARENTS RULE//
    if(rich_parents>=2  &&  individual->education_level < 3){ //individual->age <= 24 &&
        
        if(drand48()<0.4){
            individual-> education_level++;
            //printf("rich daddy");
        }
    }


    // Update the individual's state based on the total influence
    
    double edu_influence_threshold = drand48()*(NUM_EDUCATION_STATES/2.0) + NUM_EDUCATION_STATES/2;
    double eco_influence_threshold = drand48()*(NUM_SOCIOECONOMIC_STATES/2.0) + NUM_SOCIOECONOMIC_STATES/2;
    
    //printf("edu:%f\n",edu_influence_threshold);
    //printf("eco:%f\n",eco_influence_threshold);
    
    int countEdu_up=0,countEdu_down = 0, countEco_up = 0, countEco_down=0;
    
   

//inner if statements are just to enforce minimums and maximums
    if (total_EduInfluence > (edu_influence_threshold+NSE)) 
    {
        if(individual->education_level < 4){
            individual->education_level +=1;
            countEdu_up++;
        }
    }
    else if(total_EduInfluence < (edu_influence_threshold-NSE))
    {
        if (individual->education_level > 0) {
            individual->education_level -=1;
            countEdu_down++;
        }
        
    }

    if (total_EcoInfluence > (eco_influence_threshold+NSE))
    {
         if (individual->socioeconomic_status < 4) {
            individual->socioeconomic_status +=1;
            countEco_up++;
        } 
        
    }
    else if(total_EcoInfluence < (eco_influence_threshold-NSE))
    {
        
        if (individual->socioeconomic_status > 0) {
            individual->socioeconomic_status -=1;
            countEco_down++;
        }
    }

totalGlobal->countEcoDown += countEco_down;
totalGlobal->countEcoUp += countEco_up;
totalGlobal->countEduDown += countEdu_down;
totalGlobal->countEduUp += countEdu_up;


//printf("Eco_Down:%d \n Eco_UP:%d \n Edu_Down:%d \n Edu_UP:%d",countEco_down,countEco_up,countEdu_down,countEdu_up);


}

//**METHODS USED FOR MIGRATION**//

// Calculate the average education level and economic status of an individual's neighbors
NeighborAverages calculate_neighbor_averages(Grid grid, int x, int y) {
    NeighborAverages averages;
    int sum_education = 0, sum_economic = 0, count = 0;
    
    // Loop through the neighbors
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue; // Skip the individual itself
            int neighbor_x = x + i;
            int neighbor_y = y + j;
            
            // Check bounds
            if (neighbor_x >= 0 && neighbor_x < GRID_WIDTH && neighbor_y >= 0 && neighbor_y < GRID_HEIGHT) {
                sum_education += grid[neighbor_x][neighbor_y].education_level;
                sum_economic += grid[neighbor_x][neighbor_y].socioeconomic_status;
                count++;
            }
        }
    }
    
    if (count > 0) {
        averages.avg_education = sum_education / count;
        averages.avg_economic = sum_economic / count;
    }
    
    return averages;
}




//swaps 2 individuals in the grid
void swap_individuals(Grid grid, int x1, int y1, int x2, int y2) {
    // Swap the pointers in the grid
    Individual temp = grid[x1][y1];
    grid[x1][y1] = grid[x2][y2];
    grid[x2][y2] = temp;
}



//method which orchestrates the migration of an individual 
int migrate(Grid grid,Individual *individual, int x, int y){

    NeighborAverages neighbours_res;

    int r_x = (int)(drand48()*GRID_WIDTH);
    int r_y = (int)(drand48()*GRID_HEIGHT);

   
    int match = 0;
    for(int i=r_x; i<GRID_WIDTH; i++){
        
        for (int j =r_y; j<GRID_HEIGHT; j++){
            
            neighbours_res = calculate_neighbor_averages(grid, i,j );

            if (abs(individual->education_level - neighbours_res.avg_education) <= 0.2 && abs(individual->socioeconomic_status - neighbours_res.avg_economic) <=0.2){

                swap_individuals(grid,x,y,i,j);
                match = 1;
                //printf("swapped: x= %d y= %d  with x = %d  y = %d",x,y,i,j);
            }

        }

    }

//only occurs if match not found in previous nested for loop
    if(match == 0 ){

        for (int i = 0; i<r_x; i++){

            for (int j =0; j<r_y; j++){
                neighbours_res = calculate_neighbor_averages(grid, i,j );

                if (abs(individual->education_level - neighbours_res.avg_education) <= 0.2 && abs(individual->socioeconomic_status - neighbours_res.avg_economic) <=0.2){

                    swap_individuals(grid,x,y,i,j);
                    match = 1;
                }

            }
        }

    }

    if(match ==0)
    {
        //both education and socioeconomic status
       // printf("no swap");
    }

    return match; 

}

/*Here we use the difference between the individual's education level and socioeconomic status to calculate the probability of the individual
transitioning to a higher socioeconomic state using an exponential probability distribution.

*/

    //probabilities for 1,2,3: 0.2591817793182821
    //0.4511883639059736
    //0.5934303402594008
double calculateUpgradeProbability(int educationLevel, int socioeconomicStatus) {
    int difference = educationLevel - socioeconomicStatus;

    // Parameters for the exponential function
    double base = 1.2; // Base of the exponential function
    double scaleFactor = 0.3; // Adjust this factor to control how steep the curve is



    // Calculate the probability
    double finalProbability = 1 - exp(-scaleFactor * difference);

    // Ensure the probability is within 0 and 1
    if (finalProbability > 1.0) finalProbability = 1.0;
    if (finalProbability < 0.0) finalProbability = 0.0;

    return finalProbability;
}


//
void education_makes_wealth(Grid grid, int x, int y){
    Individual *individual = &grid[x][y];

    // if(drand48() < individual->education_level){


    // }

    

    //if((individual->education_level - individual->socioeconomic_status) >=2){
        
        if(drand48() < calculateUpgradeProbability(individual->education_level,individual->socioeconomic_status) && individual->socioeconomic_status<(NUM_SOCIOECONOMIC_STATES-1)){
            
            individual->socioeconomic_status++;
        }
    
}









//**Markov Transitions**/


// Function to perform Markov Chain progression for each individual
void performMarkovChain(Grid grid) {

    for (int i = 0; i < GRID_HEIGHT; ++i) {
        
        for (int j = 0; j < GRID_WIDTH; ++j) {
        
            updateState(&grid[i][j]);
        
        }
    }

}


// Function to update the state of an individual based on the Markov Matrices
void updateState(Individual *individual) {

    // Update education level
    double r = (double)drand48(); // Random number between 0 and 1
    double cumulativeProbability = 0.0;

    for (int i = 0; i < NUM_EDUCATION_STATES; ++i) {

        cumulativeProbability += education_transitionMatrix[individual->education_level][i];

        if (r <= cumulativeProbability) {
            individual->education_level = i;
            break;
        }
    }

    // Update economic status
    r = (double)drand48(); // Random number between 0 and 1, generate again for economic status
    cumulativeProbability = 0.0;

    for (int i = 0; i < NUM_SOCIOECONOMIC_STATES; ++i) {

        cumulativeProbability += socioeconomic_transitionMatrix[individual->socioeconomic_status][i];

        if (r <= cumulativeProbability) {
            individual->socioeconomic_status = i;
            break;
        }
    }

}








//**MonteCarlo**//

// Function to introduce random events using Monte Carlo sampling
void performMonteCarloSampling(Grid grid,double prob_scholarship,double prob_injury) {

        for (int x = 0; x < GRID_HEIGHT; ++x) {
        for (int y = 0; y < GRID_WIDTH; ++y) {
            simulateLifeEvents(&grid[x][y], prob_scholarship, prob_injury);
        }
    }
    
}



void simulateLifeEvents(Individual *individual,double prob_scholarship,double prob_injury) {
    // Generate a random number between 0 and 1
    double randEvent = drand48();

    
    // Check for scholarship
    if (randEvent < prob_scholarship && individual->education_level!=4) {

        individual->education_level += SCHOLARSHIP_IMPACT;
        countScholar++;
    }
    
    // Ensure education level does not go out of bounds
    

    // Reset random number for the next event
    randEvent = drand48();

    // Check for injury
    if (randEvent < prob_injury && individual->socioeconomic_status != 0) {
        individual->socioeconomic_status += INJURY_IMPACT;
        countInjury++;
    }

  
}
 
















  









//**DATA PROCESSING**//


double calculateGiniCoefficient(Grid grid) {
    // Placeholder for Gini coefficient calculation
    return 0.0;
}

double calculateAverageEducationLevel(Grid grid) {

    double totalEducation = 0;
    int numIndividuals = GRID_HEIGHT * GRID_WIDTH;
    for (int i = 0; i < GRID_HEIGHT; i++) {
        for (int j = 0; j < GRID_WIDTH; j++) {
            totalEducation += grid[i][j].education_level;
        }
    }
    return totalEducation / numIndividuals;
}

double calculateAverageSocioeconomicStatus(Grid grid) {

    int totalEcoStat = 0;
    int numIndividuals = GRID_HEIGHT * GRID_WIDTH;
    for (int i = 0; i < GRID_HEIGHT; i++) {
        for (int j = 0; j < GRID_WIDTH; j++) {
            if (grid[i][j].socioeconomic_status != 0) {
                totalEcoStat += grid[i][j].socioeconomic_status;
            }
        }
    }
    return (double)totalEcoStat / numIndividuals;
}

double calculatePovertyRate(Grid grid) {
    int povertyCount = 0;
    int numIndividuals = GRID_HEIGHT * GRID_WIDTH;
    for (int i = 0; i < GRID_HEIGHT; i++) {
        for (int j = 0; j < GRID_WIDTH; j++) {
            if (grid[i][j].socioeconomic_status == 0) {
                povertyCount++;
            }
        }
    }
    return (double)povertyCount / numIndividuals;
}





// Function to calculate all metrics for a single timestep
Metrics calculateMetrics(Grid grid) {
    Metrics m;
    // Calculate each metric
    m.giniCoefficient = calculateGiniCoefficient(grid);
    m.averageEducationLevel = calculateAverageEducationLevel(grid);
    m.averageSocioeconomicStatus = calculateAverageSocioeconomicStatus(grid);
    m.povertyRate = calculatePovertyRate(grid);
    // ... calculate other metrics ...
    return m;
}


// Function to write metrics to a file
void writeMetricsToFile(Metrics* metricsArray, int numTimesteps, FILE* outputFile) {
    for (int i = 0; i < numTimesteps; ++i) {
        fprintf(outputFile, "%d,%f,%f,%f,%f\n",
                i,
                metricsArray[i].giniCoefficient,
                metricsArray[i].averageEducationLevel,
                metricsArray[i].averageSocioeconomicStatus,
                metricsArray[i].povertyRate);
        // ... write other metrics ...
    }
}


//ouputing entire gird - for education statee
void outputGridEducationState(Grid grid, int num, const char* baseFilename) {
    // Construct the filename with the timestep
    char filename[256];
    snprintf(filename, sizeof(filename), "%s_%d.csv", baseFilename, num);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // Write the grid state to the file
    for (int i = 0; i < GRID_HEIGHT; ++i) {
        for (int j = 0; j < GRID_WIDTH; ++j) {
            fprintf(file, "%d,", grid[i][j].education_level);
            // If you need to print more than one attribute, add them here
        }
        fprintf(file, "\n"); // End of row
    }

    fclose(file);
}

//ouputing entire gird - for education statee
void outputGridEcoState(Grid grid, int num, const char* baseFilename) {
    // Construct the filename with the timestep
    char filename[256];
    snprintf(filename, sizeof(filename), "%s_%d.csv", baseFilename, num);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    // Write the grid state to the file
    for (int i = 0; i < GRID_HEIGHT; ++i) {
        for (int j = 0; j < GRID_WIDTH; ++j) {
            fprintf(file, "%d,", grid[i][j].socioeconomic_status);
            // If you need to print more than one attribute, add them here
        }
        fprintf(file, "\n"); // End of row
    }

    fclose(file);
}


// // Sigmoid function for probability calculation
// double sigmoid(double x) {
//     return 1 / (1 + exp(-x));
// }



// Function to delete files matching a pattern
void deleteFilesMatchingPattern(const char *directory, const char *pattern) {
    DIR *dir;
    struct dirent *ent;
    char filepath[512];

    if ((dir = opendir(directory)) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (strstr(ent->d_name, pattern) != NULL) {
                snprintf(filepath, sizeof(filepath), "%s/%s", directory, ent->d_name);
                if (remove(filepath) != 0) {
                    perror("Error deleting file");
                }
            }
        }
        closedir(dir);
    } else {
        perror("Error opening directory");
    }
}


// This function runs the system for a certain number of equilibration steps and discards the data
void equilibrateSystem(Grid grid, double nse_value, int equilibration_steps) {
    for (int i = 0; i < equilibration_steps; ++i) {
        evolveCellularAutomata(grid, nse_value);
        performMarkovChain(grid);
        performMonteCarloSampling(grid, 0.01, 0.01); // Use actual probabilities needed
    }
}




// Main function
int main() {

    //Delete old files
    deleteFilesMatchingPattern(".", "grid_eco_state_");
    deleteFilesMatchingPattern(".", "grid_edu_state_");

    // Seed the random number generator
    srand48(time(NULL));

    // Initialize the grid and call it society
    Grid society;
    initializeGrid(society);

    double NSE = NSE_glo;
    
    int timesteps = TIME_STEPS;
    

    Metrics metricsArray[timesteps];
    Metrics metricsArray_Scholarships[timesteps];
    FILE *outputFile = fopen("output.csv", "w");
    fprintf(outputFile, "Timestep,GiniCoefficient,AverageEducationLevel,AverageSocioeconomicStatus,PovertyRate/Unemployed\n");

    UpdateCounts* total = malloc(sizeof(UpdateCounts));
    

    total->countEduUp = 0;
    total->countEduDown = 0;
    total->countEcoUp = 0;
    total->countEcoDown = 0;

    totalGlobal = total;
    // Main simulation loop

    
    for (int i = 0; i < timesteps; ++i) {
        // Evolve the grid using Cellular Automata rules
        evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
        performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
        performMonteCarloSampling(society,SCHOLARSHIP,INJURY);

        // Record data for analysis
        metricsArray[i] = calculateMetrics(society);    
        
        if(i%INTERVAL_OF_GRID_PLOTS ==0) {
             outputGridEducationState(society,i/INTERVAL_OF_GRID_PLOTS,"grid_edu_state");
             outputGridEcoState(society,i/INTERVAL_OF_GRID_PLOTS,"grid_eco_state");
         }
    }
    


    printf("Eco_Down:%d \n Eco_UP:%d \n Edu_Down:%d \n Edu_UP:%d",total->countEcoDown,total->countEcoUp,total->countEduDown,total->countEduUp);
    writeMetricsToFile(metricsArray, timesteps, outputFile); 
    fclose(outputFile);
    printf("Injuries: %d \t Scholarships: %d",countInjury,countScholar);
    free(totalGlobal);
    


    //**Scholarship Variable**//

    
    // Open a file for writing
    FILE *file1 = fopen("schol_avgs.csv", "w");
    if (file1 == NULL) {
        perror("Error opening file");
        return -1;
    }

   
    fprintf(file1, "Scholarship Probability, Average Education, Average Socioeconomic Status\n");

    double sum_edu = 0;
    double sum_eco = 0;

    double avgEdu_schol [100];
    double avgEco_schol [100];

    for (double k =0; k<0.3; k+=0.02)
    {   
        initializeGrid(society);
        sum_edu = 0;
        sum_eco = 0;
        for (int i = 0; i < timesteps; ++i) {
            
             // Evolve the grid using Cellular Automata rules
            evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
            performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
            performMonteCarloSampling(society,k,0.005);

        // Record data for analysis
            metricsArray_Scholarships[i] = calculateMetrics(society);  
            sum_edu += metricsArray_Scholarships[i].averageEducationLevel;
            sum_eco += metricsArray_Scholarships[i].averageSocioeconomicStatus;
        }
        avgEdu_schol[(int)(k*5)] = sum_edu/((double)timesteps);
        avgEco_schol[(int)(k*5)] = sum_eco/((double)timesteps);
        fprintf(file1, "%lf, %lf, %lf\n", k, avgEdu_schol[(int)(k*5)], avgEco_schol[(int)(k*5)]);
    }
    fclose(file1);


//**INJURY VARIABLE**//

FILE *file2 = fopen("injury_avgs.csv", "w");
if (file2 == NULL) {
    perror("Error opening file");
    return -1;
}

fprintf(file2, "Injury Probability, Average Education, Average Socioeconomic Status\n");

 sum_edu = 0;
 sum_eco = 0;

double avgEdu_injury[100];
double avgEco_injury[100];

// Loop for different injury probabilities
for (double k = 0; k < 0.3; k += 0.02)
{   
    initializeGrid(society);
    sum_edu = 0;
    sum_eco = 0;
    for (int i = 0; i < timesteps; ++i) {
        
        // Evolve the grid using Cellular Automata rules
        evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
        performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
        // Adjust the performMonteCarloSampling function to use the injury probability 'k'
        performMonteCarloSampling(society, 0.005, k); // Assuming the third parameter is the injury probability

        // Record data for analysis
        Metrics metrics = calculateMetrics(society);  
        sum_edu += metrics.averageEducationLevel;
        sum_eco += metrics.averageSocioeconomicStatus;
    }
    avgEdu_injury[(int)(k * 5)] = sum_edu / ((double)timesteps);
    avgEco_injury[(int)(k * 5)] = sum_eco / ((double)timesteps);
    fprintf(file2, "%lf, %lf, %lf\n", k, avgEdu_injury[(int)(k * 5)], avgEco_injury[(int)(k * 5)]);
}
fclose(file2);



//**NSE VARIABLE**//
    //File for NSE averages
    
    FILE *file_nse = fopen("nse_avgs.csv", "w");
    if (file_nse == NULL) {
        perror("Error opening file");
        return -1;
    }

    fprintf(file_nse, "NSE Value, Average Education, Average Socioeconomic Status,Avg Edu Variance,Avg Eco Variance\n");

    // Loop for different NSE values
    int equilibration = 1000;

    for (double nse_value = 0.5; nse_value <= 1.4; nse_value += 0.1) {
        // Set the global NSE value
        

        initializeGrid(society); // Reinitialize the grid for each NSE value

        double sum_edu = 0.0;
        double sum_edu_squared = 0.0;
        double sum_eco = 0.0;
        double sum_eco_squared = 0.0;

        equilibrateSystem(society, nse_value, equilibration);

        for (int i = equilibration; i < timesteps; ++i) {
            evolveCellularAutomata(society,nse_value);
            performMarkovChain(society);
            performMonteCarloSampling(society, 0.01, 0.01); // Example probabilities

            Metrics metrics = calculateMetrics(society);
            sum_edu += metrics.averageEducationLevel;
            sum_edu_squared += metrics.averageEducationLevel * metrics.averageEducationLevel;
            sum_eco += metrics.averageSocioeconomicStatus;
            sum_eco_squared += metrics.averageSocioeconomicStatus * metrics.averageSocioeconomicStatus;
        }

        double avgEdu_nse = sum_edu / (double)(timesteps-equilibration);
        double avgEco_nse = sum_eco / (double)(timesteps-equilibration);
        double varEdu = (sum_edu_squared / (double)timesteps-equilibration) - (avgEdu_nse * avgEdu_nse);
        double varEco = (sum_eco_squared / (double)timesteps-equilibration) - (avgEco_nse * avgEco_nse);

        fprintf(file_nse, "%lf, %lf, %lf, %lf, %lf\n", nse_value, avgEdu_nse, avgEco_nse,varEdu,varEco);
    }

    //
    for (double nse_value = 1.3; nse_value <= 1.65; nse_value += 0.01) {
        // Set the global NSE value
        

        initializeGrid(society); // Reinitialize the grid for each NSE value

        double sum_edu = 0.0;
        double sum_edu_squared = 0.0;
        double sum_eco = 0.0;
        double sum_eco_squared = 0.0;
        
        equilibrateSystem(society, nse_value, equilibration);

        for (int i = equilibration; i < timesteps; ++i) {
            evolveCellularAutomata(society,nse_value);
            performMarkovChain(society);
            performMonteCarloSampling(society, 0.01, 0.01); // Example probabilities

            Metrics metrics = calculateMetrics(society);
            sum_edu += metrics.averageEducationLevel;
            sum_edu_squared += metrics.averageEducationLevel * metrics.averageEducationLevel;
            sum_eco += metrics.averageSocioeconomicStatus;
            sum_eco_squared += metrics.averageSocioeconomicStatus * metrics.averageSocioeconomicStatus;
        }

        double avgEdu_nse = sum_edu / (double)(timesteps-equilibration);
        double avgEco_nse = sum_eco / (double)(timesteps-equilibration);
        double varEdu = (sum_edu_squared / (double)(timesteps-equilibration)) - (avgEdu_nse * avgEdu_nse);
        double varEco = (sum_eco_squared / (double)(timesteps-equilibration)) - (avgEco_nse * avgEco_nse);

        fprintf(file_nse, "%lf, %lf, %lf, %lf, %lf\n", nse_value, avgEdu_nse, avgEco_nse,varEdu,varEco);
    }


        for (double nse_value = 1.65; nse_value <= 2.0; nse_value += 0.1) {
        // Set the global NSE value
        

        initializeGrid(society); // Reinitialize the grid for each NSE value

        double sum_edu = 0.0;
        double sum_edu_squared = 0.0;
        double sum_eco = 0.0;
        double sum_eco_squared = 0.0;
        
        equilibrateSystem(society, nse_value, equilibration);

            for (int i =equilibration; i < timesteps; ++i) {

                evolveCellularAutomata(society,nse_value);
                performMarkovChain(society);
                performMonteCarloSampling(society, 0.01, 0.01); // Example probabilities

                Metrics metrics = calculateMetrics(society);
                sum_edu += metrics.averageEducationLevel;
                sum_edu_squared += metrics.averageEducationLevel * metrics.averageEducationLevel;
                sum_eco += metrics.averageSocioeconomicStatus;
                sum_eco_squared += metrics.averageSocioeconomicStatus * metrics.averageSocioeconomicStatus;
            }

        double avgEdu_nse = sum_edu / (double)(timesteps-equilibration);
        double avgEco_nse = sum_eco / (double)(timesteps-equilibration);
        double varEdu = (sum_edu_squared / (double)((timesteps-equilibration))) - (avgEdu_nse * avgEdu_nse);
        double varEco = (sum_eco_squared / (double)(timesteps-equilibration)) - (avgEco_nse * avgEco_nse);

        fprintf(file_nse, "%lf, %lf, %lf, %lf, %lf\n", nse_value, avgEdu_nse, avgEco_nse,varEdu,varEco);
    }

    fclose(file_nse);


    //**ECONOMIC DEPRESSION FOR  timesteps/10 in the middle **//
    FILE *file_depression = fopen("depression.csv", "w");
    if (file_depression == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(file_depression, "Timestep,GiniCoefficient,AverageEducationLevel,AverageSocioeconomicStatus,PovertyRate/Unemployed\n");
    for (int i = 0; i < timesteps; ++i) {
            
             // Evolve the grid using Cellular Automata rules
            evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
            performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
            performMonteCarloSampling(society,SCHOLARSHIP,INJURY);

        // Record data for analysis
            metricsArray[i] = calculateMetrics(society);  
        }

    for (int i = timesteps/2; i < timesteps/2 + timesteps/10; ++i) {
            
             // Evolve the grid using Cellular Automata rules
            evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
            performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
            performMonteCarloSampling(society,SCHOLARSHIP,0.2);

        // Record data for analysis
        metricsArray[i] = calculateMetrics(society);  
        }

    for (int i = timesteps/2 + timesteps/10; i < timesteps; ++i) {
            
             // Evolve the grid using Cellular Automata rules
            evolveCellularAutomata(society,NSE_glo);

        // Perform Markov Chain progression for each individual
            performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
            performMonteCarloSampling(society,SCHOLARSHIP,INJURY);

        // Record data for analysis
            metricsArray[i] = calculateMetrics(society);  
        }
        
        writeMetricsToFile(metricsArray, timesteps, file_depression); 
        fclose(file_depression);

    return 0;
}



















