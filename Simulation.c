/*
1)Consider how social dynamics change over time. For example, as individuals age, the influence of "parent" cells might decrease.

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Define constants for the grid size, number of states, etc.
#define GRID_WIDTH 100
#define GRID_HEIGHT 100
#define NUM_EDUCATION_STATES 5
#define NUM_SOCIOECONOMIC_STATES 5
#define NUM_EVENTS 3

#define NUM_ITERATIONS 100

#define POVERTY_THRESHOLD 1058
#define MINIMUM_WAGE 4067


#define NEIGHBORHOOD_SIZE 8


// Define a structure for an individual in the grid
typedef struct {
    int education_level;
    int socioeconomic_status;
    int age;
    // Add more attributes as needed
} Individual;

// Define a grid of individuals.
typedef Individual Grid[GRID_HEIGHT][GRID_WIDTH];




// Function prototypes
void initializeGrid(Grid grid);
void evolveCellularAutomata(Grid grid);
void performMarkovChain(Grid grid);
void performMonteCarloSampling(Grid grid);
void recordData(Grid grid);
void analyzeData();


// Define a struct to hold all metrics for a single timestep
typedef struct {
    double giniCoefficient;
    double averageEducationLevel;
    double employmentRate;
    double povertyRate;
    // ... other metrics ...
} Metrics;

// Function prototypes for metric calculations
Metrics calculateMetrics(Grid grid);
void writeMetricsToFile(Metrics* metricsArray, int numTimesteps, FILE* outputFile);



// Transition matrix where transitionMatrix[i][j] is the probability of moving from state i to state j
double education_transitionMatrix[NUM_EDUCATION_STATES][NUM_EDUCATION_STATES] = {
    // These values should be filled in based on your data or model
    {0.1, 0.6, 0.1, 0.1, 0.1},
    {0.1, 0.1, 0.6, 0.1, 0.1},
    {},
    {},
    {}

    // ... and so on for each state
};

double socioeconomic_transitionMatrix[NUM_SOCIOECONOMIC_STATES][NUM_SOCIOECONOMIC_STATES] = {
    // These values should be filled in based on your data or model
    {0.1, 0.6, 0.1, 0.1, 0.1},
    {0.1, 0.1, 0.6, 0.1, 0.1},
    {},
    {},
    {}

    // ... and so on for each state
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


// Main function
int main() {
    // Seed the random number generator
    srand48(time(NULL));

    // Initialize the grid and call it society
    Grid society;
    initializeGrid(society);

    // Main simulation loop
    for (int timestep = 0; timestep < 1000; ++timestep) {
        // Evolve the grid using Cellular Automata rules
        evolveCellularAutomata(society);

        // Perform Markov Chain progression for each individual
        performMarkovChain(society);

        // Introduce random events using Monte Carlo sampling
        performMonteCarloSampling(society);

        // Record data for analysis
        recordData(society);
    }

    // Analyze the recorded data
    analyzeData();

    return 0;
}

// Function to initialize the grid with individuals
void initializeGrid(Grid grid) {
    for (int i = 0; i < GRID_HEIGHT; ++i) {
        for (int j = 0; j < GRID_WIDTH; ++j) {
            grid[i][j].education_level = (int)(drand48() *  NUM_EDUCATION_STATES);
            grid[i][j].socioeconomic_status = (int)(drand48() * NUM_SOCIOECONOMIC_STATES);
            grid[i][j].age = (int)(drand48() * 100);
            // Initialize other attributes...
        }
    }
}

// Function to evolve the grid using Cellular Automata rules
void evolveCellularAutomata(Grid grid) {
    // Implement the evolution rules for the Cellular Automata


    //updating based off of neighbours
    for (int x = 0; x < GRID_HEIGHT; ++x) {
        for (int y = 0; y < GRID_WIDTH; ++y) {
            updateCellBasedOnNeighbors(grid, x, y);
        }
    }
}


// Function to perform Markov Chain progression for each individual
void performMarkovChain(Grid grid) {

    for (int i = 0; i < GRID_HEIGHT; ++i) {
        
        for (int j = 0; j < GRID_WIDTH; ++j) {
        
            updateState(&grid[i][j]);
        
        }
    }

}

// Function to introduce random events using Monte Carlo sampling
void performMonteCarloSampling(Grid grid) {
    // Implement the Monte Carlo sampling for random events
}

// Function to record data for analysis
void recordData(Grid grid) {

//  


}

// Function to analyze the recorded data
void analyzeData() {
    // Implement data analysis (this might be done in Python as per the proposal)
}
 





//HELPER FUNCTIONS



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






// Function to calculate the influence of neighbors
void updateCellBasedOnNeighbors(Grid grid, int x, int y) {

    Individual *individual = &grid[x][y];
    double total_EduInfluence = 0.0;
    double total_EcoInfluence = 0.0;
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

        }
    }

    // Update the individual's state based on the total influence
    // This is a placeholder for your update logic
    double edu_influence_threshold = drand48()*(NUM_EDUCATION_STATES/2) + NUM_EDUCATION_STATES/2;
    double eco_influence_threshold = drand48()*(NUM_SOCIOECONOMIC_STATES/2) + NUM_SOCIOECONOMIC_STATES/2;


    if (total_EduInfluence > edu_influence_threshold+1) 
    {
        individual->education_level +=1;
    }
    else if(total_EduInfluence < edu_influence_threshold-1)
    {
        individual->education_level -=1;

    }

    if (total_EcoInfluence > eco_influence_threshold+1)
    {
        individual->socioeconomic_status +=1;
    }
    else if(total_EcoInfluence < eco_influence_threshold-1)
    {
        individual->socioeconomic_status -=1;

    }

  
}







// Function to calculate all metrics for a single timestep
Metrics calculateMetrics(Grid grid) {
    Metrics m;
    // Calculate each metric
    m.giniCoefficient = calculateGiniCoefficient(grid);
    m.averageEducationLevel = calculateAverageEducationLevel(grid);
    m.employmentRate = calculateEmploymentRate(grid);
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
                metricsArray[i].employmentRate,
                metricsArray[i].povertyRate);
        // ... write other metrics ...
    }
}




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

double calculateEmploymentRate(Grid grid) {
    int employedCount = 0;
    int numIndividuals = GRID_HEIGHT * GRID_WIDTH;
    for (int i = 0; i < GRID_HEIGHT; i++) {
        for (int j = 0; j < GRID_WIDTH; j++) {
            if (grid[i][j].socioeconomic_status == EMPLOYED) {
                employedCount++;
            }
        }
    }
    return (double)employedCount / numIndividuals;
}

double calculatePovertyRate(Grid grid) {
    int povertyCount = 0;
    int numIndividuals = GRID_HEIGHT * GRID_WIDTH;
    for (int i = 0; i < GRID_HEIGHT; i++) {
        for (int j = 0; j < GRID_WIDTH; j++) {
            if (grid[i][j].income < POVERTY_THRESHOLD) {
                povertyCount++;
            }
        }
    }
    return (double)povertyCount / numIndividuals;
}



// // Sigmoid function for probability calculation
// double sigmoid(double x) {
//     return 1 / (1 + exp(-x));
// }