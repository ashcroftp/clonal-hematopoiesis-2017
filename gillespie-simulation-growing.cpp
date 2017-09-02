// gillespie-simulation-growing.cpp
// Author: Peter Ashcroft, ETH Zürich 

// Stochastic simulation algorithm (SSA) [Gillespie, J Stat Phys, (1977)]
// for invetsigating clonal expansion of hematopoietic stem cells
// which exist in two anatomical compartments.
// We record time to reach 4% clonality, and simulate the expansion of
// the hematopoietic system during maturation 
// Designed for use on a HPC cluster

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

using namespace std;

// Class to store all parameters for the model
// and the Gillespie algorithm
class GillespieAlgorithm
{
private:
  // System parameters
  unsigned nReactions;       // Number of reaction channels
  unsigned N;                // System size parameter (carrying capacity for total # niches)
  unsigned sStar,nStar;      // Neutral equilibrium parameters
  unsigned S2;               // Initial dose of donor cells
  double lt;                 // Lifetime of a cell in PB

  // Reaction rates
  double beta,delta,a,d,alpha;                    // Neutral reaction parameters
  double gBeta,gDelta,gA,gD,epsilon;              // Selective reaction parameters
  double beta1,beta2,delta1,delta2,a1,a2,d1,d2;   // Full reaction parameters
  double rNiche;                                  // Niche expansion rate
  
    
  // Analysis paramters
  unsigned paramIndex,simIndex;     // Parameter regime and simulation indices
  unsigned nRuns;                   // Number of times to run the algorithm
  vector<int> seeds;                // Vector to store seeds used for each run
  ofstream Output;                  // Ofstream for results
  char Filename[100];               // Results filename

public:
  // Constructor (Empty)
  GillespieAlgorithm()
    : paramIndex(0), simIndex(0), nRuns(1), nReactions(11)
  {  
    seeds.resize(nRuns);                                                    // Resize seed vector
    sprintf(Filename, "output_%02i_%03i.dat", paramIndex, simIndex);        // Define filename
  };

  // Constructor
  GillespieAlgorithm(unsigned& paramIndex_, unsigned& simIndex_, unsigned& nRuns_)
    : paramIndex(paramIndex_), simIndex(simIndex_), nRuns(nRuns_), nReactions(11)
  {
    seeds.resize(nRuns);                                                    // Resize seed vector
    sprintf(Filename, "output_%02i_%03i.dat", paramIndex, simIndex);          // Define filename
  };

  
  // Member functions
  void initialise_seeds();                                                  // Seed initialisation
  void set_parameters();                                                    // Assign parameters
  void output_parameters();                                                 // Output parameters
  void run();                                                               // Run the nRuns simulations
  double propensity_functions(unsigned& reaction, vector<unsigned>& x);     // Propensity function declaration
  void population_update(unsigned& reaction, vector<unsigned>& x);          // Population update declaration
  void compute(unsigned runIndex);                                          // Execute Gillespie algorithm and output sample data
  
}; // End of GillespieAlgorithm class



// Seed initialisation
void GillespieAlgorithm::initialise_seeds()
{
  random_device rd;                                                         // Define random device
  mt19937 my_seed_gen(rd());                                                // Use local machine randomness to initialise the seed generator
  uniform_int_distribution<int> seed_dist;                                  // Assign distribution of seeds (uniform ints)
  for(unsigned i = 0;i<1000;i++) int a = seed_dist(my_seed_gen);            // Burn-in 1000 random numbers
  for(unsigned i = 0; i < nRuns; i++) seeds[i] = seed_dist(my_seed_gen);    // Fill seeds vector
}// End of GillespieAlgorithm::initialise_seeds()



// Assign model parameters
void GillespieAlgorithm::set_parameters()
{
  N = 10000 * pow(10,paramIndex);           // System parameters
  sStar = floor(0.01 * (double)N);
  nStar = 9900;
  lt = 60.0/1440.0;
  S2 = 1;                                   // Here we use paramIndex to vary the dose of donor cells
  alpha = 0.0;                              // No death allowed in the bone marrow niches

  rNiche = 0.3 / 365.25;
  
  
  beta = 1.0/(7.0 * 40.0);                          // Neutral parameters
  delta = beta * (double)nStar / ( (double)sStar + alpha * (double)nStar );
  d = ( ( ( 1.0 / lt ) * (double)sStar ) / (double)nStar ) - beta;
  a = ( ( 1.0 / lt ) - ( ( beta * (double)nStar ) / ( (double)sStar + alpha * (double)nStar) ) ) * ( (double)N / ( (double)N - (double)nStar ) );

  gBeta = 1.0;                              // Selective parameters: Only change reproductive rate beta
  gA = 0.0;
  gDelta = 0.0;
  gD = 0.0;
  epsilon = 0.0;                            // epsilon = 0 gives the neutral model

  beta1 = beta;                             // Reaction parameters: Host
  a1 = a;
  delta1 = delta;
  d1 = d;

  beta2 = beta * (1.0 + epsilon * gBeta);   // Reaction parameters: Donor
  a2 = a * (1.0 + epsilon * gA);
  delta2 = delta * (1.0 + epsilon * gDelta);
  d2 = d * (1.0+ epsilon * gD);
}// End of GillespieAlgorithm::set_parameters()


// Run simulations
void GillespieAlgorithm::run()
{
  for(unsigned i=0; i<nRuns; i++) compute(i);
}// End of GillespieAlgorithm::run()



// Output model parameters
void GillespieAlgorithm::output_parameters()
{
  Output.open(Filename);                                                                              // Open output file
  Output << N << "\t" << sStar << "\t" << nStar << "\t" << lt << "\t" << S2 << "\t" << alpha << endl; // System size et al
  Output << beta << "\t" << delta << "\t" << a << "\t" << d << endl;                                  // Neutral parameters
  Output << gBeta << "\t" << gDelta << "\t" << gA << "\t" << gD << "\t" << epsilon << endl;           // Selective parameters
  Output.close();                                                                                     // Close output file
}// End of GillespieAlgorithm::output_parameters()



// Propensity functions for all possible reactions
double GillespieAlgorithm::propensity_functions(unsigned& reaction, vector<unsigned>& x)
{ // x=(s1,s2,n1,n2,N)
  switch(reaction)
    {
    case 0 : // n1 -> n1+s1
      return beta1 * (double)x[2];
      break;
     
    case 1 : // n2 -> n2+s2
      return beta2 * (double)x[3];
      break;

    case 2 : // s1 -> 0
      return delta1 * (double)x[0];
      break;

    case 3 : // s2 -> 0
      return delta2 * (double)x[1];
      break;
      
    case 4 : // n1 -> s1
      return d1 * (double)x[2];
      break;

    case 5 : // n2 -> s2
      return d2 * (double)x[3];
      break;

    case 6 : // s1+(N-n1-n2) -> n1
      return a1 * (double)x[0] * ( ( (double)x[4] - (double)x[2] - (double)x[3] ) / (double)x[4] );
      break;

    case 7 : // s2+(N-n1-n2) -> n2
      return a2 * (double)x[1] * ( ( (double)x[4] - (double)x[2] - (double)x[3] ) / (double)x[4] );
      break;

    case 8 : // n1 -> 0
      return alpha * delta1 * (double)x[2];
      break;

    case 9 : // n2 -> 0
      return alpha * delta2 * (double)x[3];
      break;

    case 10 : // N -> N+1
      return rNiche * (double)x[4] * ( 1.0 - (double)x[4]/(double)N );
      break;
      
    default :
      return 0.0;
    }// End of switch
}// End of GillespieAlgorithm:propensity_functions()



// Population update for all possible reactions
void GillespieAlgorithm::population_update(unsigned& reaction, vector<unsigned>& x)
{ // x=(s1,s2,n1,n2,N)
  switch(reaction)
    {   
    case 0 : // n1 -> n1+s1
      x[0]++;
      break;
	   
    case 1 : // n2 -> n2+s2
      x[1]++;
      break;
      
    case 2 : // s1 -> 0
      x[0]--;
      break;
      
    case 3 : // s2 -> 0
      x[1]--;
      break;
      
    case 4 : // n1 -> s1
      x[0]++;x[2]--;
      break;
      
    case 5 : // n2 -> s2
      x[1]++;x[3]--;
      break;
      
    case 6 : // s1+(N-n1-n2) -> n1
      x[0]--;x[2]++;
      break;
      
    case 7 : // s2+(N-n1-n2) -> n2
      x[1]--;x[3]++;
      break;

    case 8 : // n1 -> 0
      x[2]--;
      break;
      
    case 9 : // n2 -> 0
      x[3]--;
      break;

    case 10 : // N -> N+1
      x[4]++;
      break;
      
    default :
      x[0]++;x[0]--;
    }// End of switch
}// End of GillespieAlgorithm::population_update()



// Execute Gillespie algorithm and output data
void GillespieAlgorithm::compute(unsigned runIndex)
{
  // Initialisation
  int seed = seeds[runIndex];                         // Get seed from vector
  mt19937 mt_rand(seed);                              // Initialise RNG using seed
  uniform_real_distribution<double> dist(1.0e-9,1.0); // We want RNs in (0,1]
  
  vector<double> aa(nReactions);                      // Propensity value storage for each reaction
  double aa0;                                         // Propensity sum
 
  long double t;                                      // Time storage (long double because long time...)
  vector<unsigned> x(5);                              // Population storage vector: x=(s1,s2,n1,n2,N)
  unsigned y;                                         // y=min[(s1+n1),(s2+n2)]: used to test for fixation/extinction

  vector<double> r(2);                                // Variables to store random numbers, time step and fired reaction
  double tau,sum;
  unsigned mu;

  double chimerism;                                   // BM chimerism of donor cells
  
  // Initial condition
  t = 0.0;
  x[4] = floor( 0.05 * (double)N); // Start at 1/20th of the size
  x[0] = floor( (double)sStar * (double)x[4] / (double)N ); x[1] = S2;
  x[2] = floor( (double)nStar * (double)x[4] / (double)N ); x[3] = 0;
  y = min( x[0]+x[2] , x[1]+x[3] );
  chimerism = (double)x[3] / ((double)x[2] + (double)x[3]);
  
  while( y > 0 && chimerism < 0.04)                                      // Loop over timesteps until fixation
    {      
      r[0] = dist(mt_rand); r[1] = dist(mt_rand);     // Generate 2 uniform random numbers in (0,1]
      
      aa0 = 0.0;                                      // Calculate propensity functions
      for(unsigned m=0; m < nReactions; m++)
	{
	  aa[m] = propensity_functions(m,x);
	  aa0 += aa[m];
	}
      tau = (1.0/aa0) * log(1.0/r[0]);                // Calculate time step

      sum = 0.0;                                      // Determine which reaction channel has fired
      for(unsigned m=0; m < nReactions; m++)
	{
	  sum += aa[m];
	  if(sum > aa0*r[1]){mu=m;break;}
	}

     population_update(mu,x); t += tau;               // Update population and time

     chimerism = (double)x[3] / ((double)x[2] + (double)x[3]);
     
     y = min( x[0]+x[2] , x[1]+x[3] );
    }// End of loop over timesteps


  // Output result
  Output.open(Filename,ios::app);                     // Append to file after each run
  // Choice of output
  // Output system state at fixation
  Output << t << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << x[4] << endl;
  Output.close(); // Close output file

}// End of GillespieAlgorithm::compute()



// Executable takes arguments in array argv
int main(int argc, char* argv[])
{
  if(argc != 4)    // Argument check
    {
      cout << "ERROR: There should be three arguments!" << endl;
      return 1;
    }
  
  unsigned paramIndex = atoi(argv[1]);     // Index of parameters
  unsigned simIndex = atoi(argv[2]);       // Index of which block of simulations (we will split 1000 sims into groups of 10-100)
  unsigned nRuns = atoi(argv[3]);          // Number of runs per node
  
  GillespieAlgorithm gillespie(paramIndex,simIndex,nRuns);      // Construct the GillespieAlgorithm class with these indices
  gillespie.initialise_seeds();                                 // Initialise seeds vector using local machine randomness
  gillespie.set_parameters();                                   // Assign the model parameters according to paramIndex
  gillespie.output_parameters();                                // Output model parameters to file
  gillespie.run();                                              // Execute this simulation block --
                                                                // results are written to file in real time
  return 0;
}
