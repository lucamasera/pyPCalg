/**
	Author: Luca Masera.
    Edited from the PC++ project https://bitbucket.org/francesco-asnicar/gene_network_expansion.
	Copyright (C) 2017, all rights reserved

	This file (pypcalg.cpp) is part of the PC++ project.

	PyPCalg is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <boost/version.hpp>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <iostream>

namespace p = boost::python;
using namespace std;

#define M_SQRT1_2 0.70710678118654757273731092936941422522068023681640625 // definition for M_SQRT1_2

// Type definition for the pair of integers.
typedef std::pair<int,int> intpair;

class Graph {
public:
    Graph(const int, const int, const bool); //constructor : initialize the matrix and the numNeighbours
    ~Graph(void); //destructor : free the allocated memory
    bool** matrix; //matrix that represent the graph
    bool** cutMap; //matrix that represent the graph
    bool** double_directed; //keep track of double directed arcs
    int nRows; //represent the number of rows of matrix (and columns), bioData, means, standardDeviations and numNeighbours
    int nCols; //represent the number of columns of bioData
    double** bioData; //matrix that will contains the data to compute Pearson coefficient for the d-separation test
    double* means; //array that contains the means for each node in the graph
    std::string* probeIDs; //array that contains the name of each probe taken in account
    double* standardDeviations; //array that contains the standard deviations for each node in the graph
    int* numNeighbours; //represents the number of adjacents for each node (thought as column vector)
    double** rho; //represents the correlation matrix
    int*** sepSet; // contains the separation set for each pairs of nodes
    int** lenSepSet; // contains the lenght of the separation set for each pairs of nodes
    void computeStandardDeviations(void); //compute the standard deviation for each node in the graph
    void computeCorrelations(void); //compute the correlation coefficient of the base case, and store it in rho
    void initializeCutMap(); //initialize the boolean matrix to 'true', but the diagonal, setted to 'false'

private:
    void initializeMatrix(bool**, const int); //initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
    void initializeNeighbours(int*, const int); //initialize the array numNeighbours with the value dim-1, since the initial graph is connected
    void initializeZero(double*, const int); //initialize the given array till dim to 0.0
    void initializeMatrixZero(int**, const int); //initialize the int matrix to 0
    void initializeMatrixNull(int***, const int); //Initialize the int 3D matrix to NULL
};

// Compares two pairs of the type <probe identifier, lookup index>.
bool comparator (const intpair &, const intpair &);

// Reads the file TILE modifying both sizes and data.
bool readTile(const std::string, int* &, intpair** &, int &);

// Reads the file CGN saving the biological data that will be used to compute the correlation coefficients.
bool readCGN(const std::string, const intpair*, Graph* &);

// Computes the continous density function.
double comulativeNormalDistribution(const double);

// Finds the correlation coefficient when l is greater than 1 (formally, when l > 1).
double correlations(const int, const int, const Graph*, const int*, const int*, const int, double **);

// Checks if a given string (of the form array of chars) whether representing a float number or not.
bool isFloat(const char*);

// Prints the uncutted edges of the graph in a .csv file.
void fprintEdgesCSV(Graph*, const intpair*, const std::string, const int);

// Counts the number of (uncutted) edges in the graph.
int countArcs(bool**, const int, const int, const int);

// 
void testAndRemove(const int*, const int*, double, Graph* &, const int, const int, const int, const double, const bool);

// 
void remove(Graph* &);

// 
double getCorrelationCoefficient(const int*, const int*, const int, Graph* &, const int, const int, double**);

// 
void iterativeComb(int*, const int, const int, Graph* &, const int, const int, const double, double**, int*, const bool, const bool);

// 
void findAllSubsets(Graph* &, const int, const int, const int, const double, int*, double**, int*, const bool, const bool);

//
void skeleton(Graph* &, const double, const bool, const bool);

/**int main(int argc, char **argv)
{
  Py_Initialize();
}**/

p::tuple skeleton_wrapper(const p::list& expression_data, const float alpha, const bool return_sepset){
  int n_rows = boost::python::len(expression_data);
  int n_cols = boost::python::len(p::extract<p::list>((expression_data)[0]));

  Graph* g = new Graph(n_rows, n_cols, return_sepset);

  g->bioData = new double*[n_rows];

  for (int i = 0; i < g->nRows; i++) {
    g->bioData[i] = new double[n_cols];
  }

  for(int i = 0; i < n_rows; i++){
    p::list row = p::extract<p::list>((expression_data)[i]);
    for (int j = 0; j < n_cols; ++j){
      float x = p::extract<float>( (row)[j] );
      g->bioData[i][j] = x;
    }
  }

  //compute the standard deviations
  g->computeStandardDeviations();
    
  //compute the correlations coefficients
  g->computeCorrelations();

  skeleton(g, alpha, false, return_sepset);

  p::tuple retval;

  p::list adj;
  for(int i = 0; i < n_rows; i++){
    p::list row;
    for (int j = 0; j < n_rows; ++j){
      row.append(g->matrix[i][j]);
    }
    adj.append(row);
  }
  
  if (return_sepset) {
    p::list sepsets;
    for(int i = 1; i < n_rows; i++){
      p::list row;
      for (int j = 0; j < i; ++j){
        p::list sepset;
        for (int k = 0; k < g->lenSepSet[j][i]; k++){
          sepset.append(g->sepSet[j][i][k]);
        }
        row.append(sepset);
      }
      sepsets.append(row);
    }

    retval = p::make_tuple(adj, sepsets);
  } else {
    retval = p::make_tuple(adj, p::object());
  }

  delete g;

  return retval;
}

BOOST_PYTHON_MODULE(pypcalg)
{
    p::def("skeleton", skeleton_wrapper);
}

/**
 *  Construction that take as parameter the dimension with which the matrix is build
 */
Graph::Graph(const int n_rows, const int n_cols, const bool direct) {
  nRows = n_rows;
  nCols = n_cols;

  numNeighbours = new int[nRows];

  //create the bool matrix
  matrix = new bool*[nRows];

  for (int i = 0; i < nRows; i++) {
    matrix[i] = new bool[nRows];
  }

  cutMap = new bool*[nRows];

  for (int i = 0; i < nRows; i++) {
    cutMap[i] = new bool[nRows];
  }

  means = new double[nRows];
  probeIDs = new std::string[nRows];
  standardDeviations = new double[nRows];

  //create rho matrix
  rho = new double*[nRows];

  for (int i = 0; i < nRows; i++) {
    rho[i] = new double[nRows];
  }

  // create the sepSet basic structure
  if (direct) {
    sepSet = new int**[nRows];

    for (int i = 0; i < nRows; i++) {
      sepSet[i] = new int*[nRows];
    }

    // create the lenSepSet matrix
    lenSepSet = new int*[nRows];

    for (int i = 0; i < nRows; i++) {
      lenSepSet[i] = new int[nRows];
    }

    initializeMatrixNull(sepSet, nRows);
    initializeMatrixZero(lenSepSet, nRows);
  }

  //initialize matrix, numNeighbours and l
  initializeMatrix(matrix, nRows);
  initializeCutMap();
  initializeNeighbours(numNeighbours, nRows);
  initializeZero(means, nRows);
  initializeZero(standardDeviations, nRows);
  
}

/**
 *
 */
Graph::~Graph(void) {
  //empty the memory for matrix
  for (int i = 0; i < nRows; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;

  //empty the memory for cutMap
  for (int i = 0; i < nRows; i++) {
    delete[] cutMap[i];
  }
  delete[] cutMap;

  //empty the memory for bioData
  for (int i = 0; i < nRows; i++) {
    delete[] bioData[i];
  }
  delete[] bioData;

  //empty the memory for means
  delete[] means;

  //empty the memory for probeIDs
  delete[] probeIDs;

  //empty the memory for standardDeviations
  delete[] standardDeviations;

  //empty the memory for numNeighbours
  delete[] numNeighbours;

  //empty the memory for rho
  for (int i = 0; i < nRows; i++) {
    delete[] rho[i];
  }
  delete[] rho;
}


/**
 *  Compute the standard deviations for each node in the graph
 */
void Graph::computeStandardDeviations(void) {
  for (int r = 0; r < nRows; r++) {
    for (int c = 0; c < nCols; c++) {
      standardDeviations[r] += pow((bioData[r][c] - means[r]), 2);
    }

    standardDeviations[r] /= (double) nCols;
    standardDeviations[r] = sqrt(standardDeviations[r]);
  }
}

/**
 *  Compute the correlation coefficient of the base case, and store it in rho
 */
void Graph::computeCorrelations(void) {
  double covariance = 0.0;

  for (int i = 0; i < nRows; i++) {
    for (int j = 0; j < nRows; j++) {
      covariance = 0.0;

      for (int k = 0; k < nCols; k++) {
        covariance += (bioData[i][k] - means[i]) * (bioData[j][k] - means[j]);
      }

      //divide covariance by nCols
      covariance /= nCols;

      //covariance(i, j) / sigma(i) * sigma(j)
      rho[i][j] = covariance / (standardDeviations[i] * standardDeviations[j]);
    }
  }
}

/**
 *  Initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
 */
void Graph::initializeMatrix(bool** matrix, const int dim) {
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      if (i != j) {
        matrix[i][j] = true;
      } else {
        matrix[i][j] = false;
      }
    }
  }
}

/**
 *  Initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
 */
void Graph::initializeCutMap() {
  for (int i = 0; i < Graph::nRows; i++) {
    for (int j = 0; j < Graph::nRows; j++) {
      Graph::cutMap[i][j] = false;
    }
  }
}

/**
 *  Initialize the array numNeighbours with the value dim - 1, since the initial graph is connected
 */
void Graph::initializeNeighbours(int* numNeighbours, const int dim) {
  for (int i = 0; i < dim; i++) {
    numNeighbours[i] = dim - 1;
  }
}

/**
 *  Initialize the given array till dim to 0.0
 */
void Graph::initializeZero(double* array, const int dim) {
  for (int i = 0; i < dim; i++) {
    array[i] = 0.0;
  }
}

/**
 *  Initialize the int matrix to 0
 */
void Graph::initializeMatrixZero(int** matrix, const int dim) {
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      matrix[i][j] = 0;
    }
  }
}

/**
 *  Initialize the int 3D matrix to NULL
 */
void Graph::initializeMatrixNull(int*** matrix, const int dim) {
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      matrix[i][j] = NULL;
    }
  }
}


void testAndRemove(const int* neighbours, const int* subset, double correlationCoefficient, Graph* &g, const int r, const int c, const int l,
  const double alpha, const bool star, const bool directed) {
  double pVal;
  double const cutAt = 0.9999999;
  bool NAdelete = true;

  if (isnan(correlationCoefficient)) {
    correlationCoefficient = 0;
  }

  correlationCoefficient = min(cutAt, max(-cutAt, correlationCoefficient));

  pVal = sqrt(g->nCols - l - 3.0) * 0.5 * log((1 + correlationCoefficient) / (1 - correlationCoefficient));

  if (isnan(pVal)) {
    pVal = 0.0;
  }

  pVal = 2 * (1 - comulativeNormalDistribution(abs(pVal)));

  if (isnan(pVal)) {
    pVal = NAdelete ? 1.0 : 0.0;
  }

  //test d-separation
  if (pVal >= alpha) {
    
    // output_format = 0 means undirect, hence no separation set
    if (directed) {
      // save the size of the sepset at position r, c
      g->lenSepSet[r][c] = l;

      if (l > 0) {
        //create the vector for the sepset at position r, c
        g->sepSet[r][c] = new int[l];
        
        //save the sepset at position r, c
        for (int pos = 0; pos < l; pos++){
          //do not save on position c, r for saving memory
          g->sepSet[r][c][pos] = neighbours[subset[pos]];
        }
      }
    }

    if (star){
      // mark edges
      g->cutMap[r][c] = g->cutMap[c][r] = true;
    } else {
      //remove edges
      g->matrix[r][c] = g->matrix[c][r] = false;

      //decrement neighbours
      g->numNeighbours[r]--;
      g->numNeighbours[c]--;
    }
    
  }
}

/** Cuts all the edge we marked as to cut.
  
  @param Graph* &g
  The reference of the Graph object representing the gene network.
 */
void remove(Graph* &g) {
    for (int r = 0; r < g->nRows - 1; r++) {
        for (int c = r + 1; c < g->nRows; c++) {
            if (g->cutMap[r][c]) {
                g->numNeighbours[r]--;
                g->numNeighbours[c]--;
                g->matrix[r][c] = g->matrix[c][r] = false;
            }
        }
    }
}

/**
 *
 */
double getCorrelationCoefficient(const int* neighbours, const int* subset, const int l, Graph* &g, const int r, const int c, double** p) {
  double correlationCoefficient;

  if (l == 2) {
    double rhoRC, rhoRH1, rhoCH1;
    int h1 = neighbours[subset[0]];
    int h2 = neighbours[subset[1]];

    rhoRC = (g->rho[r][c] - (g->rho[r][h2] * g->rho[c][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[c][h2], 2)));
    rhoRH1 = (g->rho[r][h1] - (g->rho[r][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));
    rhoCH1 = (g->rho[c][h1] - (g->rho[c][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[c][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));

    correlationCoefficient = (rhoRC - (rhoRH1 * rhoCH1)) / sqrt((1 - pow(rhoRH1, 2)) * (1 - pow(rhoCH1, 2)));
  } else if (l == 3) {
    double rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3, rhoRC_H3, rhoRH1_H3, rhoRH2_H3, rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3;
    int h1 = neighbours[subset[0]];
    int h2 = neighbours[subset[1]];
    int h3 = neighbours[subset[2]];

    rhoRC_H3 = (g->rho[r][c] - (g->rho[r][h3] * g->rho[c][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[c][h3], 2)));
    rhoRH1_H3 = (g->rho[r][h1] - (g->rho[r][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
    rhoRH2_H3 = (g->rho[r][h2] - (g->rho[r][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
    rhoCH1_H3 = (g->rho[c][h1] - (g->rho[c][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
    rhoCH2_H3 = (g->rho[c][h2] - (g->rho[c][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
    rhoH1H2_H3 = (g->rho[h1][h2] - (g->rho[h1][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[h1][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));

    rhoRC_H2H3 = (rhoRC_H3 - (rhoRH2_H3 * rhoCH2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoCH2_H3, 2)));
    rhoRH1_H2H3 = (rhoRH1_H3 - (rhoRH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));
    rhoCH1_H2H3 = (rhoCH1_H3 - (rhoCH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoCH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));

    correlationCoefficient = (rhoRC_H2H3 - (rhoRH1_H2H3 * rhoCH1_H2H3)) / sqrt((1 - pow(rhoRH1_H2H3, 2)) * (1 - pow(rhoCH1_H2H3, 2)));
  } else {
    correlationCoefficient = correlations(r, c, g, neighbours, subset, l, p);
  }

  return correlationCoefficient;
}

/**
 *
 * Thanks to Ilya Bursov, http://stackoverflow.com/questions/19327847/n-choose-k-for-large-n-and-k
 */
void iterativeComb(int* neighbours, const int neighboursDim, const int l, Graph* &g, const int r, const int c, const double alpha, double** p,
  int* currentCombination, const bool star, const bool directed) {
  double coeff;

  for (int i = 0; i < neighboursDim; i++) {
    currentCombination[i] = i;
  }

  currentCombination[l - 1] = l - 1 - 1;

  do {
    if (currentCombination[l - 1] == (neighboursDim - 1)) {
      int i = l - 1 - 1;

      while (currentCombination[i] == (neighboursDim - l + i)) {
        i--;
      }

      currentCombination[i]++;

      for (int j = i + 1; j < l; j++) {
        currentCombination[j] = currentCombination[i] + j - i;
      }
    } else {
      currentCombination[l - 1]++;
    }

    coeff = getCorrelationCoefficient(neighbours, currentCombination, l, g, r, c, p);
    testAndRemove(neighbours, currentCombination, coeff, g, r, c, l, alpha, star, directed);
  } while (!g->cutMap[r][c] && g->matrix[r][c] && !((currentCombination[0] == (neighboursDim - l)) && (currentCombination[l - 1] == (neighboursDim - 1))));
}

/** Finds all the subsets adj(i)\{j} with cardinality equals to l (formally, |adj(i)\{j}| = l).
  
  @param Graph* &g
  The reference of the Graph object representing the gene network.

  @param const int i
  Index of the selected row, that also represents the "departure" node i of an edge i->j.

  @param const int j
  Index of the selected column, that also represents the "arrival" node j of an edge i->j.

  @param const int l
  Reached dimension actually taken into account for the subset cardinality.
  In this case l is greater than 1 (formally, l > 1).

  @param bool &hasWorked


  @param const double alpha
  Complement of the confidence interval

  @param int* neighbours
  Set of the indexes of the nearby nodes of the given edge i->j

  @param double** p
  Rhos.

  @return nothing.
*/
void findAllSubsets(Graph* &g, const int i, const int j, const int l, const double alpha, int* neighbours, double** p,
  int* currentCombination, const bool star, const bool directed) {
  int neighboursDim = 0;

  if (l == 0) {
    testAndRemove(NULL, NULL, g->rho[i][j], g, i, j, l, alpha, star, directed);
  } else if (l == 1) {
    //find neighbours (when l > 0) of i without considering j
    for (int k = 0; (g->matrix[i][j]) && !(g->cutMap[i][j]) && (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
      if (g->matrix[i][k] && (k != j)) {
        neighboursDim++;
        double correlationCoefficient = (g->rho[i][j] - (g->rho[i][k] * g->rho[j][k])) / sqrt((1 - pow(g->rho[i][k], 2)) * (1 - pow(g->rho[j][k], 2)));
        int pos[1];
        currentCombination[0] = k;
        pos[0] = 0;
        testAndRemove(currentCombination, pos, correlationCoefficient, g, i, j, l, alpha, star, directed);
      }
    }
  } else {
    //find neighbours (when l > 1) of i without considering j
    for (int k = 0; (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
      if ((g->matrix[i][k]) && (k != j)) {
        neighbours[neighboursDim++] = k;
      }
    }

    //look for all subset of length l
    iterativeComb(neighbours, neighboursDim, l, g, i, j, alpha, p, currentCombination, star, directed);
  }
}

/**
 *
 */
void skeleton(Graph* &g, const double alpha, const bool star, const bool directed) {
  int l = -1;
  bool hasWorked = true; //boolean to see that there is at least an arc i,j s.t. |ad j(C, i)\{j}| >= l TO CHECK FROM THE TEXT
  int* neighbours = new int[g->nRows]; //alloc the neighbours array (save time)
  int* currentCombination = new int[g->nRows]; //alloc an array for saving time
  double** p = new double*[g->nRows]; //alloc the rho for correlations() (save time)

  for (int i = 0; i < g->nRows; i++) {
    p[i] = new double[g->nRows];
  }

  //PC-algorithm
  while ((hasWorked) && (l < g->nRows)) {
    l++;
    hasWorked = false;

    for (int i = 0; i < g->nRows; i++) {
      for (int j = 0; j < g->nRows; j++) {
        //check if exists the arc between i and j
        if (g->matrix[i][j] && (g->numNeighbours[i] > l)) {
          hasWorked = true;
          findAllSubsets(g, i, j, l, alpha, neighbours, p, currentCombination, star, directed);
        }
      }
    }

    if (star) {
            remove(g);
            g->initializeCutMap();
        }
  }

  // free the memory    
  delete[] neighbours;

  for (int i = 0; i < g->nRows; i++) {
    delete[] p[i];
  }
  delete[] p;

  delete[] currentCombination;
}

bool comparator (const intpair &l, const intpair &r) {
  return l.first < r.first; 
}

/** Computes the continous density function.
  M_SQRT1_2 takes value 1/sqrt(2).
  
  @param const double value
  Value for which it will be computed its cumulative normal distribution.

  @return The cumulative normal distribution decimal value for the passed parameter.
*/
double comulativeNormalDistribution(const double value) {
  return 0.5 * erfc(-value * M_SQRT1_2);
}

/** Finds the correlation coefficient when l is greater than 1 (formally, when l > 1).
  
  @param const int a
  Index of the selected row, that also represents the "departure" node i of an edge i->j.

  @param const int b
  Index of the selected column, that also represents the "arrival" node j of an edge i->j.

  @param const Graph* g
  The Graph object representing the gene network.

  @param const int* neighbours
  Set of the nearby nodes.

  @param const int* subset
  Subset of the lookup indexes of the nearby nodes.
  Note that this subset has cardinality l.

  @param const int l
  Reached dimension actually taken into account for the subset cardinality.
  In this case l is greater than 1 (formally, l > 1).

  @param double ** p
  Rho value.
  
  @return The decimal value of the computed correlation for the edge a->b depending on the given neighbours' subset.
*/
double correlations(const int a, const int b, const Graph* g, const int* neighbours, const int* subset, const int l, double ** p) {
  int dim = l + 2;

  //initialization of p (looks like rho)
  for (int i = 0; i < dim - 1; i++) {
    for (int j = i + 1; j < dim; j++) {
      int first, second;

      if (i == 0) {
        first = a;
      } else if (i == 1) {
        first = b;
      } else {
        first = neighbours[subset[i - 2]];
      }

      if (j == 1) {
        second = b;
      } else {
        second = neighbours[subset[j - 2]];
      }

      p[i][j] = p[j][i] = g->rho[first][second];
    }
  }

  for (int k = 1; k <= l; k++) {
    for (int i = 0; i <= (l - k); i++) {
      for (int j = i + 1; j < (dim - k); j++) {
        p[i][j] = p[j][i] = (p[i][j] - p[i][dim - k] * p[j][dim - k]) / (sqrt((1 - pow(p[i][dim - k], 2)) * (1 - pow(p[j][dim - k], 2))));
      }
    }
  }

  return p[0][1];
}

/** Checks if a given string (of the form array of chars) whether representing a float number or not.
  
  @param const char* number
  String (or, rather, array of characters) that should represent a decimal number.

  @return TRUE if the string follows the correct format for representing a float number. FALSE otherwise.
*/
bool isFloat(const char* number) {
  bool noMorePoint = false;

  for (int i = 0; number[i] != '\0'; i++) {
    if ((number[i] < '0') || (number[i] > '9')) {
      if (number[i] == '.') {
        if (!noMorePoint) {
          noMorePoint = true;
        } else {
          return false;
        }
      } else {
        return false;
      }
    }
  }

  return true;
}

/** Counts the number of (uncutted) edges in the graph.

  @param bool** matrix
  Matrix of booleans representing a tabular form of the presence and absence of all the edges in the graph.
  The boolean located in the cell of i-th row and j-th column represents the presence/absence of the edge i->j.

  @param const int rows
  The number of rows in the matrix.

  @param const int cols
  The number of rows in the matrix.

  @return The number of uncutted edges.
*/
int countArcs(bool** matrix, const int rows, const int cols, const int mode) {
  int counter = -1;

  if (mode == 0) {
    counter = 0;

    for (int i = 0; i < rows; i++) {
      for (int j = i + 1; j < cols; j++) {
        if (matrix[i][j]) {
          counter++;
        }
      }
    } 
  } else {
    counter = 0;
    
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        if (matrix[i][j]) {
          counter++;
        }
      }
    }
  }

  return counter;
}
