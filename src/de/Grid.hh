#ifndef GRID_HH
#define GRID_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// ROOT includes
#include "TPolyLine.h"

// Additional includes
#include "Matrix.hh"
#include "Vertex.hh"

/** @brief Data structure for a grid of points and their mapping onto a scalar field.
 */
class Grid {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Grid();

  /** @brief Vertex constructor, constructs the grid from a vertex vector
   *
   *  @param	vertices	Array of vertices and their mapping
   */
  Grid(const std::vector<Vertex>& vertices);

  /** @brief Matrix constructor, sets matrix of vertices and their mapping on a grid
   * 
   *  @param	mat	nxN-matrix of vertices, n = dimension, N = number of vertices
   *  @param	vec	N-vector of the mapping, one per vertex, same ordering
   */
  Grid(const Matrix<double>& mat,
       const std::vector<double>& vec);

  /** @brief Function constructor, constructs the grid from a function
   *
   *  @param	fcn	Function of vector field that returns a scalar field
   *  @param	number	Number of steps between min and max values
   *  @param	lower	Lower bounds of the grid in each dimension
   *  @param	upper	Upper bounds of the grid in each dimension
   */
  Grid(double (*fcn)(const std::vector<double>& val),
       const std::vector<size_t>& number,
       const std::vector<double>& lower,
       const std::vector<double>& upper);

  /** @brief Random constructor, constructs an N_1*...*N_n random grid between lower and upper
   *
   *  @param	number		Number of steps between min and max values
   *  @param	lower		Lower bounds of the grid in each dimension
   *  @param	upper		Upper bounds of the grid in each dimension
   *  @param	algo		Alogrithm for random number generation
   */
  Grid(const std::vector<size_t>& number,
       const std::vector<double>& lower,
       const std::vector<double>& upper,
       const std::string algo="flat");

  /** @brief Copy constructor */
  Grid(const Grid& grid);

  /** @brief Equality operator */
  Grid& operator=(const Grid& grid);

  /** @brief Not operator */
  bool operator!() const				{ return !bool(_dim); }

  /** @brief Returns the vertex  */
  Vertex& operator[](const size_t id)			{ return _vertices[id]; }

  /** @brief Returns the vertex  */
  const Vertex& operator[](const size_t id) const	{ return _vertices[id]; }

  /** @brief Subscript operator */
  Vertex& operator[](const std::vector<size_t> indices)	{ return _vertices[ID(indices)]; }

  /** @brief Subscript operator */
  const Vertex& operator[](const std::vector<size_t> indices) const
							{ return _vertices[ID(indices)]; }

  /** @brief Destructor */
  ~Grid();

  /** @brief Returns the dimension of the vector field in which the vertex is defined */
  const size_t size() const				{ return _vertices.size(); }

  /** @brief Initializes the grid */
  void Initialize();

  /** @brief Initializes the grid, provided a matrix of vertices and their mapping
   *
   *  @param	mat	nxN-matrix of vertices, n = dimension, N = number of vertices
   *  @param	vec	N-vector of the mapping, one per vertex, same ordering
   */
  void InitializeFromMatrix(const Matrix<double>& mat, const std::vector<double>& vec);

  /** @brief Initializes the grid, provided a scalar function defined within boundaries
   *
   *  @param	fcn	Function of vector field that returns a scalar field
   */
  void InitializeFromFunction(double (*fcn)(const std::vector<double>& val));

  /** @brief Initializes the grid to random values of scalar field
   *
   *  @param	algo		Alogrithm for random number generation
   */
  void InitializeRandom(const std::string& algo);

  /** @brief Find the closest cell of the grid to a given point of coordinates v
   *
   *  @param	v	Vector of coordinates of the n-point
   *
   *  @return		Cell of surrounding vertices
   */
  Grid GetBestCell(const std::vector<double>& v) const;

  /** @brief Find the IDs of the points on the grid that directly preceed v */
  std::vector<size_t> GetLowerIndices(const std::vector<double>& v) const;

  /** @brief Returns a cell of the grid given the lower bound indicies */
  Grid GetCell(const std::vector<size_t>& indices) const;

  /** @brief Returns all the cells of the grid in an array */
  std::vector<Grid> GetCellArray() const;

  /** @brief Returns the volume of one of the cell of the grid */
  double GetCellVolume() const;

  /** @brief Returns the dimension of the space in which the grid is constructed */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the array of number of points in each dimension */
  const std::vector<size_t>& GetNpointsArray() const	{ return _number; }

  /** @brief Returns the number of points in dimension i */
  const size_t& GetNpoints(const size_t i) const	{ return _number[i]; }

  /** @brief Returns the array of step sizes in each dimension */
  const std::vector<size_t>& GetStepArray() const	{ return _step; }

  /** @brief Returns the step size in dimension i */
  const size_t& GetStep(const size_t i) const		{ return _step[i]; }

  /** @brief Returns the array of period size in each dimension */
  const std::vector<double>& GetPeriodArray() const	{ return _period; }

  /** @brief Returns the period size in dimension i */
  const double& GetPeriod(const size_t i) const		{ return _period[i]; }

  /** @brief Returns the array of lower bounds in each dimension */
  const std::vector<double>& GetLowerBoundArray() const	{ return _lower; }

  /** @brief Returns the lower bound in dimension i */
  const double& GetLowerBound(const size_t i) const	{ return _lower[i]; }

  /** @brief Returns the array of upper bounds in each dimension */
  const std::vector<double>& GetUpperBoundArray() const	{ return _upper; }

  /** @brief Returns the upper bound in dimension i */
  const double& GetUpperBound(const size_t i) const	{ return _upper[i]; }

  /** @brief Sets the coordinates of the vertex */
  void SetVertex(const std::vector<size_t>& indices, const Vertex& vertex)
							{ _vertices[ID(indices)] = vertex; }

  /** @brief Returns the array of all the vertices */
  const std::vector<Vertex>& GetVertexArray() const	{ return _vertices; }

  /** @brief Returns the array of all the vertices, allows modifications */
  std::vector<Vertex>& GetVertexArray()			{ return _vertices; }

  /** @brief Returns the position of a point of index id */
  const std::vector<double>& GetPositionArray(const size_t id) const
							{ return _vertices[id].GetCoordinates(); }

  /** @brief Returns the position in dimension i for a point index id in this dimension */
  double GetPosition(const size_t i, const size_t id) const;

  /** @brief Returns the position in dimensions vi for a point index ids in these dimension */
  std::vector<double> GetPosition(const std::vector<size_t>& vi,
				  const std::vector<size_t>& ids) const;

  /** @brief Returns the volume occupied by the grid */
  double GetVolume() const;

  /** @brief Returns the id of the element within the vertex matrix */
  size_t ID(const std::vector<size_t>& indices) const;

  /** @brief Returns the indices of the element within the vertex matrix */
  std::vector<size_t> Indices(size_t id) const;

  /** @brief Returns the indices of the element within the vertex matrix */
  bool IsInside(const std::vector<double>& v) const;

  /** @brief Returns an array of polygons to be drawn in ROOT (only 2D) */
  std::vector<TPolyLine*> Polygons() const;

 private:
  size_t		_dim;		///< Dimension of the space
  std::vector<Vertex>	_vertices;	///< N_1*...*N_n-vector of vertices (Vertex object)
  std::vector<size_t>	_number;	///< n-vector, number of steps on each axis
  std::vector<size_t>	_step;		///< n-vector, stepping of each axis
  std::vector<double>	_period;	///< n-vector, length of a period of each axis
  std::vector<double>	_lower;		///< n-vector, lower bounds of the grid
  std::vector<double>	_upper;		///< n-vector, upper bounds of the grid
};

#endif
