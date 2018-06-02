#ifndef VERTEX_HH
#define VERTEX_HH

// C++ includes
#include <stdlib.h>
#include <vector>
#include <cmath>

/** @brief Defines an n-dimensional vertex and its mapping onto a sclar field.
 *
 *	   If the verte is on a Grid, it stores the Grid indices.
 */
class Vertex {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Vertex();

  /** @brief Normal constructor, sets vector of coordinates and its mapping 
   *  @param	coords		Coordinates of the vertex
   *  @param	value		Mapping of the scalar field in the vertex
   *  @param	indices		Set of indices that specifies the position of the vertex on the grid
   */
  Vertex(const std::vector<double>& coords,
	 const double& value=1.,
	 const std::vector<size_t>& indices=std::vector<size_t>());

  /** @brief Size constructor, sets vector size */
  Vertex(const size_t n);

  /** @brief Copy constructor */
  Vertex(const Vertex& vertex);

  /** @brief Equality operator */
  Vertex& operator=(const Vertex& vertex);

  /** @brief Returns the id^th component of the vector field at which the vertex is defined */
  double& operator[](const size_t id)				{ return _coords[id]; }

  /** @brief Returns the id^th component of the vector field at which the vertex is defined */
  const double& operator[](const size_t id) const		{ return _coords[id]; }

  /** @brief Destructor */
  ~Vertex();

  /** @brief Returns the dimension of the vector field in which the vertex is defined */
  const size_t size() const					{ return _coords.size(); }

  /** @brief Returns the coordinates of the vertex */
  const std::vector<double>& GetCoordinates() const		{ return _coords; }

  /** @brief Returns the coordinates of the vertex */
  void SetCoordinates(const std::vector<double>& coords)	{ _coords = coords; }

  /** @brief Returns the value of the field at the vertex */
  const double& GetValue() const				{ return _value; }

  /** @brief Returns the value of the field at the vertex */
  void SetValue(const double& value)				{ _value = value; }

  /** @brief Returns the dimension of space in which the vertex lives */
  const size_t GetDimension()	const				{ return _coords.size(); }

  /** @brief Returns the value of the field at the vertex */
  const std::vector<size_t>& GetIndices() const			{ return _indices; }

  /** @brief Returns the value of the field at the vertex */
  void SetIndices(const std::vector<size_t>& indices)		{ _indices = indices; }

  /** @brief Returns the L2 vector norm of the position vector */
  double GetNorm() const;

 private:
  std::vector<double>	_coords;	///< n-vector, position of the vertex
  double		_value;		///< Value of the scalar field in the vertex
  std::vector<size_t>	_indices;	///< n-vector, indices of the point within a grid
};

#endif
