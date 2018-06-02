template<typename T> void IsNonZero(const std::string& location,
	       			    const std::string& what,
 	       			    const T& value) {

  if ( !value )
      throw(Exceptions::Exception(Exceptions::nonRecoverable, what+" is zero", location));
}

template<typename T> void IsGreater(const std::string& location,
	       			    const std::string& what,
 	       			    const T& value,
				    const T& limit,
				    const bool strictly) {

  if ( !strictly && value < limit ) {
   throw(Exceptions::Exception(Exceptions::nonRecoverable,
	 what+" is strictly lower than "+std::to_string(limit),
         location));
  } else if ( strictly && value <= limit ) {
   throw(Exceptions::Exception(Exceptions::nonRecoverable,
	 what+" is lower or equal to "+std::to_string(limit),
         location));
  }
}

template<typename T> void IsLower(const std::string& location,
	       			  const std::string& what,
 	       			  const T& value,
				  const T& limit,
				  const bool strictly) {

  if ( !strictly && value > limit ) {
   throw(Exceptions::Exception(Exceptions::nonRecoverable,
	 what+" is strictly greater than "+std::to_string(limit),
         location));
  } else if ( strictly && value >= limit ) {
   throw(Exceptions::Exception(Exceptions::nonRecoverable,
	 what+" is greater or equal to "+std::to_string(limit),
         location));
  }
}

template<typename T> void IsEqual(const std::string& location,
	       			  const std::string& what,
 	       			  const T& a,
 	       			  const T& b) {

  if ( a != b )
      throw(Exceptions::Exception(Exceptions::nonRecoverable, what+" are not equal", location));
}

template<typename T> void IsNotEmpty(const std::string& location,
				     const std::string& what,
				     const std::vector<T>& vector) {

  if ( !vector.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable, what+" is empty", location));
}

template<typename T> void HasAtLeast(const std::string& location,
				     const std::string& what,
				     const std::vector<T>& vector,
				     const size_t n) {

  if ( vector.size() < n )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" contains less than "+std::to_string(n)+" elements",
	    location));
}

template<typename T> void HasAtMost(const std::string& location,
				    const std::string& what,
				    const std::vector<T>& vector,
				    const size_t n) {

  if ( vector.size() > n )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" contains more than "+std::to_string(n)+" elements",
	    location));
}

template<typename T> void SameSize(const std::string& location,
				   const std::string& what,
				   const std::vector<T>& va,
				   const std::vector<T>& vb) {

  if ( va.size() != vb.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" have different dimensions",
	    location));
}

template<typename T> void SameSize(const std::string& location,
				   const std::string& what,
				   const std::vector<std::vector<T>>& ma,
				   const std::vector<std::vector<T>>& mb) {

  if ( ma.size() != mb.size() || ma[0].size() != mb[0].size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" have different dimensions",
	    location));
}

template<typename T> void Multiplicable(const std::string& location,
				   	const std::string& what,
				   	const std::vector<std::vector<T>>& ma,
				   	const std::vector<std::vector<T>>& mb) {

  if ( ma[0].size() != mb.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" have incompatible dimensions to be multiplied",
	    location));
}

template<typename T> void IsSquare(const std::string& location,
				   const std::string& what,
				   const std::vector<std::vector<T>>& matrix) {

  if ( matrix.size() != matrix[0].size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" is not square",
	    location));
}

template<typename T> void IsSymmetric(const std::string& location,
				      const std::string& what,
				      const std::vector<std::vector<T>>& matrix) {

  try {
    Assert::IsSquare(location, what, matrix);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  what+" cannot be symmetric"+std::string(e.what()),
	  location));
  }

  size_t i, j;
  for (i = 0; i < matrix.size(); i++)
    for (j = i+1; j < matrix.size(); j++)
      if ( matrix[i][j] != matrix[j][i] )
          throw(Exceptions::Exception(Exceptions::nonRecoverable,
		what+" is not symmetric",
		location));
}

template<typename T1, typename T2> void Contains(const std::string& location,
				      		 const std::string& what,
				      		 const std::map<T1, T2>& map,
				    		 const T1& key) {

  if ( map.find(key) == map.end() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
				  what+" does not contain "+key,
				  location));
}
