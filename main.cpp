#include <algorithm>
#include <iostream>
#include <list>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using num_type = double;
using vector = typename ::boost::numeric::ublas::vector<num_type>;
using matrix = typename ::boost::numeric::ublas::matrix<num_type>;

int main(void)
{
	unsigned dims = 0;
//	::std::cout << "Set number of dimensions: ";
	::std::cin >> dims;

	::std::vector<vector> points;
	do {
		vector v(dims);
		for (unsigned i = 0; i < dims; ++i) {
			num_type tmp;
			::std::cin >> tmp;
			v(i) = tmp;
		}
		if (::std::any_of(v.begin(), v.end(),
		                  [](const num_type n){ return n != 0; }))
			points.push_back(v);
		else
			break;
	} while (1);

	matrix A(dims, points.size());

	for (unsigned i = 0; i < points.size(); ++i) {
		::std::cout << points[i] << ::std::endl;
		for (unsigned j = 0; j < dims; ++j) {
			A(j,i) = points[i][j];
		}
	}
	::std::cout << A << ::std::endl;
	return 0;
}
