#include <algorithm>
#include <iostream>
#include <list>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using num_type = double;
using vector = typename ::boost::numeric::ublas::vector<num_type>;

int main(void)
{
	unsigned dims = 0;
//	::std::cout << "Set number of dimensions: ";
	::std::cin >> dims;

	static const ::boost::numeric::ublas::zero_vector<num_type> zero_v(dims);
	::std::list<vector> points;
	vector v(dims);
	do {
		v.clear();
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

	for (const auto &v: points)
		::std::cout << v << ::std::endl;

	return 0;
}
