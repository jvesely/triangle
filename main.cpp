#include <algorithm>
#include <iostream>
#include <list>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using num_type = double;
using vector = typename ::boost::numeric::ublas::vector<num_type>;
using scalar_vector = typename ::boost::numeric::ublas::scalar_vector<num_type>;
using matrix = typename ::boost::numeric::ublas::matrix<num_type>;

template<typename T>
bool is_zero(const T &);

template<>
bool is_zero<vector>(const vector &v)
{
	return !::std::any_of(v.begin(), v.end(),
		              [](const num_type n){ return n != 0; });
}

static inline matrix diag(const vector &v)
{
	matrix m(v.size(), v.size());
	m.clear();
	for (unsigned i = 0; i < v.size(); ++i)
		m(i,i) = v(i);
	return m;
}

static int find_pivot(const matrix &AD, const vector &p_prime)
{
	for (unsigned i = 0; i < AD.size2(); ++i) {
		const vector &candidate = column(AD, i);
		::std::cout << "pivot candidate: " << candidate << ::std::endl;
		if (norm_2(candidate) <= norm_2(candidate - p_prime))
			return i;
	}
	return -1;
}

int main(void)
{
	unsigned dims = 0;
	::std::cin >> dims;

	::std::vector<vector> points;
	do {
		vector v(dims);
		for (unsigned i = 0; i < dims; ++i) {
			num_type tmp;
			::std::cin >> tmp;
			v(i) = tmp;
		}
		if (!is_zero(v))
			points.push_back(v);
		else
			break;
	} while (1);

	matrix A(dims, points.size());

	for (unsigned i = 0; i < points.size(); ++i) {
		::std::cout << points[i] << ::std::endl;
		column(A, i) = points[i];
	}
	::std::cout << "A: " << A << ::std::endl;

	static const num_type n = points.size();
	static const scalar_vector e_over_n(n, 1.0 / n);
	vector d = e_over_n;

	{
		/* Step 1 compute D, and p'.
		 * Test if p' is close enough to origin.
		 */
		matrix D = diag(d) * n;
		matrix AD = prod(A,  D);
		vector p_prime = prod(AD, e_over_n);

		::std::cout << "d: " << d << ::std::endl;
		::std::cout << "D: " << D << ::std::endl;
		::std::cout << "AD: " << AD << ::std::endl;
		::std::cout << "p': " << p_prime << ::std::endl;
		if (is_zero(p_prime)) {
			::std::cout << "Origin IS in the convex hull!\n";
			return 0;
		}

		/* Step 1.5: Find pivot (vertex closer to origin than p') */
		const int pivot_index = find_pivot(AD, p_prime);
		if (pivot_index == -1) {
			::std::cout << "Origin IS NOT in the convex hull!\n";
			return 1;
		}
		vector pivot = column(AD, pivot_index);
		::std::cout << "pivot: " << pivot << ::std::endl;
	}
	return 0;
}
