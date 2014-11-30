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
using unit_vector = typename ::boost::numeric::ublas::unit_vector<num_type>;
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

static num_type multi_prod(const vector &v, const matrix &m, const vector &u)
{
	vector tmp = prod(v, m);
	matrix upgrade(1, tmp.size());
	row(upgrade, 0) = tmp;
	vector tmp2 = prod(upgrade, u);
	return tmp2(0);
}

template<typename T>
num_type find_minimum(const T& f, num_type low, num_type high)
{
	//TODO actually find minimum
	return (low + high) / 2;
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

	::std::cout << "Polytope points:\n";
	for (unsigned i = 0; i < points.size(); ++i) {
		::std::cout << points[i] << ::std::endl;
		column(A, i) = points[i];
	}
	::std::cout << "A: " << A << ::std::endl;

	static const num_type n = points.size();
	static const scalar_vector e_over_n(n, 1.0 / n);
	vector d = e_over_n;
	unsigned it = 0;

	{
		::std::cout << "\nIteration: " << it++ << ::std::endl;
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
		::std::cout << "p' distance form O: " << norm_2(p_prime)
		            << ::std::endl;

		/* Step 1.5: Find pivot (vertex closer to origin than p') */
		const int pivot_index = find_pivot(AD, p_prime);
		if (pivot_index == -1) {
			::std::cout << "Origin IS NOT in the convex hull!\n";
			return 1;
		}
		vector pivot = column(AD, pivot_index);
		::std::cout << "pivot(" << pivot_index << "): " << pivot
		            << ::std::endl;

		/* Step 2: create distance function. */
		const vector u = unit_vector(n, pivot_index) - e_over_n;
		auto y = [&](num_type a){ return e_over_n + (u * a); };

		::std::cout << "y(0): " << y(0) << ::std::endl;
		::std::cout << "y(1): " << y(1) << ::std::endl;

		auto f = [&](num_type a) -> vector {
			return prod(D, y(a)) /
			       multi_prod(scalar_vector(n, 1), D, y(a));
		};
		auto f_norm_sq = [&](num_type a) {
			auto tmp = norm_2(prod(A, f(a)));
			return tmp * tmp;
		};

		/* Step 2.5: find closest point. */
		num_type a = find_minimum(f_norm_sq, 0.0, 1.0);
		::std::cout << "Found a*: " << a << ::std::endl;
		::std::cout << "Norm at a*: " << f_norm_sq(a) << ::std::endl;

		/* Step 3: set new d */
		d = f(a);
		::std::cout << "d': " << d << ::std::endl;
	}
	return 0;
}
