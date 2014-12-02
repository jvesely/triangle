/*
 * Copyright (c) 2014 Jan Vesely
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * - The name of the author may not be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cfloat>
#include <cmath>
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


/** Check whether vector is all 0
 * @param v Vector to check
 * @retval True, if all elements in the vector are 0
 *         False, otherwise.
 */
bool is_zero(const vector &v)
{
	return !::std::any_of(v.begin(), v.end(),
		              [](const num_type n){ return n != 0; });
}

/** Function creates diagonal matrix from vector v
 * @param v Vector to be put on diagonal
 * @return Square matrix with elements of v on diagonal, other elements are 0.
 */
static inline matrix diag(const vector &v)
{
	matrix m(v.size(), v.size());
	/* Constructor creates uninitialized matrix, we need to clear it */
	m.clear();
	for (unsigned i = 0; i < v.size(); ++i)
		m(i,i) = v(i);
	return m;
}

/** Look for vertices that are closer to O than p_prime.
 * @param AD matrix, columns are vertices.
 * @param p_prime reference point.
 * @return index of column holding the vertex with greatest distance difference
 *         >= 0
 *         -1, if there is no vertex with distance to O greater or equal to
 *         distance to p_prime.
 */
static int find_pivot(const matrix &AD, const vector &p_prime)
{
	double norm_diff = 0;
	int index = -1;
	for (unsigned i = 0; i < AD.size2(); ++i) {
		const vector &candidate = column(AD, i);
		const double diff =
		    norm_2(candidate - p_prime) - norm_2(candidate);
		::std::cout << "pivot candidate: " << candidate << "(" << diff
		            << ")" << ::std::endl;
		if (diff >= norm_diff) {
			norm_diff = diff;
			index = i;
		}
	}
	return index;
}

/** Helper function to do vector * matrix * vector multiplication.
 * @param v the first operand
 * @param m the second operand
 * @param u the third operand
 * @return V x M x U
 */
static num_type multi_prod(const vector &v, const matrix &m, const vector &u)
{
	vector tmp = prod(v, m);
	matrix upgrade(1, tmp.size());
	row(upgrade, 0) = tmp;
	vector tmp2 = prod(upgrade, u);
	return tmp2(0);
}

/** Scan the range and find minimum
 * @param f function, we use template to work around lack of lambda-function
 *                    type
 * @param low lower bound
 * @param high upper bound
 * @return Point that is closes to real point where @p f attains its minimum
 *
 * This function assumes that @p f is monotonic between bounds and min.
 */
template<typename T>
num_type find_minimum(const T& f, num_type low, num_type high)
{
	while (1) {
		assert(high >= low);
		if (::std::nextafter(low, DBL_MAX) >= high)
			return f(low) < f(high) ? low : high;
		// Or some other nice value
		const num_type slice = (high - low) / 4;
#if 0
		::std::cout << "Finding minimum (slice = " << slice << ") in ("
		            << low << ", " << high << "): (" << f(low) << ", "
		            << f(high) << ")\n";
#endif
		num_type prev_res = f(low);
		num_type it = low + slice;
		while (prev_res > f(it) && it < high) {
			prev_res = f(it);
			it += slice;
			assert(it <= high);
		}
		/* Now 'it' is the first increase sample. It can be: */
		/* a) Right at the beginning,
		 * the minimum is in the first slice */
		if (it - slice <= low) {
			high = it;
			continue;
		}
		/* b) At the end, the minimum is in the last slice */
		if (it + slice > high) {
			low = it - slice;
			continue;
		}
		/* c) Somewhere in the middle, the minimum is within
		 * w slices of the previous value. */
		low = it - 2 * slice;
		high = it;
//		assert(f(it - slice) <= f(low));
//		assert(f(it - slice) <= f(high));

	}
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

	/* Keep limit as safeguard */
	while (it < 1200) {
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
		::std::cout << "d': " << f(a) << ::std::endl;
		if (is_zero(vector(d - f(a)))) {
			::std::cout << "Can't get new point, exiting\n";
			::std::cout << "We should be close enough to say "
			            << "that origin IS in the convex hull\n";
			break;
		}
		d = f(a);
#if 0
		::std::cout << "y(a): " << y(a) << ::std::endl;
		::std::cout << "Dy(a): " << prod(D,y(a)) << ::std::endl;
		::std::cout << "eTDy(a): " << multi_prod(scalar_vector(n, 1), D, y(a)) << ::std::endl;
		::std::cout << "f(a): " << prod(D,y(a))/multi_prod(scalar_vector(n, 1), D, y(a)) << ::std::endl;
#endif
	}
	return 0;
}
