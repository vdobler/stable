// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Package sort provides primitives for sorting slices and user-defined
// collections.
package sort

// A type, typically a collection, that satisfies sort.Interface can be
// sorted by the routines in this package.  The methods require that the
// elements of the collection be enumerated by an integer index.
type Interface interface {
	// Len is the number of elements in the collection.
	Len() int
	// Less returns whether the element with index i should sort
	// before the element with index j.
	Less(i, j int) bool
	// Swap swaps the elements with indexes i and j.
	Swap(i, j int)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Insertion sort
func insertionSort(data Interface, a, b int) {
	for i := a + 1; i < b; i++ {
		for j := i; j > a && data.Less(j, j-1); j-- {
			data.Swap(j, j-1)
		}
	}
}

// siftDown implements the heap property on data[lo, hi).
// first is an offset into the array where the root of the heap lies.
func siftDown(data Interface, lo, hi, first int) {
	root := lo
	for {
		child := 2*root + 1
		if child >= hi {
			break
		}
		if child+1 < hi && data.Less(first+child, first+child+1) {
			child++
		}
		if !data.Less(first+root, first+child) {
			return
		}
		data.Swap(first+root, first+child)
		root = child
	}
}

func heapSort(data Interface, a, b int) {
	first := a
	lo := 0
	hi := b - a

	// Build heap with greatest element at top.
	for i := (hi - 1) / 2; i >= 0; i-- {
		siftDown(data, i, hi, first)
	}

	// Pop elements, largest first, into end of data.
	for i := hi - 1; i >= 0; i-- {
		data.Swap(first, first+i)
		siftDown(data, lo, i, first)
	}
}

// Quicksort, following Bentley and McIlroy,
// ``Engineering a Sort Function,'' SP&E November 1993.

// medianOfThree moves the median of the three values data[a], data[b], data[c] into data[a].
func medianOfThree(data Interface, a, b, c int) {
	m0 := b
	m1 := a
	m2 := c
	// bubble sort on 3 elements
	if data.Less(m1, m0) {
		data.Swap(m1, m0)
	}
	if data.Less(m2, m1) {
		data.Swap(m2, m1)
	}
	if data.Less(m1, m0) {
		data.Swap(m1, m0)
	}
	// now data[m0] <= data[m1] <= data[m2]
}

func swapRange(data Interface, a, b, n int) {
	for i := 0; i < n; i++ {
		data.Swap(a+i, b+i)
	}
}

func doPivot(data Interface, lo, hi int) (midlo, midhi int) {
	m := lo + (hi-lo)/2 // Written like this to avoid integer overflow.
	if hi-lo > 40 {
		// Tukey's ``Ninther,'' median of three medians of three.
		s := (hi - lo) / 8
		medianOfThree(data, lo, lo+s, lo+2*s)
		medianOfThree(data, m, m-s, m+s)
		medianOfThree(data, hi-1, hi-1-s, hi-1-2*s)
	}
	medianOfThree(data, lo, m, hi-1)

	// Invariants are:
	//	data[lo] = pivot (set up by ChoosePivot)
	//	data[lo <= i < a] = pivot
	//	data[a <= i < b] < pivot
	//	data[b <= i < c] is unexamined
	//	data[c <= i < d] > pivot
	//	data[d <= i < hi] = pivot
	//
	// Once b meets c, can swap the "= pivot" sections
	// into the middle of the slice.
	pivot := lo
	a, b, c, d := lo+1, lo+1, hi, hi
	for {
		for b < c {
			if data.Less(b, pivot) { // data[b] < pivot
				b++
			} else if !data.Less(pivot, b) { // data[b] = pivot
				data.Swap(a, b)
				a++
				b++
			} else {
				break
			}
		}
		for b < c {
			if data.Less(pivot, c-1) { // data[c-1] > pivot
				c--
			} else if !data.Less(c-1, pivot) { // data[c-1] = pivot
				data.Swap(c-1, d-1)
				c--
				d--
			} else {
				break
			}
		}
		if b >= c {
			break
		}
		// data[b] > pivot; data[c-1] < pivot
		data.Swap(b, c-1)
		b++
		c--
	}

	n := min(b-a, a-lo)
	swapRange(data, lo, b-n, n)

	n = min(hi-d, d-c)
	swapRange(data, c, hi-n, n)

	return lo + b - a, hi - (d - c)
}

func quickSort(data Interface, a, b, maxDepth int) {
	for b-a > 7 {
		if maxDepth == 0 {
			heapSort(data, a, b)
			return
		}
		maxDepth--
		mlo, mhi := doPivot(data, a, b)
		// Avoiding recursion on the larger subproblem guarantees
		// a stack depth of at most lg(b-a).
		if mlo-a < b-mhi {
			quickSort(data, a, mlo, maxDepth)
			a = mhi // i.e., quickSort(data, mhi, b)
		} else {
			quickSort(data, mhi, b, maxDepth)
			b = mlo // i.e., quickSort(data, a, mlo)
		}
	}
	if b-a > 1 {
		insertionSort(data, a, b)
	}
}

// Sort sorts data.
// It makes one call to data.Len to determine n, and O(n*log(n)) calls to
// data.Less and data.Swap. The sort is not guaranteed to be stable.
func Sort(data Interface) {
	// Switch to heapsort if depth of 2*ceil(lg(n+1)) is reached.
	n := data.Len()
	maxDepth := 0
	for i := n; i > 0; i >>= 1 {
		maxDepth++
	}
	maxDepth *= 2
	quickSort(data, 0, n, maxDepth)
}

type reverse struct {
	// This embedded Interface permits Reverse to use the methods of
	// another Interface implementation.
	Interface
}

// Less returns the opposite of the embedded implementation's Less method.
func (r reverse) Less(i, j int) bool {
	return r.Interface.Less(j, i)
}

// Reverse returns the reverse order for data.
func Reverse(data Interface) Interface {
	return &reverse{data}
}

// IsSorted reports whether data is sorted.
func IsSorted(data Interface) bool {
	n := data.Len()
	for i := n - 1; i > 0; i-- {
		if data.Less(i, i-1) {
			return false
		}
	}
	return true
}

// Convenience types for common cases

// IntSlice attaches the methods of Interface to []int, sorting in increasing order.
type IntSlice []int

func (p IntSlice) Len() int           { return len(p) }
func (p IntSlice) Less(i, j int) bool { return p[i] < p[j] }
func (p IntSlice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// Sort is a convenience method.
func (p IntSlice) Sort() { Sort(p) }

// Float64Slice attaches the methods of Interface to []float64, sorting in increasing order.
type Float64Slice []float64

func (p Float64Slice) Len() int           { return len(p) }
func (p Float64Slice) Less(i, j int) bool { return p[i] < p[j] || isNaN(p[i]) && !isNaN(p[j]) }
func (p Float64Slice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// isNaN is a copy of math.IsNaN to avoid a dependency on the math package.
func isNaN(f float64) bool {
	return f != f
}

// Sort is a convenience method.
func (p Float64Slice) Sort() { Sort(p) }

// StringSlice attaches the methods of Interface to []string, sorting in increasing order.
type StringSlice []string

func (p StringSlice) Len() int           { return len(p) }
func (p StringSlice) Less(i, j int) bool { return p[i] < p[j] }
func (p StringSlice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

// Sort is a convenience method.
func (p StringSlice) Sort() { Sort(p) }

// Convenience wrappers for common cases

// Ints sorts a slice of ints in increasing order.
func Ints(a []int) { Sort(IntSlice(a)) }

// Float64s sorts a slice of float64s in increasing order.
func Float64s(a []float64) { Sort(Float64Slice(a)) }

// Strings sorts a slice of strings in increasing order.
func Strings(a []string) { Sort(StringSlice(a)) }

// IntsAreSorted tests whether a slice of ints is sorted in increasing order.
func IntsAreSorted(a []int) bool { return IsSorted(IntSlice(a)) }

// Float64sAreSorted tests whether a slice of float64s is sorted in increasing order.
func Float64sAreSorted(a []float64) bool { return IsSorted(Float64Slice(a)) }

// StringsAreSorted tests whether a slice of strings is sorted in increasing order.
func StringsAreSorted(a []string) bool { return IsSorted(StringSlice(a)) }

// -------------------------------------------------------------------------
// Stable sorting

// MergeSort sorts data by merging consecutive blocks by calling merge.
// Cut determines when to switch from recursive merge sort to insertion sort.
// The maximum recursion depth is returned.
func MergeSort(data Interface, merge func(data Interface, a, m, b int, d int) int, cut int, a, b int, depth int) (maxRd int) {
	if b-a <= 1 {
		return depth
	}

	if b-a < cut {
		insertionSort(data, a, b)
		return depth
	} else {
		m := a + (b-a)/2
		ld := MergeSort(data, merge, cut, a, m, depth+1)
		rd := MergeSort(data, merge, cut, m, b, depth+1)
		md := merge(data, a, m, b, depth+1)
		if ld > rd {
			maxRd = ld
		} else {
			maxRd = rd
		}
		if md > maxRd {
			maxRd = md
		}
		return maxRd
	}
}

// -------------------------------------------------------------------------
// Block rotation for stable sorting

// Rotate exchanges the two consecutive blocks data[a,m) and data[m,b).
func rotate(data Interface, a, m, b int) {
	// RotateSpeed(data, a, m, b)
	RotateSwap(data, a, m, b)
	// rotateOptimal(data, a, m, b)
	// rotateSimple(data, a, m, b)
	// RotateJuggling(data, a, m, b)
}

// RotateSimple exchanges the two consecutive blocks data[a,m) and data[m,b).
// It performs $\approx 2(b-a)$ swap operations.
func rotateSimple(data Interface, a, m, b int) {
	reverseBlock(data, a, m)
	reverseBlock(data, m, b)
	reverseBlock(data, a, b)
}
func RotateSimple(data Interface, a, m, b int) { rotateSimple(data, a, m, b) }

// ReverseBlock reverses the elements of the subsequence data[a,b).
//    ..., data[a-1], data[a], data[a+1], ..., data[b-2], data[b-1], data[b], ...
// is turned into
//    ..., data[a-1], data[b-1], data[b-2], ..., data[a+1], data[a], data[b], ...
func reverseBlock(data Interface, a, b int) {
	b--
	for a < b {
		data.Swap(a, b)
		a += 1
		b -= 1
	}
}

// RotateOptimal exchanges the two consecutive blocks data[a,m) and data[m,b).
// It performs $O((b-a) \log(\min(m-a,b-b)))$ swap operations ?
// The implementation is a direct translation of gcc 4.6.3 libstdc++ __rotate()
// function for bidirectional iterators in stl_algo.h
// (See http://gcc.gnu.org/onlinedocs/gcc-4.6.3/libstdc++/api/a01045_source.html)
func RotateOptimal(data Interface, a, m, b int) { rotateOptimal(data, a, m, b) }
func rotateOptimal(data Interface, a, m, b int) {
	if a == m || m == b {
		return
	}
	n := b - a
	k := m - a
	p := a

	if k == n-k {
		// Both blocks are of length k (=n/2) and rotating is easy.
		for i := 0; i < k; i++ {
			data.Swap(a+i, m+i)
		}
		return
	}

	for {
		if k < n-k {
			if k == 1 {
				for i := p; i < p+n-1; i++ {
					data.Swap(i, i+1)
				}
				return
			}
			q := p + k
			for i := 0; i < n-k; i++ {
				data.Swap(p, q)
				p++
				q++
			}
			n %= k
			if n == 0 {
				return
			}
			n, k = k, n
			k = n - k
		} else {
			k = n - k
			if k == 1 {
				for i := p + n; i < p; i-- {
					data.Swap(i-1, i)
				}
			}
			q := p + n
			p = q - k
			for i := 0; i < n-k; i++ {
				p--
				q--
				data.Swap(p, q)
			}
			n %= k
			if n == 0 {
				return
			}
			n, k = k, n
		}
	}
}

func RotateSwap(data Interface, a, m, b int) {
	i := m - a
	if i == 0 {
		return
	}
	j := b - m
	if j == 0 {
		return
	}

	if i == j {
		for m < b {
			data.Swap(a, m)
			a++
			m++
		}
		return
	}

	p := a + i
	for i != j {
		if i > j {
			SwapBlock(data, p-i, p-i+j, p)
			i -= j
		} else {
			SwapBlock(data, p-i, p, p+j-i)
			j -= i
		}
	}
	SwapBlock(data, p-i, p, p)
}

// Swap two same length blocks uv to vu.
func SwapBlock(data Interface, first, end, second int) {
	n := end - first
	for i := 0; i < n; i++ {
		data.Swap(first+i, second+i)
	}
	/*
		for first < end {
			data.Swap(first, second)
			first++
			second++
		}
	*/
}

func floatHole(data Interface, first1, last1, first2, first3 int) {
	for first1 != last1 {
		data.Swap(first1, first2) // ABC --> BAC
		data.Swap(first2, first3) // BAC --> BCA
		first1++
		first2++
		first3++
	}
}

// See A. Kutzner.
func RotateSpeed(data Interface, a, m, b int) {
	deltaL := m - a
	if deltaL == 0 {
		return
	}
	deltaR := b - m
	if deltaR == 0 {
		return
	}

	for {
		// println(deltaL, deltaR)
		if deltaL < deltaR {
			if 2*deltaL > deltaR {
				m = a + (deltaR - deltaL)
				if a != m {
					// float_hole__(first, pivot, first + deltal, first + deltar);
					floatHole(data, a, m, a+deltaL, a+deltaR)
				}
				SwapBlock(data, m, a+deltaL, m+deltaR)
				a = m
				deltaR -= deltaL
				deltaL -= deltaR
			} else {
				m = a + deltaL
				floatHole(data, a, m, m, a+deltaR)
				a = m
				deltaR = deltaR - 2*deltaL
			}
		} else {
			if 2*deltaR > deltaL {
				m = a + deltaR
				floatHole(data, m+deltaR, m+deltaL, m, m-(deltaL-deltaR))
				SwapBlock(data, a, m-(deltaL-deltaR), a+deltaL)
				a = m
				deltaL -= deltaR
				deltaR -= deltaL
			} else {
				m = a + deltaR
				floatHole(data, a+deltaL, m+deltaL, m, a)
				a = m + deltaR
				deltaL = deltaL - 2*deltaR
			}
		}
		if deltaL <= 0 || deltaR <= 0 {
			break
		}
	}
}

func Gcd(m, n int) int {
	for n != 0 {
		m, n = n, m%n
	}
	return m
}

func RotateJuggling(data Interface, a, m, b int) {
	n := b - a
	k := m - a
	l := n - k
	if k == 0 {
		return
	} else if k == 1 {
		b--
		for a < b {
			data.Swap(a, a+1)
			a++
		}
		return
	} else if k == l {
		SwapBlock(data, a, m, m)
		return
	} else if l == 1 {
		for m > a {
			data.Swap(m, m-1)
			m--
		}
		return
	}

	d := Gcd(n, k)
	if k < l {
		for i := 0; i < d; i++ {
			p := a
			for j := 0; j < l/d; j++ {
				if p > a+l {
					data.Swap(p, p-l)
					p -= l
				}
				data.Swap(p, p+k)
				p += k
			}
			a++
		}
	} else {
		for i := 0; i < d; i++ {
			p := a
			for j := 0; j < k/d-1; j++ {
				if p < b-k {
					data.Swap(p, p+k)
					p += k
				}
				data.Swap(p, p-l)
				p -= l
			}
			a++
		}
	}
	return
}

// -------------------------------------------------------------------------
// Sym Merge Sort

// SymMergeSort performs stable in-place sorting of data.
func SymMergeSort(data Interface) {
	MergeSort(data, SymMerge, 7, 0, data.Len(), 0)
}

// SymMerge merges the two sorted subsequences data[first1,first2) and
// data[first2,last) using the SymMerge algorithm from: Pok-Son Kim and
// Arne Kutzner, "Stable Minimum Storage Merging by Symmetric Comparisons"
// (2004?)
// It needs $O(m \log(n/m+1))$ comparisons and its recursion depth is
// bound by $\lceil \log(n+m) \rceil$ where $m \leq n$ are the length
// of the both sequences to merge.
func SymMerge(data Interface, first1, first2, last int, depth int) int {
	if first1 >= first2 || first2 >= last {
		return depth
	}

	m := (first1 + last) / 2
	n := m + first2
	start := 0
	if first2 > m {
		start = n - last
		r, p := m, n-1
		for start < r {
			c := (start + r) / 2
			if !data.Less(p-c, c) {
				start = c + 1
			} else {
				r = c
			}
		}
	} else {
		start = first1
		r, p := first2, n-1
		for start < r {
			c := (start + r) / 2
			if !data.Less(p-c, c) {
				start = c + 1
			} else {
				r = c
			}
		}
	}
	end := n - start
	rotate(data, start, first2, end)
	ld := SymMerge(data, first1, start, m, depth+1)
	rd := SymMerge(data, m, end, last, depth+1)
	if ld > rd {
		return ld
	}
	return rd
}

func bsearch(data Interface, l, r, p int) int {
	for l < r {
		m := (l + r) / 2
		if !data.Less(p-m, m) {
			l = m + 1
		} else {
			r = m
		}
	}
	return l
}

// -------------------------------------------------------------------------
// Split Merge Sort

// SplitMergeSort performs stable in-place sorting of data.
func SplitMergeSort(data Interface) {
	MergeSort(data, SplitMerge, 7, 0, data.Len(), 0)
}

// SplitMerge merges the two sorted subsequences data[first1,first2) and
// data[first2,last) using the SymMerge algorithm from: Pok-Son Kim and
// Arne Kutzner, "A Simple Algorith for Stable Minimum Storage Merging",
// (2006?)
// It needs $O(m \log(n/m+1))$ comparisons and its recursion depth is
// bound by $m - 1$ where $m \leq n$ are the length of the both sequences
// to merge.
func SplitMerge(data Interface, first1, first2, last int, depth int) int {
	// u is in data[first1:first2], v is in data[first2:last]
	if first1 >= first2 || first2 >= last {
		return depth
	}

	l, r := first1, first2
	ld, rd := first2, last
	m, md := 0, 0

	for {
		if l < r {
			m = (l + r) / 2
		}
		if ld < rd {
			md = (ld + rd) / 2
		}

		if !data.Less(md, m) {
			l, rd = m+1, md
		} else {
			ld, r = md+1, m
		}

		if l >= r && ld >= rd {
			break
		}
	}

	rotate(data, r, first2, ld)
	ldep := SplitMerge(data, first1, r, r+rd-first2, depth+1)
	rdep := SplitMerge(data, l+ld-first2, ld, last, depth+1)
	if ldep > rdep {
		return ldep
	}
	return rdep

}

// -------------------------------------------------------------------------
// Inplace Stable Sort

// InplaceStableSort performs stable in-place sorting of data.
func InplaceStableSort(data Interface) {
	MergeSort(data, MergeWithoutBuffer, 12, 0, data.Len(), 0)
}

// merge the two sorted blocks data[first,middle) and data[middle,last).
// See http://gcc.gnu.org/onlinedocs/gcc-4.6.3/libstdc++/api/a01045_source.html#l03015
func MergeWithoutBuffer(data Interface, first, middle, last int, depth int) int {
	len1, len2 := middle-first, last-middle
	if len1 == 0 || len2 == 0 {
		return depth
	}
	if len1+len2 == 2 {
		if data.Less(middle, first) {
			data.Swap(first, middle)
		}
		return depth
	}

	first_cut := first
	second_cut := middle
	len11 := 0
	len22 := 0
	if len1 > len2 {
		len11 = len1 / 2
		first_cut += len11
		second_cut = lowerBound(data, middle, last, first_cut)
		len22 = second_cut - middle
	} else {
		len22 = len2 / 2
		second_cut += len22
		first_cut = upperBound(data, first, middle, second_cut)
		len11 = first_cut - first
	}
	rotate(data, first_cut, middle, second_cut)
	newMiddle := first_cut + (second_cut - middle)
	ld := MergeWithoutBuffer(data, first, first_cut, newMiddle, depth+1)
	rd := MergeWithoutBuffer(data, newMiddle, second_cut, last, depth+1)
	if ld > rd {
		return ld
	}
	return rd

}

func lowerBound(data Interface, first, last, cut int) int {
	len := last - first
	for len > 0 {
		half := len >> 1
		middle := first + len/2
		if data.Less(middle, cut) {
			first = middle
			first++
			len = len - half - 1
		} else {
			len = half
		}
	}
	return first
}

func upperBound(data Interface, first, last, cut int) int {
	len := last - first
	for len > 0 {
		half := len >> 1
		middle := first + half
		if data.Less(cut, middle) {
			len = half
		} else {
			first = middle
			first++
			len = len - half - 1
		}
	}
	return first
}
