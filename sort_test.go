// Copyright 2009 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package sort_test

import (
	"fmt"
	"math"
	"math/rand"
	"reflect"
	. "sort"
	"strconv"
	"testing"
)

var ints = [...]int{74, 59, 238, -784, 9845, 959, 905, 0, 0, 42, 7586, -5467984, 7586}
var float64s = [...]float64{74.3, 59.0, math.Inf(1), 238.2, -784.0, 2.3, math.NaN(), math.NaN(), math.Inf(-1), 9845.768, -959.7485, 905, 7.8, 7.8}
var strings = [...]string{"", "Hello", "foo", "bar", "foo", "f00", "%*&^*&^&", "***"}

func TestSortIntSlice(t *testing.T) {
	data := ints
	a := IntSlice(data[0:])
	Sort(a)
	if !IsSorted(a) {
		t.Errorf("sorted %v", ints)
		t.Errorf("   got %v", data)
	}
}

func TestSortFloat64Slice(t *testing.T) {
	data := float64s
	a := Float64Slice(data[0:])
	Sort(a)
	if !IsSorted(a) {
		t.Errorf("sorted %v", float64s)
		t.Errorf("   got %v", data)
	}
}

func TestSortStringSlice(t *testing.T) {
	data := strings
	a := StringSlice(data[0:])
	Sort(a)
	if !IsSorted(a) {
		t.Errorf("sorted %v", strings)
		t.Errorf("   got %v", data)
	}
}

func TestInts(t *testing.T) {
	data := ints
	Ints(data[0:])
	if !IntsAreSorted(data[0:]) {
		t.Errorf("sorted %v", ints)
		t.Errorf("   got %v", data)
	}
}

func TestFloat64s(t *testing.T) {
	data := float64s
	Float64s(data[0:])
	if !Float64sAreSorted(data[0:]) {
		t.Errorf("sorted %v", float64s)
		t.Errorf("   got %v", data)
	}
}

func TestStrings(t *testing.T) {
	data := strings
	Strings(data[0:])
	if !StringsAreSorted(data[0:]) {
		t.Errorf("sorted %v", strings)
		t.Errorf("   got %v", data)
	}
}

func TestSortLarge_Random(t *testing.T) {
	n := 1000000
	if testing.Short() {
		n /= 100
	}
	data := make([]int, n)
	for i := 0; i < len(data); i++ {
		data[i] = rand.Intn(100)
	}
	if IntsAreSorted(data) {
		t.Fatalf("terrible rand.rand")
	}
	Ints(data)
	if !IntsAreSorted(data) {
		t.Errorf("sort didn't sort - 1M ints")
	}
}

func TestReverseSortIntSlice(t *testing.T) {
	data := ints
	data1 := ints
	a := IntSlice(data[0:])
	Sort(a)
	r := IntSlice(data1[0:])
	Sort(Reverse(r))
	for i := 0; i < len(data); i++ {
		if a[i] != r[len(data)-1-i] {
			t.Errorf("reverse sort didn't sort")
		}
		if i > len(data)/2 {
			break
		}
	}
}

func BenchmarkSortString1K(b *testing.B) {
	b.StopTimer()
	for i := 0; i < b.N; i++ {
		data := make([]string, 1<<10)
		for i := 0; i < len(data); i++ {
			data[i] = strconv.Itoa(i ^ 0x2cc)
		}
		b.StartTimer()
		Strings(data)
		b.StopTimer()
	}
}

func BenchmarkSortInt1K(b *testing.B) {
	b.StopTimer()
	for i := 0; i < b.N; i++ {
		data := make([]int, 1<<10)
		for i := 0; i < len(data); i++ {
			data[i] = i ^ 0x2cc
		}
		b.StartTimer()
		Ints(data)
		b.StopTimer()
	}
}

func BenchmarkSortInt64K(b *testing.B) {
	b.StopTimer()
	for i := 0; i < b.N; i++ {
		data := make([]int, 1<<16)
		for i := 0; i < len(data); i++ {
			data[i] = i ^ 0xcccc
		}
		b.StartTimer()
		Ints(data)
		b.StopTimer()
	}
}

const (
	_Sawtooth = iota
	_Rand
	_Stagger
	_Plateau
	_Shuffle
	_NDist
)

const (
	_Copy = iota
	_Reverse
	_ReverseFirstHalf
	_ReverseSecondHalf
	_Sorted
	_Dither
	_NMode
)

type testingData struct {
	desc        string
	t           *testing.T
	data        []int
	maxswap     int // number of swaps allowed
	ncmp, nswap int
}

func (d *testingData) Len() int { return len(d.data) }
func (d *testingData) Less(i, j int) bool {
	d.ncmp++
	return d.data[i] < d.data[j]
}
func (d *testingData) Swap(i, j int) {
	if d.nswap >= d.maxswap {
		d.t.Errorf("%s: used %d swaps sorting slice of %d", d.desc, d.nswap, len(d.data))
		d.t.FailNow()
	}
	d.nswap++
	d.data[i], d.data[j] = d.data[j], d.data[i]
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func lg(n int) int {
	i := 0
	for 1<<uint(i) < n {
		i++
	}
	return i
}

func testBentleyMcIlroy(t *testing.T, sort func(Interface)) {
	sizes := []int{100, 1023, 1024, 1025}
	if testing.Short() {
		sizes = []int{100, 127, 128, 129}
	}
	dists := []string{"sawtooth", "rand", "stagger", "plateau", "shuffle"}
	modes := []string{"copy", "reverse", "reverse1", "reverse2", "sort", "dither"}
	var tmp1, tmp2 [1025]int
	for _, n := range sizes {
		for m := 1; m < 2*n; m *= 2 {
			for dist := 0; dist < _NDist; dist++ {
				j := 0
				k := 1
				data := tmp1[0:n]
				for i := 0; i < n; i++ {
					switch dist {
					case _Sawtooth:
						data[i] = i % m
					case _Rand:
						data[i] = rand.Intn(m)
					case _Stagger:
						data[i] = (i*m + i) % n
					case _Plateau:
						data[i] = min(i, m)
					case _Shuffle:
						if rand.Intn(m) != 0 {
							j += 2
							data[i] = j
						} else {
							k += 2
							data[i] = k
						}
					}
				}

				mdata := tmp2[0:n]
				for mode := 0; mode < _NMode; mode++ {
					switch mode {
					case _Copy:
						for i := 0; i < n; i++ {
							mdata[i] = data[i]
						}
					case _Reverse:
						for i := 0; i < n; i++ {
							mdata[i] = data[n-i-1]
						}
					case _ReverseFirstHalf:
						for i := 0; i < n/2; i++ {
							mdata[i] = data[n/2-i-1]
						}
						for i := n / 2; i < n; i++ {
							mdata[i] = data[i]
						}
					case _ReverseSecondHalf:
						for i := 0; i < n/2; i++ {
							mdata[i] = data[i]
						}
						for i := n / 2; i < n; i++ {
							mdata[i] = data[n-(i-n/2)-1]
						}
					case _Sorted:
						for i := 0; i < n; i++ {
							mdata[i] = data[i]
						}
						// Ints is known to be correct
						// because mode Sort runs after mode _Copy.
						Ints(mdata)
					case _Dither:
						for i := 0; i < n; i++ {
							mdata[i] = data[i] + i%5
						}
					}

					desc := fmt.Sprintf("n=%d m=%d dist=%s mode=%s", n, m, dists[dist], modes[mode])
					d := &testingData{desc: desc, t: t, data: mdata[0:n], maxswap: n * lg(n) * 12 / 10}
					sort(d)
					// Uncomment if you are trying to improve the number of compares/swaps.
					// t.Logf("%s: ncmp=%d, nswp=%d", desc, d.ncmp, d.nswap)

					// If we were testing C qsort, we'd have to make a copy
					// of the slice and sort it ourselves and then compare
					// x against it, to ensure that qsort was only permuting
					// the data, not (for example) overwriting it with zeros.
					//
					// In go, we don't have to be so paranoid: since the only
					// mutating method Sort can call is TestingData.swap,
					// it suffices here just to check that the final slice is sorted.
					if !IntsAreSorted(mdata) {
						t.Errorf("%s: ints not sorted", desc)
						t.Errorf("\t%v", mdata)
						t.FailNow()
					}
				}
			}
		}
	}
}

func TestSortBM(t *testing.T) {
	testBentleyMcIlroy(t, Sort)
}

func TestHeapsortBM(t *testing.T) {
	testBentleyMcIlroy(t, Heapsort)
}

// This is based on the "antiquicksort" implementation by M. Douglas McIlroy.
// See http://www.cs.dartmouth.edu/~doug/mdmspe.pdf for more info.
type adversaryTestingData struct {
	data      []int
	keys      map[int]int
	candidate int
}

func (d *adversaryTestingData) Len() int { return len(d.data) }

func (d *adversaryTestingData) Less(i, j int) bool {
	if _, present := d.keys[i]; !present {
		if _, present := d.keys[j]; !present {
			if i == d.candidate {
				d.keys[i] = len(d.keys)
			} else {
				d.keys[j] = len(d.keys)
			}
		}
	}

	if _, present := d.keys[i]; !present {
		d.candidate = i
		return false
	}
	if _, present := d.keys[j]; !present {
		d.candidate = j
		return true
	}

	return d.keys[i] >= d.keys[j]
}

func (d *adversaryTestingData) Swap(i, j int) {
	d.data[i], d.data[j] = d.data[j], d.data[i]
}

func TestAdversary(t *testing.T) {
	const size = 100
	data := make([]int, size)
	for i := 0; i < size; i++ {
		data[i] = i
	}

	d := &adversaryTestingData{data, make(map[int]int), 0}
	Sort(d) // This should degenerate to heapsort.
}

// -------------------------------------------------------------------
// Stable sorting

func TestGcd(t *testing.T) {
	for i, tc := range []struct{ m, n, g int }{
		{15, 12, 3},
		{15, 16, 1},
		{12, 16, 4},
		{8, 16, 8},
		{16, 8, 8},
		{7, 13, 1},
		{7, 4, 1},
	} {
		if got := Gcd(tc.m, tc.n); got != tc.g {
			t.Errorf("%d: gcd(%d,%d)=%d, got %d", i, tc.m, tc.n, tc.g, got)
		}
	}
}

func testRotation(t *testing.T, algo func(data Interface, a, m, b int), name string) {
	data := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}

	algo(IntSlice(data), 0, 5, 10) // |u| = |v|
	if !reflect.DeepEqual(data, []int{5, 6, 7, 8, 9, 0, 1, 2, 3, 4}) {
		t.Fatalf("%s\nwant %v\n got %v", name, []int{5, 6, 7, 8, 9, 0, 1, 2, 3, 4}, data)
	}

	algo(IntSlice(data), 2, 6, 9) // |u| > |v|
	if !reflect.DeepEqual(data, []int{5, 6, 1, 2, 3, 7, 8, 9, 0, 4}) {
		t.Fatalf("%s\nwant %v\n got %v", name, []int{5, 6, 1, 2, 3, 7, 8, 9, 0, 4}, data)
	}

	algo(IntSlice(data), 1, 3, 8) // |u| < |v|
	if !reflect.DeepEqual(data, []int{5, 2, 3, 7, 8, 9, 6, 1, 0, 4}) {
		t.Fatalf("%s\nwant %v\n got %v", name, []int{5, 2, 3, 7, 8, 9, 6, 1, 0, 4}, data)
	}

	algo(IntSlice(data), 4, 5, 7) // |u| == 1
	if !reflect.DeepEqual(data, []int{5, 2, 3, 7, 9, 6, 8, 1, 0, 4}) {
		t.Fatalf("%s\nwant %v\n got %v", name, []int{5, 2, 3, 7, 9, 6, 8, 1, 0, 4}, data)
	}

	algo(IntSlice(data), 2, 5, 6) // |v| == 1
	if !reflect.DeepEqual(data, []int{5, 2, 6, 3, 7, 9, 8, 1, 0, 4}) {
		t.Fatalf("%s\nwant %v\n got %v", name, []int{5, 2, 6, 3, 7, 9, 8, 1, 0, 4}, data)
	}

	t.Logf("%s: ok", name)
}

func TestJuggling(t *testing.T) {
	data := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	RotateJuggling(IntSlice(data), 0, 2, 5)
	fmt.Printf("data = %v\n", data)
}

func TestRotation(t *testing.T) {
	testRotation(t, RotateSimple, "rotateSimple")
	testRotation(t, RotateOptimal, "rotateOptimal")
	testRotation(t, RotateSwap, "rotateSwap")
	testRotation(t, RotateSpeed, "rotateSpeed")
	testRotation(t, RotateJuggling, "rotateJuggling")
}

type intPairs []struct {
	a, b int
}

// IntPairs compare on a only.
func (d intPairs) Len() int           { return len(d) }
func (d intPairs) Less(i, j int) bool { return d[i].a < d[j].a }
func (d intPairs) Swap(i, j int)      { d[i], d[j] = d[j], d[i] }
func (d intPairs) initB() {
	for i := range d {
		d[i].b = i
	}
}

// check if the bs to one a are in order.
func (d intPairs) inOrder() bool {
	lastA, lastB := -1, 0
	for i := 0; i < len(d); i++ {
		if lastA != d[i].a {
			lastA = d[i].a
			lastB = d[i].b
			continue
		}
		if d[i].b <= lastB {
			return false
		}
		lastB = d[i].b
	}
	return true
}

//
// Ints
//
func testInts(t *testing.T, algo func(data Interface), name string) {
	data := ints
	algo(IntSlice(data[0:]))
	if !IntsAreSorted(data[0:]) {
		t.Errorf("%s\nsorted %v\n   got %v", name, ints, data)
	}
}

func TestSymMergeSortInts(t *testing.T)      { testInts(t, SymMergeSort, "SymMergeSort") }
func TestSplitMergeSortInts(t *testing.T)    { testInts(t, SplitMergeSort, "SplitMergeSort") }
func TestInplaceStableSortInts(t *testing.T) { testInts(t, InplaceStableSort, "InplaceStableSort") }

//
// Random
//
func testRandom(t *testing.T, algo func(data Interface), name string) {
	n := 1000000
	if testing.Short() {
		n /= 100
	}
	data := make([]int, n)
	for i := 0; i < len(data); i++ {
		data[i] = rand.Intn(100)
	}
	if IntsAreSorted(data) {
		t.Fatalf("terrible rand.rand")
	}
	SymMergeSort(IntSlice(data))
	if !IntsAreSorted(data) {
		t.Errorf("%s sort didn't sort - 1M ints", name)
	}
}

func TestSymMergeSortRandom(t *testing.T)      { testRandom(t, SymMergeSort, "SymMergeSort") }
func TestSplitMergeSortRandom(t *testing.T)    { testRandom(t, SplitMergeSort, "SplitMergeSort") }
func TestInplaceStableSortRandom(t *testing.T) { testRandom(t, InplaceStableSort, "InplaceStableSort") }

//
// Stability
//
func testStability(t *testing.T, algo func(data Interface), name string) {
	n, m := 1000000, 1000
	if testing.Short() {
		n, m = 10000, 100
	}
	data := make(intPairs, n)

	// random distribution
	for i := 0; i < len(data); i++ {
		data[i].a = rand.Intn(m)
	}
	if IsSorted(data) {
		t.Fatalf("terrible rand.rand")
	}
	data.initB()
	algo(data)
	if !IsSorted(data) {
		t.Errorf("%s didn't sort - %d ints", name, n)
	}
	if !data.inOrder() {
		t.Errorf("%s wasn't stable - %d ints", name, n)
	}

	// already sorted
	data.initB()
	algo(data)
	if !IsSorted(data) {
		t.Errorf("%s didn't sort - %d ints", name, n)
	}
	if !data.inOrder() {
		t.Errorf("%s wasn't stable - %d ints", name, n)
	}

	// sorted reversed
	for i := 0; i < len(data); i++ {
		data[i].a = len(data) - i
	}
	data.initB()
	algo(data)
	if !IsSorted(data) {
		t.Errorf("%s didn't sort - %d ints", name, n)
	}
	if !data.inOrder() {
		t.Errorf("sym merge wasn't stable - 1M ints")
	}
}

func TestSymMergeSortStability(t *testing.T) {
	testStability(t, SymMergeSort, "SymMergeSort")
}
func TestSplitMergeSortStability(t *testing.T) {
	testStability(t, SplitMergeSort, "SplitMergeSort")
}
func TestInplaceStableSortStability(t *testing.T) {
	testStability(t, InplaceStableSort, "InplaceStableSort")
}

//
// count swaps, compares and recursion depth
//
func countOps(t *testing.T, merge func(Interface, int, int, int, int) int, name string) {
	sizes := []int{1 << 14, 1 << 16, 1 << 18, 1 << 20, 1<<22}
	for _, dist := range []string{"RND", "rnd", "inc", "dec", "saw", "was"} {
		for _, n := range sizes {
			td := testingData{
				desc:    name,
				t:       t,
				data:    make([]int, n),
				maxswap: 1 << 32,
			}

			for i := 0; i < n; i++ {
				switch dist {
				case "rnd":
					td.data[i] = rand.Intn(n/100)
				case "RND":
					td.data[i] = rand.Intn(10 * n)
				case "inc":
					td.data[i] = i
				case "dec":
					td.data[i] = n - i
				case "saw":
					td.data[i] = i % 256
				case "was":
					td.data[i] = n - i%256
				default:
					panic(dist)
				}
			}
			rd := MergeSort(&td, merge, 0, 0, n, 1)
			s, c := float64(td.nswap), float64(td.ncmp)
			logN := math.Log(float64(n))
			N := float64(n)
			t.Logf("%s %7d %s: %9d swaps, %9d cmps, %2d recd; s/NlogN=%.2f, s/Nlog^2N=%.2f, c/N=%5.2f, c/NlogN=%.2f",
				name, n, dist, td.nswap, td.ncmp, rd, s/(N*logN), 10*s/(N*logN*logN), c/N, c/(N*logN))
		}
	}
}

func TestSplitMergeOps(t *testing.T)    { countOps(t, SplitMerge, "SplitMergeSort") }
func TestSymMergeOps(t *testing.T) { countOps(t, SymMerge, "SymMergeSort") }
func TestInplaceStableOps(t *testing.T) { countOps(t, MergeWithoutBuffer, "InplaceStableSort") }

//
// Benchmarks
//

func benchmarkInt(b *testing.B, size int, algo func(Interface), name string) {
	b.StopTimer()
	var xor int
	switch size {
	case 1 << 10:
		xor = 0x2cc
	case 1 << 16:
		xor = 0xcccc
	case 1 << 22:
		xor = 0x2ccccc
	default:
		panic(size)
	}
	data := make([]int, size)
	for i := 0; i < b.N; i++ {
		for i := 0; i < len(data); i++ {
			data[i] = i ^ xor
		}
		b.StartTimer()
		algo(IntSlice(data))
		b.StopTimer()
		if !IsSorted(IntSlice(data)) {
			b.Errorf("%s did not sort %d ints", name, size)
		}
	}

}

func benchmarkString(b *testing.B, algo func(Interface), name string) {
	b.StopTimer()
	data := make([]string, 1<<10)
	for i := 0; i < b.N; i++ {
		for i := 0; i < len(data); i++ {
			data[i] = strconv.Itoa(i ^ 0x2cc)
		}
		b.StartTimer()
		algo(StringSlice(data))
		b.StopTimer()
		if !IsSorted(StringSlice(data)) {
			b.Errorf("%s did not sort strings", name)
		}
	}
}

func BenchmarkSymMergeSortInt1K(b *testing.B) {
	benchmarkInt(b, 1<<10, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortInt1K(b *testing.B) {
	benchmarkInt(b, 1<<10, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortInt1K(b *testing.B) {
	benchmarkInt(b, 1<<10, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortInt1K(b *testing.B) {
	benchmarkInt(b, 1<<10, Sort, "Sort")
}

func BenchmarkSymMergeSortInt64K(b *testing.B) {
	benchmarkInt(b, 1<<16, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortInt64K(b *testing.B) {
	benchmarkInt(b, 1<<16, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortInt64K(b *testing.B) {
	benchmarkInt(b, 1<<16, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortInt64K(b *testing.B) {
	benchmarkInt(b, 1<<16, Sort, "Sort")
}

func BenchmarkSymMergeSortInt4M(b *testing.B) {
	benchmarkInt(b, 1<<22, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortInt4M(b *testing.B) {
	benchmarkInt(b, 1<<22, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortInt4M(b *testing.B) {
	benchmarkInt(b, 1<<22, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortInt4M(b *testing.B) { benchmarkInt(b, 1<<22, Sort, "Sort") }

func BenchmarkSymMergeSortString1K(b *testing.B) {
	benchmarkString(b, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortString1K(b *testing.B) {
	benchmarkString(b, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortString1K(b *testing.B) {
	benchmarkString(b, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortString1K(b *testing.B) {
	benchmarkString(b, Sort, "Sort")
}

func benchmarkSorted(b *testing.B, size int, reversed bool, algo func(Interface), name string) {
	b.StopTimer()
	data := make([]int, size)
	for i := 0; i < b.N; i++ {
		for i := 0; i < len(data); i++ {
			if reversed {
				data[i] = size - i
			} else {
				data[i] = i
			}
		}
		b.StartTimer()
		algo(IntSlice(data))
		b.StopTimer()
		if !IsSorted(IntSlice(data)) {
			b.Errorf("%s bad sorting", name, size)
		}
	}
}

func BenchmarkSymMergeSortedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, false, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, false, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, false, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, false, Sort, "Sort")
}
func BenchmarkSymMergeReversedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, true, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeReversedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, true, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableReversedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, true, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortReversedInt64K(b *testing.B) {
	benchmarkSorted(b, 1<<16, true, Sort, "Sort")
}

func benchmarkRandom(b *testing.B, subset bool, algo func(data Interface), name string) {
	b.StopTimer()
	n := 1000000
	m := 10 * n
	if subset {
		m = n / 1000
	}
	data := make([]int, n)
	for j := 0; j < b.N; j++ {
		for i := 0; i < len(data); i++ {
			data[i] = rand.Intn(m)
		}
		b.StartTimer()
		algo(IntSlice(data))
		b.StopTimer()
	}
}
func BenchmarkSymMergeSortUniqe(b *testing.B) {
	benchmarkRandom(b, false, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortUniqe(b *testing.B) {
	benchmarkRandom(b, false, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortUniqe(b *testing.B) {
	benchmarkRandom(b, false, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortUniqe(b *testing.B) {
	benchmarkRandom(b, false, Sort, "Sort")
}
func BenchmarkSymMergeSortSubset(b *testing.B) {
	benchmarkRandom(b, true, SymMergeSort, "SymMergeSort")
}
func BenchmarkSplitMergeSortSubset(b *testing.B) {
	benchmarkRandom(b, true, SplitMergeSort, "SplitMergeSort")
}
func BenchmarkInplaceStableSortSubset(b *testing.B) {
	benchmarkRandom(b, true, InplaceStableSort, "InplaceStableSort")
}
func BenchmarkStandardSortSubset(b *testing.B) {
	benchmarkRandom(b, true, Sort, "Sort")
}