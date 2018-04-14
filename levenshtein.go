// Levenshtein
// A module for fast computation of:
// - Levenshtein (edit) distance and edit sequence manipulation
// - string similarity
// - approximate median strings, and generally string averaging
// - string sequence and set similarity
//
// Levenshtein has a some overlap with difflib (SequenceMatcher).  It
// supports only strings, not arbitrary sequence types, but on the
// other hand it's much faster.
//
// It supports both normal and Unicode strings, but can't mix them, all
// arguments to a function (method) have to be of the same type (or its subclasses)
//

// The part contain:
// Basic stuff, Levenshtein distance.
// Other simple distances: Hamming, Jaro, Jaro-Winkler.
// Generalized medians, the greedy algorithm, and greedy improvements.
// Quick (voting) medians.
// Set medians.
// Set, sequence distances.

package levenshtein

import (
	"math"
	"sort"
)

const (
	INSERT_DELETE_COST = 1
	REPLACE_COST       = 1
	REPLACE_COST2      = 2
	STD_PREFIX_WEIGHT  = 0.1

	LEV_EPSILON  = 1e-14
	LEV_INFINITY = 1e100
)

// Edit opration type
// DON'T CHANGE! used ad arrays indices and the bits are occasionally used as flags
const (
	LEV_EDIT_KEEP = iota
	LEV_EDIT_REPLACE
	LEV_EDIT_INSERT
	LEV_EDIT_DELETE
	LEV_EDIT_LAST // sometimes returned when an error occurs
)

// ////////////////////////////////////////////////////
//
// Basic stuff, Levenshtein distance
//
// ////////////////////////////////////////////////////

type RuneSlice []rune

func (c RuneSlice) Len() int {
	return len(c)
}
func (c RuneSlice) Swap(i, j int) {
	c[i], c[j] = c[j], c[i]
}
func (c RuneSlice) Less(i, j int) bool {
	return c[i] < c[j]
}

func btoi(b bool) int {
	if b {
		return 1
	}
	return 0
}

func itob(i int) bool {
	return i != 0
}

func max2int(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min2int(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func min3int(a, b, c int) int {
	if a < b {
		if a < c {
			return a
		}
	} else {
		if b < c {
			return b
		}
	}
	return c
}

func min2float64(a, b float64) float64 {
	if a < b {
		return a
	}
	return b
}

func min3float64(a, b, c float64) float64 {
	if a < b {
		if a < c {
			return a
		}
	} else {
		if b < c {
			return b
		}
	}
	return c
}

// lev_edit_distance:
// xcost: If nonzero, the replace operation has weight 2, otherwise all
//         edit operations have equal weights of 1.
//
// Computes Levenshtein edit distance of two strings.
//
// Returns: The edit distance.
//
func lev_edit_distance(rs1, rs2 []rune, xcost int) int {
	len1 := len(rs1)
	len2 := len(rs2)
	var bi, ei int
	// strip common prefix
	for bi = 0; bi < len1 && bi < len2; bi++ {
		if rs1[bi] != rs2[bi] {
			break
		}
	}
	// strip common suffix
	for ei = 0; ei < len1-bi && ei < len2-bi; ei++ {
		if rs1[len1-1-ei] != rs2[len2-1-ei] {
			break
		}
	}
	len1 = len1 - ei - bi
	len2 = len2 - ei - bi
	// catch trivial cases
	if len1 == 0 {
		return len2
	}
	if len2 == 0 {
		return len1
	}
	rs1new := rs1[bi:(bi + len1)]
	rs2new := rs2[bi:(bi + len2)]
	// make the inner cycle (i.e. strings2) the longer one
	if len1 > len2 {
		len1, len2 = len2, len1
		rs1new, rs2new = rs2new, rs1new
	}
	// check len1 == 1 separately
	if len1 == 1 {
		var extra int
		for _, c := range rs2new {
			if c == rs1new[0] {
				extra = 1
				break
			}
		}
		if xcost != 0 {
			return len2 + 1 - 2*extra
		} else {
			return len2 - extra
		}
	}
	len1++
	len2++
	half := len1 >> 1
	// initalize first row
	row := make([]int, len2)
	endi := len2 - 1
	if xcost != 0 {
		for i := 0; i < len2; i++ {
			row[i] = i
		}
		for i := 1; i < len1; i++ {
			D := i
			x := i
			p2i := 0
			for pi := 1; pi <= endi; pi++ {
				if rs1new[i-1] == rs2new[p2i] {
					D--
					x = D
				} else {
					x++
				}
				p2i++
				D = row[pi]
				D++
				if x > D {
					x = D
				}
				row[pi] = x
			}
		}
	} else {
		// in this case we don't have to scan two corner triangles (of size len1/2)
		// in the matrix because no best path can go throught them. note this
		// breaks when len1 == len2 == 2 so the memchr() special case above is
		// necessary
		for i := 0; i < len2-half; i++ {
			row[i] = i
		}
		row[0] = len1 - half - 1
		for i := 1; i < len1; i++ {
			var pi, p2i int
			var D, x int
			// skip the upper triangle
			if i >= len1-half {
				offset := i - (len1 - half)
				p2i = offset
				pi = offset
				c3 := row[pi] + btoi(rs1new[i-1] != rs2new[p2i])
				pi++
				p2i++
				x = row[pi]
				x++
				D = x
				if x > c3 {
					x = c3
				}
				row[pi] = x
				pi++
			} else {
				pi = 1
				p2i = 0
				D, x = i, i
			}
			// skip the lower triangle
			if i <= half+1 {
				endi = len2 + i - half - 2
			}
			// main
			for pi <= endi {
				D--
				c3 := D + btoi(rs1new[i-1] != rs2new[p2i])
				p2i++
				x++
				if x > c3 {
					x = c3
				}
				D = row[pi]
				D++
				if x > D {
					x = D
				}
				row[pi] = x
				pi++
			}
			// lower triangle sentinel
			if i <= half {
				D--
				c3 := D + btoi(rs1new[i-1] != rs2new[p2i])
				p2i++
				x++
				if x > c3 {
					x = c3
				}
				row[pi] = x
			}
		}
	}
	return row[endi]
}

// Compute absolute Levenshtein distance of two strings.
//
// distance(string1, string2)
//
// Examples (it's hard to spell Levenshtein correctly):
//
// >>> distance('Levenshtein', 'Lenvinsten')
// 4
// >>> distance('Levenshtein', 'Levensthein')
// 2
// >>> distance('Levenshtein', 'Levenshten')
// 1
// >>> distance('Levenshtein', 'Levenshtein')
// 0
//
func Distance(s1, s2 string) int {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_edit_distance(rs1, rs2, 0)
}

// Compute similarity of two strings.
//
// ratio(string1, string2)
//
// The similarity is a number between 0 and 1, it's usually equal or
// somewhat higher than difflib.SequenceMatcher.ratio(), because it's
// based on real minimal edit distance.
//
// Examples:
//
// >>> ratio('Hello world!', 'Holly grail!')
// 0.583333
// >>> ratio('Brian', 'Jesus')
// 0.0
//
func Ratio(s1, s2 string) float64 {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	len1 := len(rs1)
	len2 := len(rs2)
	if len1 == 0 || len2 == 0 {
		if len1 == 0 && len2 == 0 {
			return 1.0
		}
		return 0.0
	}
	ldist := lev_edit_distance(rs1, rs2, 1)
	lensum := len1 + len2
	return float64(lensum-ldist) / float64(lensum)
}

// ////////////////////////////////////////////////////
//
// Other simple distances: Hamming, Jaro, Jaro-Winkler
//
// ////////////////////////////////////////////////////

// lev_hamming_distance:
//
// Computes Hamming distance of two strings.
//
// The strings must have the same length.
//
// Returns: The Hamming distance.
//
func lev_hamming_distance(rs1, rs2 []rune) int {
	len1 := len(rs1)
	len2 := len(rs2)
	if len1 > len2 {
		rs1, rs2 = rs2, rs1
		len1, len2 = len2, len1
	}
	dist := 0
	for i := 0; i < len1; i++ {
		if rs1[i] != rs2[i] {
			dist++
		}
	}
	return dist + (len2 - len1)
}

// lev_jaro_ratio:
//
// Computes Jaro string similarity metric of two strings.
//
// Returns: The Jaro metric of @string1 and @string2.
//
func lev_jaro_ratio(rs1, rs2 []rune) float64 {
	len1 := len(rs1)
	len2 := len(rs2)
	if len1 == 0 || len2 == 0 {
		if len1 == 0 && len2 == 0 {
			return 1.0
		}
		return 0.0
	}
	// make len1 always shorter (or equally long)
	if len1 > len2 {
		rs1, rs2 = rs2, rs1
		len1, len2 = len2, len1
	}
	// from Python-Levensthtein
	/* The literature about Jaro metric is confusing as the method of assigment
	 * of common characters is nowhere specified.  There are several possible
	 * deterministic mutual assignments of common characters of two strings.
	 * We use earliest-position method, which is however suboptimal (e.g., it
	 * gives two transpositions in jaro("Jaro", "Joaro") because of assigment
	 * of the first `o').  No reasonable algorithm for the optimal one is
	 * currently known to me. */
	idx := make([]int, len1)
	halflen := (len1 + 1) / 2
	var match int
	// the part with allowed range overlapping left
	for i := 0; i < halflen; i++ {
		for j := 0; j < i+halflen; j++ {
			if rs1[j] == rs2[i] && idx[j] == 0 {
				match++
				idx[j] = match
				break
			}
		}
	}
	// the part with allowed range overlapping right
	to := min2int(len1+halflen, len2)
	for i := halflen; i < to; i++ {
		for j := i - halflen; j < len1; j++ {
			if rs1[j] == rs2[i] && idx[j] == 0 {
				match++
				idx[j] = match
				break
			}
		}
	}
	if match == 0 {
		return 0.0
	}
	trans := 0
	i := 0
	for j := 0; j < len1; j++ {
		if idx[j] != 0 {
			i++
			if idx[j] != i {
				trans++
			}
		}
	}
	md := float64(match)
	return (md/float64(len1) + md/float64(len2) + 1.0 - float64(trans)/md/2.0) / 3.0
}

// lev_jaro_winkler_ratio:
// pfweight: Prefix weight, i.e., how much weight should be given to a
//            common prefix.
//
// Computes Jaro-Winkler string similarity metric of two strings.
//
// The formula is J + pfweight*P*(1-J), where J is Jaro metric and P is the
// length of common prefix.
//
// Returns: The Jaro-Winkler metric of @string1 and @string2.
//
func lev_jaro_winkler_ratio(rs1, rs2 []rune, pfweight float64) float64 {
	len1 := len(rs1)
	len2 := len(rs2)
	m := min2int(len1, len2)
	var p int
	for p = 0; p < m; p++ {
		if rs1[p] != rs2[p] {
			break
		}
	}
	r := lev_jaro_ratio(rs1, rs2)
	r = r + (1.0-r)*float64(p)*pfweight
	if r > 1.0 {
		return 1.0
	}
	return r
}

// Compute Hamming distance of two strings.
//
// hamming(string1, string2)
//
// The Hamming distance is simply the number of differing characters.
// That means the length of the strings must be the same.
//
// Examples:
// >>> hamming('Hello world!', 'Holly grail!')
// 7
// >>> hamming('Brian', 'Jesus')
// 5
//
func Hamming(s1, s2 string) int {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_hamming_distance(rs1, rs2)
}

// Compute Jaro string similarity metric of two strings.
//
// jaro(string1, string2)
//
// The Jaro string similarity metric is intended for short strings like
// personal last names.  It is 0 for completely different strings and
// 1 for identical strings.
//
// Examples:
// >>> jaro('Brian', 'Jesus')
// 0.0
// >>> jaro('Thorkel', 'Thorgier')
// 0.779761
// >>> jaro('Dinsdale', 'D')
// 0.708333
//
func Jaro(s1, s2 string) float64 {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_jaro_ratio(rs1, rs2)
}

// Compute Jaro string similarity metric of two strings.
//
// jaro_winkler(string1, string2[, prefix_weight])
//
// The Jaro-Winkler string similarity metric is a modification of Jaro
// metric giving more weight to common prefix, as spelling mistakes are
// more likely to occur near ends of words.
//
// The prefix weight is inverse value of common prefix length sufficient
// to consider the strings *identical*.  If no prefix weight is
// specified, 1/10 is used.
//
// Examples:
//
// >>> jaro_winkler('Brian', 'Jesus')
// 0.0
// >>> jaro_winkler('Thorkel', 'Thorgier')
// 0.867857
// >>> jaro_winkler('Dinsdale', 'D')
// 0.7375
// >>> jaro_winkler('Thorkel', 'Thorgier', 0.25)
// 1.0
//
func Jaro_winkler(s1, s2 string) float64 {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_jaro_winkler_ratio(rs1, rs2, STD_PREFIX_WEIGHT)
}

// ////////////////////////////////////////////////////
//
// Generalized medians, the greedy algorithm, and greedy improvements
//
// ////////////////////////////////////////////////////

// compute the sets of symbols each string contains, and the set of symbols
// in any of them (symset).  meanwhile, count how many different symbols
// there are (used below for symlist).
func make_symlist(rstrlis [][]rune) []rune {
	symlist := make([]rune, 0)
	symmap := make(map[rune]bool)
	for _, rs := range rstrlis {
		for _, c := range rs {
			if _, ok := symmap[c]; !ok {
				symmap[c] = true
				symlist = append(symlist, c)
			}
		}
	}
	rune_slice := RuneSlice(symlist)
	sort.Sort(rune_slice)
	copy(symlist, []rune(rune_slice))
	return symlist
}

// get (optional) weights, use 1 if none specified.
func extract_weightlist(wlist []float64, n int) []float64 {
	if wlist == nil {
		if n <= 0 {
			return nil
		}
		weights := make([]float64, n)
		for i, _ := range weights {
			weights[i] = 1.0
		}
		return weights
	}
	nw := len(wlist)
	wlen := max2int(nw, n)
	weights := make([]float64, wlen)
	for i := 0; i < nw; i++ {
		if wlist[i] < 0.0 {
			return nil
		}
		weights[i] = wlist[i]
	}
	for i := nw; i < wlen; i++ {
		weights[i] = 1.0
	}
	return weights
}

func get_rstrlis_maxLen(rstrlis [][]rune) int {
	var maxlen int
	for _, rs := range rstrlis {
		leni := len(rs)
		if maxlen < leni {
			maxlen = leni
		}
	}
	return maxlen
}

// lev_greedy_median:
// strlis: An array of strings, that may contain NUL characters.
// weights: The string weights (they behave exactly as multiplicities, though
//            any positive value is allowed, not just integers).
//
// Finds a generalized median string of strlis using the greedy algorithm.
//
// Note it's considerably more efficient to give a string with weight 2
// than to store two identical strings in strlis (with weights 1).
//
// Returns: The generalized median, as a newly allocated string.
//
func lev_greedy_median(rstrlis [][]rune, weights []float64) string {
	if len(rstrlis) == 0 {
		return ""
	}
	// find all symbols
	symlist := make_symlist(rstrlis)
	if len(symlist) == 0 {
		return ""
	}
	n := len(rstrlis)
	maxlen := get_rstrlis_maxLen(rstrlis)
	rows := make([][]int, n)
	for i, rs := range rstrlis {
		rows[i] = make([]int, len(rs)+1)
		for j, _ := range rows[i] {
			rows[i][j] = j
		}
	}
	// build up the approximate median string symbol by symbol
	// XXX: we actually exit on break below, but on the same condition
	stoplen := 2*maxlen + 1
	row := make([]int, stoplen+1)
	median := make([]rune, stoplen)
	mediandist := make([]float64, stoplen+1)
	mediandist[0] = 0.0
	for i := 0; i < n; i++ {
		mediandist[0] += float64(len(rstrlis[i])) * weights[i]
	}
	var minminsum float64
	for ln := 1; ln < stoplen+1; ln++ {
		minminsum = LEV_INFINITY
		row[0] = ln
		// iterate over all symbols we may want to add
		for _, symbol := range symlist {
			var totaldist, minsum float64
			// sum Levenshtein distances from all the strings, with given weights
			for i, rs := range rstrlis {
				leni := len(rs)
				minl := ln
				x := ln
				// compute how another row of Levenshtein matrix would look for median
				// string with this symbol added
				for k := 0; k < leni; k++ {
					D := rows[i][k] + btoi(symbol != rs[k])
					x++
					x = min3int(x, D, rows[i][k+1]+1)
					if x < minl {
						minl = x
					}
				}
				minsum += float64(minl) * weights[i]
				totaldist += float64(x) * weights[i]
			}
			// is this symbol better than all the others?
			if minsum < minminsum {
				minminsum = minsum
				mediandist[ln] = totaldist
				median[ln-1] = symbol
			}
		}
		// stop the iteration if we no longer need to recompute the matrix rows
		// or when we are over maxlen and adding more characters doesn't seem useful
		if ln == stoplen || (ln > maxlen && mediandist[ln] > mediandist[ln-1]) {
			stoplen = ln
			break
		}
		// now the best symbol is known, so recompute all matrix rows for this one
		symbol := median[ln-1]
		for i, rs := range rstrlis {
			leni := len(rs)
			oldrow := rows[i]
			// compute a row of Levenshtein matrix
			for k := 1; k < leni+1; k++ {
				c1 := oldrow[k] + 1
				c2 := row[k-1] + 1
				c3 := oldrow[k-1] + btoi(symbol != rs[k-1])
				row[k] = min3int(c3, c2, c1)
			}
			copy(oldrow, row[0:(leni+1)])
		}
	}
	// find the string with minimum total distance
	bestlen := 0
	for ln := 1; ln <= stoplen; ln++ {
		if mediandist[ln] < mediandist[bestlen] {
			bestlen = ln
		}
	}
	return string(median[0:bestlen])
}

// Knowing the distance matrices up to some row, finish the distance
// computations.
//
// rs1, len1 are already shortened.
//
func finish_distance_computations(len1 int, rs1 []rune, rstrlis [][]rune, weights []float64, rows [][]int, row []int) float64 {
	var distsum float64
	var offset int
	// catch trivia case
	if len1 == 0 {
		for j, rs := range rstrlis {
			distsum += float64(rows[j][len(rs)]) * weights[j]
		}
		return distsum
	}
	// iterate through the strings and sum the distances
	for j, rs := range rstrlis {
		leni := len(rs)
		ln := len1
		// strip common suffix (prefix CAN'T be stripped)
		for ln != 0 && leni != 0 && rs[leni-1] == rs1[ln-1] {
			ln--
			leni--
		}
		// catch trivial cases
		if ln == 0 {
			distsum += float64(rows[j][leni]) * weights[j]
			continue
		}
		offset = rows[j][0]
		if leni == 0 {
			distsum += float64(offset+ln) * weights[j]
			continue
		}
		// complete the matrix
		copy(row, rows[j][0:(leni+1)])
		endi := leni
		for i := 1; i <= ln; i++ {
			pi := 1
			rsi := 0
			x := i + offset
			D := x
			for pi <= endi {
				D--
				c3 := D + btoi(rs1[i-1] != rs[rsi])
				rsi++
				x++
				if x > c3 {
					x = c3
				}
				D = row[pi]
				D++
				if x > D {
					x = D
				}
				row[pi] = x
				pi++
			}
		}
		distsum += weights[j] * float64(row[endi])
	}
	return distsum
}

// lev_median_improve:
// s: The approximate generalized median string to be improved.
// strlis: An array of strings, that may contain NUL characters.
// weights: The string weights (they behave exactly as multiplicities, though
//          any positive value is allowed, not just integers).
//
// Tries to make @s a better generalized median string of @strings with
// small perturbations.
//
// It never returns a string with larger SOD than s; in the worst case, a
// string identical to s is returned.
//
// Returns: The improved generalized median, as a newly allocated string.
//
func lev_median_improve(rs []rune, rstrlis [][]rune, weights []float64) string {
	n := len(rstrlis)
	maxlen := get_rstrlis_maxLen(rstrlis)
	// find all symbols
	symlist := make_symlist(rstrlis)
	//symlist = []rune("Lehinstv")
	symlistlen := len(symlist)
	if symlistlen == 0 {
		return ""
	}
	rows := make([][]int, n)
	for i, tmp_rs := range rstrlis {
		leni := len(tmp_rs)
		rows[i] = make([]int, leni+1)
		for j := 0; j <= leni; j++ {
			rows[i][j] = j
		}
	}
	stoplen := 2*maxlen + 1
	row := make([]int, stoplen+2)
	// initialize median to given string
	median := make([]rune, stoplen+1)
	// we need -1st element for insertions a pos 0
	mi := 1
	medlen := len(rs)
	copy(median[mi:], rs)
	minminsum := finish_distance_computations(medlen,
		median[mi:(mi+medlen)], rstrlis, weights, rows, row)
	// sequentially try perturbations on all positions
	var operation, pos int
	var symbol, orig_symbol rune
	for pos <= medlen {
		symbol = median[mi+pos]
		operation = LEV_EDIT_KEEP
		// IF pos < medlength: FOREACH symbol: try to replace the symbol
		// at pos, if some lower the total distance, chooste the best
		if pos < medlen {
			orig_symbol = median[mi+pos]
			for _, tmp_symbol := range symlist {
				if tmp_symbol == orig_symbol {
					continue
				}
				median[mi+pos] = tmp_symbol
				sum := finish_distance_computations(medlen-pos,
					median[(mi+pos):(mi+medlen)], rstrlis, weights, rows, row)
				if sum < minminsum {
					minminsum = sum
					symbol = tmp_symbol
					operation = LEV_EDIT_REPLACE
				}
			}
			median[mi+pos] = orig_symbol
		}
		// FOREACH symbol: try to add it at pos, if some lower the total
		// distance, chooste the best (increase medlength)
		// We simulate insertion by replacing the character at pos-1
		orig_symbol = median[mi+pos-1]
		for _, tmp_symbol := range symlist {
			median[mi+pos-1] = tmp_symbol
			sum := finish_distance_computations(medlen-pos+1,
				median[(mi+pos-1):(mi+medlen)], rstrlis, weights, rows, row)
			if sum < minminsum {
				minminsum = sum
				symbol = tmp_symbol
				operation = LEV_EDIT_INSERT
			}
		}
		median[mi+pos-1] = orig_symbol
		// IF pos < medlength: try to delete the symbol at pos, if it lowers
		// the total distance remember it (decrease medlength)
		if pos < medlen {
			sum := finish_distance_computations(medlen-pos-1,
				median[(mi+pos+1):(mi+medlen)], rstrlis, weights, rows, row)
			if sum < minminsum {
				minminsum = sum
				operation = LEV_EDIT_DELETE
			}
		}
		// actually perform the operation
		switch operation {
		case LEV_EDIT_REPLACE:
			median[mi+pos] = symbol
		case LEV_EDIT_INSERT:
			copy(median[(mi+pos+1):], median[(mi+pos):(mi+medlen)])
			median[mi+pos] = symbol
			medlen++
		case LEV_EDIT_DELETE:
			copy(median[(mi+pos):], median[(mi+pos+1):(mi+medlen)])
			medlen--
		default:
			//
		}
		if !(medlen <= stoplen) {
			// error
		}
		// now the result is known, so recompute all matrix rows and move on
		if operation != LEV_EDIT_DELETE {
			symbol = median[mi+pos]
			row[0] = pos + 1
			for i, tmp_rs := range rstrlis {
				leni := len(tmp_rs)
				// compute a row of Levenshtein matrix
				for k := 1; k <= leni; k++ {
					c1 := rows[i][k] + 1
					c2 := row[k-1] + 1
					c3 := rows[i][k-1] + btoi(symbol != tmp_rs[k-1])
					row[k] = min3int(c1, c2, c3)
				}
				copy(rows[i], row[0:(leni+1)])
			}
			pos++
		}
	}
	return string(median[mi:(mi + medlen)])
}

func median_common(strlis []string, wlist []float64, engine func([][]rune, []float64) string) string {
	weights := extract_weightlist(wlist, len(strlis))
	n := len(strlis)
	rstrlis := make([][]rune, n)
	for i, s := range strlis {
		rstrlis[i] = []rune(s)
	}
	return engine(rstrlis, weights)
}

func median_improve_common(s string, strlis []string, wlist []float64, engine func([]rune, [][]rune, []float64) string) string {
	n := len(strlis)
	weights := extract_weightlist(wlist, n)
	rs := []rune(s)
	rstrlis := make([][]rune, n)
	for i, s := range strlis {
		rstrlis[i] = []rune(s)
	}
	return engine(rs, rstrlis, weights)
}

// Find an approximate generalized median string using greedy algorithm.
//
// median(string_sequence[, weight_sequence])
//
// You can optionally pass a weight for each string as the second
// argument.  The weights are interpreted as item multiplicities
// although any non-negative real numbers are accepted.  Use them to
// improve computation speed when strings often appear multiple times
// in the sequence.
//
// Examples:
//
// >>> median(['SpSm', 'mpamm', 'Spam', 'Spa', 'Sua', 'hSam'])
// 'Spam'
// >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
// ...          'Leveshtei', 'Lenshtein', 'Lvenstein',
// ...          'Levenhtin', 'evenshtei']\n" \
// >>> median(fixme)
// 'Levenshtein'
//
func Median(strlis []string, wlist []float64) string {
	return median_common(strlis, wlist, lev_greedy_median)
}

// Improve an approximate generalized median string by perturbations.
//
// median_improve(string, string_sequence[, weight_sequence])
//
// The first argument is the estimated generalized median string you
// want to improve, the others are the same as in median().  It returns
// a string with total distance less or equal to that of the given string.
//
// Note this is much slower than median().  Also note it performs only
// one improvement step, calling median_improve() again on the result
// may improve it further, though this is unlikely to happen unless the
// given string was not very similar to the actual generalized median.
//
// Examples:
//
// >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
// ...          'Leveshtei', 'Lenshtein', 'Lvenstein',
// ...          'Levenhtin', 'evenshtei']
// >>> median_improve('spam', fixme)
// 'enhtein'
// >>> median_improve(median_improve('spam', fixme), fixme)
// 'Levenshtein'
//
// It takes some work to change spam to Levenshtein.
//
func Median_improve(s string, strlis []string, wlist []float64) string {
	return median_improve_common(s, strlis, wlist, lev_median_improve)
}

// ////////////////////////////////////////////////////
//
// Quick (voting) medians
//
// ////////////////////////////////////////////////////

// compute the sets of symbols each string contains, and the set of symbols
// in any of them (symset).  meanwhile, count how many different symbols
// there are (used below for symlist).
// the symset is passed as an argument to avoid its allocation and
// deallocation when it's used in the caller too
//
func make_symlistset(rstrlis [][]rune) []rune {
	symset := make(map[rune]bool)
	symlist := make([]rune, 0)
	for _, rs := range rstrlis {
		for _, c := range rs {
			if _, ok := symset[c]; !ok {
				symset[c] = true
				symlist = append(symlist, c)
			}
		}
	}
	rune_slice := RuneSlice(symlist)
	sort.Sort(rune_slice)
	copy(symlist, []rune(rune_slice))
	return symlist
}

func lev_quick_median(rstrlis [][]rune, weights []float64) string {
	// first check whether the result would be an empty string
	// and compute resulting string length
	var ml, wl float64
	for i, rs := range rstrlis {
		ml += weights[i] * float64(len(rs))
		wl += weights[i]
	}
	if wl == 0.0 {
		return ""
	}
	ml = float64(math.Floor(float64(ml/wl + 0.499999)))
	ln := int(ml)
	if ln == 0 {
		return ""
	}
	symlist := make_symlistset(rstrlis)
	median := make([]rune, ln)
	symset := make(map[rune]float64)
	for j := 0; j < ln; j++ {
		for key, _ := range symset {
			delete(symset, key)
		}
		// let all strings vote
		for i, rs := range rstrlis {
			leni := len(rs)
			start := float64(leni) / ml * float64(j)
			end := start + float64(leni)/ml
			istart := int(math.Floor(float64(start)))
			iend := int(math.Ceil(float64(end)))
			// rounding errors can overflow the buffer
			if iend > leni {
				iend = leni
			}
			for k := istart + 1; k < iend; k++ {
				if _, ok := symset[rs[k]]; !ok {
					symset[rs[k]] = 0.0
				}
				symset[rs[k]] += weights[i]
			}
			if _, ok := symset[rs[istart]]; !ok {
				symset[rs[istart]] = 0.0
			}
			symset[rs[istart]] += weights[i] * float64(float64(1+istart)-start)
			symset[rs[iend-1]] -= weights[i] * float64(float64(iend)-end)
		}
		// find the elected symbol
		k := symlist[0]
		var wt float64
		for t_k, t_wt := range symset {
			if t_wt > wt {
				k, wt = t_k, t_wt
			}
		}
		median[j] = k
	}
	return string(median)
}

// Find a very approximate generalized median string, but fast.
//
// quickmedian(string[, weight_sequence])
//
// See median() for argument description.
//
// This method is somewhere between setmedian() and picking a random
// string from the set; both speedwise and quality-wise.
//
// Examples:
//
// >>> fixme = ['Levnhtein', 'Leveshein', 'Leenshten',
// ...          'Leveshtei', 'Lenshtein', 'Lvenstein',
// ...          'Levenhtin', 'evenshtei']
// >>> quickmedian(fixme)
// 'Levnshein'
//
func Quickmedian(strlis []string, wlist []float64) string {
	return median_common(strlis, wlist, lev_quick_median)
}

// ////////////////////////////////////////////////////
//
// Set medians
//
// ////////////////////////////////////////////////////

// lev_set_median_index:
// weights: The string weights (they behave exactly as multiplicities, though
//           any positive value is allowed, not just integers).
//
// Finds the median string of a string set @strings.
//
// Returns: An index in @strings pointing to the set median, -1 in case of failure.
//
func lev_set_median_index(rstrlis [][]rune, weights []float64) int {
	n := len(rstrlis)
	if n == 0 {
		return -1
	}
	if n == 1 {
		return 0
	}
	minidx := 0
	var mindist float64 = LEV_INFINITY
	distances := make([]int, n*(n-1)/2)
	for i, _ := range distances {
		distances[i] = -1
	}
	var j, d int
	var dist float64
	for i := 0; i < n; i++ {
		rs := rstrlis[i]
		dist = 0.0
		// below diagonal
		for j = 0; j < i && dist < mindist; j++ {
			dindex := (i-1)*(i-2)/2 + j
			if distances[dindex] >= 0 {
				d = distances[dindex]
			} else {
				d = lev_edit_distance(rstrlis[j], rs, 0)
			}
			if d < 0 {
				return -1
			}
			dist += weights[j] * float64(d)
		}
		// no need to compare item with itself
		j++
		// above diagonal
		for ; j < n && dist < mindist; j++ {
			dindex := (j-1)*(j-2)/2 + i
			distances[dindex] = lev_edit_distance(rstrlis[j], rs, 0)
			if distances[dindex] < 0 {
				return -1
			}
			dist += weights[j] * float64(distances[dindex])
		}
		if dist < mindist {
			mindist = dist
			minidx = i
		}
	}
	return minidx
}

// lev_set_median:
// weights: The string weights (they behave exactly as multiplicities, though
//           any positive value is allowed, not just integers).
//
// Finds the median string of a string set @strings.
//
// Returns: The set median as a newly allocate string. NULL in the case of failure.
//
func lev_set_median(rstrlis [][]rune, weights []float64) string {
	minidx := lev_set_median_index(rstrlis, weights)
	if minidx == -1 {
		return ""
	}
	rs := rstrlis[minidx]
	medlength := len(rs)
	if medlength == 0 {
		return ""
	}
	return string(rs[0:medlength])
}

// Find set median of a string set (passed as a sequence).
//
// setmedian(string[, weight_sequence])
//
// See median() for argument description.
//
// The returned string is always one of the strings in the sequence.
//
// Examples:
//
// >>> setmedian(['ehee', 'cceaes', 'chees', 'chreesc',
// ...            'chees', 'cheesee', 'cseese', 'chetese'])
// 'chees'
//
func Setmedian(strlis []string, wlist []float64) string {
	return median_common(strlis, wlist, lev_set_median)
}

// ////////////////////////////////////////////////////
//
// Set, sequence distances
//
// ////////////////////////////////////////////////////

// lev_edit_seq_distance:
// strlis1: An array of strings that may contain NUL characters.
// strlis2: An array of strings that may contain NUL characters.
//
// Finds the distance between string sequences @strings1 and @strings2.
//
// In other words, this is a double-Levenshtein algorithm.
//
// The cost of string replace operation is based on string similarity: it's
// zero for identical strings and 2 for completely unsimilar strings.
//
// Returns: The distance of the two sequences.
//
func lev_edit_seq_distance(rstrlis1, rstrlis2 [][]rune) float64 {
	n1 := len(rstrlis1)
	n2 := len(rstrlis2)
	var bi, ei int
	// strip common prefix
	for bi = 0; bi < n1 && bi < n2; bi++ {
		if !isSameRuneArr(rstrlis1[bi], rstrlis2[bi]) {
			break
		}
	}
	// strip common suffix
	for ei = 0; ei < n1-bi && ei < n2-bi; ei++ {
		if !isSameRuneArr(rstrlis1[n1-1-ei], rstrlis2[n2-1-ei]) {
			break
		}
	}
	n1 = n1 - ei - bi
	n2 = n2 - ei - bi
	if n1 == 0 {
		return float64(n2)
	}
	if n2 == 0 {
		return float64(n1)
	}
	rstrlis1new := rstrlis1[bi:(bi + n1)]
	rstrlis2new := rstrlis2[bi:(bi + n2)]
	// make the inner cycle (i.e. strings2) the longer one
	if n1 > n2 {
		n1, n2 = n2, n1
		rstrlis1new, rstrlis2new = rstrlis2new, rstrlis1new
	}
	n1++
	n2++
	// initalize first row
	row := make([]float64, n2)
	for i := 0; i < n2; i++ {
		row[i] = float64(i)
	}
	endi := n2 - 1
	// go through the matrix and compute the costs.  yes, this is an extremely
	// obfuscated version, but also extremely memory-conservative and relatively fast.
	for i := 1; i < n1; i++ {
		pi := 1
		rs1 := rstrlis1new[i-1]
		leni1 := len(rs1)
		s2i := 0
		var D, x, q float64
		D = float64(i) - 1.0
		x = float64(i)
		for pi <= endi {
			ln := leni1 + len(rstrlis2new[s2i])
			q = 0.0
			if ln == 0 {
				q = D
			} else {
				d := lev_edit_distance(rs1, rstrlis2new[s2i], 1)
				s2i++
				if d < 0 {
					return -1.0
				}
				q = D + 2.0/float64(ln)*float64(d)
			}
			x += 1.0
			D = row[pi]
			x = min3float64(x, q, D+1.0)
			row[pi] = x
			pi++
		}
	}
	return row[endi]
}

func isSameRuneArr(rs1, rs2 []rune) bool {
	if len(rs1) != len(rs2) {
		return false
	}
	n := len(rs1)
	for i := 0; i < n; i++ {
		if rs1[i] != rs2[i] {
			return false
		}
	}
	return true
}

// Munkers-Blackman algorithm.
//
func munkers_blackman(n1, n2 int, dists []float64) []int {
	covc := make([]int, n1)
	zstarc := make([]int, n1)
	covr := make([]int, n2)
	zstarr := make([]int, n2)
	zprimer := make([]int, n2)
	// step 0 (subtract minimal distance) and step 1 (find zeroes)
	for j := 0; j < n1; j++ {
		var minidx int
		min := dists[j]
		pi := n1
		for i := 1; i < n2; i++ {
			if min > dists[j+pi] {
				minidx = i
				min = dists[j+pi]
			}
			pi += n1
		}
		// subtract
		pi = 0
		for i := 0; i < n2; i++ {
			dists[j+pi] -= min
			if dists[j+pi] < LEV_EPSILON {
				dists[j+pi] = 0.0
			}
			pi += n1
		}
		// star the zero, if possible
		if zstarc[j] == 0 && zstarr[minidx] == 0 {
			zstarc[j] = minidx + 1
			zstarr[minidx] = j + 1
		} else {
			// otherwise try to find some other
			pi = 0
			for i := 0; i < n2; i++ {
				if i != minidx && dists[j+pi] == 0.0 && zstarc[j] == 0 && zstarr[i] == 0 {
					zstarc[j] = i + 1
					zstarr[i] = j + 1
					break
				}
				pi += n1
			}
		}
	}
	// main
	var min float64
	var gi, nc, jump_tag int
	for {
		// step 2 (cover columns containing z*)
		nc = 0
		for j := 0; j < n1; j++ {
			if zstarc[j] != 0 {
				covc[j] = 1
				nc++
			}
		}
		if nc == n1 {
			break
		}
		// step 3 (find uncovered zeroes)
		for {
			// step_3:
			// search uncovered matrix entries
			for j := 0; j < n1; j++ {
				jump_tag = 0
				pi := j
				if covc[j] != 0 {
					continue
				}
				for gi = 0; gi < n2; gi++ {
					if covr[gi] == 0 && dists[pi] == 0.0 {
						// when a zero is found, prime it
						zprimer[gi] = j + 1
						if zstarr[gi] != 0 {
							// if there's a z* in the same row,
							// uncover the column, cover the row and redo
							covr[gi] = 1
							covc[zstarr[gi]-1] = 0
							//goto step_3
							jump_tag = 3
							break
						} else {
							// if there's no z*,
							// we are at the end of our path an can convert z' to z*
							//goto step_4
							jump_tag = 4
							break
						}
					}
					pi += n1
				}
				if jump_tag == 3 || jump_tag == 4 {
					break
				}
			}
			if jump_tag == 3 {
				// goto step_3
				continue
			} else if jump_tag == 4 {
				// goto step_4
				break
			}
			// step 5 (new zero manufacturer)
			// we can't get here, unless no zero is found at all
			// find the smallest uncovered entry
			min = LEV_INFINITY
			for j := 0; j < n1; j++ {
				pi := j
				if covc[j] != 0 {
					continue
				}
				for i := 0; i < n2; i++ {
					if covr[i] == 0 && min > dists[pi] {
						min = dists[pi]
					}
					pi += n1
				}
			}
			// add it to all covered rows
			for i := 0; i < n2; i++ {
				pi := i * n1
				if covr[i] == 0 {
					continue
				}
				for j := 0; j < n1; j++ {
					dists[pi] += min
					pi++
				}
			}
			// subtract if from all uncovered columns
			for j := 0; j < n1; j++ {
				pi := j
				if covc[j] != 0 {
					continue
				}
				for i := 0; i < n2; i++ {
					dists[pi] -= min
					if dists[pi] < LEV_EPSILON {
						dists[pi] = 0.0
					}
					pi += n1
				}
			}
		}
		// step 4 (increment the number of z*)
		// i is the row number (we get it from step 3)
		// step_4:
		gi++
		for {
			x := gi
			gi--
			j := zprimer[gi] - 1
			zstarr[gi] = j + 1
			gi = zstarc[j]
			zstarc[j] = x
			if gi == 0 {
				break
			}
		}
		for j := 0; j < n2; j++ {
			zprimer[j] = 0
			covr[j] = 0
		}
		for j := 0; j < n1; j++ {
			covc[j] = 0
		}
	}
	for j := 0; j < n1; j++ {
		zstarc[j]--
	}
	return zstarc
}

// lev_set_distance:
// strlis1: An array of strings that may contain NUL characters.
// strlis2: An array of strings that may contain NUL characters.
//
// Finds the distance between string sets strlis1 and strlis2.
//
// The difference from lev_edit_seq_distance() is that order doesn't matter.
// The optimal association of @strings1 and @strings2 is found first and
// the similarity is computed for that.
//
// Uses sequential Munkers-Blackman algorithm.
//
// Returns: The distance of the two sets.
//
func lev_set_distance(rstrlis1, rstrlis2 [][]rune) float64 {
	n1 := len(rstrlis1)
	n2 := len(rstrlis2)
	if n1 == 0 {
		return float64(n2)
	}
	if n2 == 0 {
		return float64(n1)
	}
	// make the number of columns (n1) smaller than the number of rows
	if n1 > n2 {
		rstrlis1, rstrlis2 = rstrlis2, rstrlis1
		n1, n2 = n2, n1
	}
	// compute distances from each to each
	dists := make([]float64, n1*n2)
	di := 0
	for _, rs2 := range rstrlis2 {
		leni2 := len(rs2)
		for _, rs1 := range rstrlis1 {
			ln := leni2 + len(rs1)
			if ln > 0 {
				d := lev_edit_distance(rs2, rs1, 1)
				if d < 0 {
					return -1.0
				}
				dists[di] = float64(d) / float64(ln)
			} else {
				dists[di] = 0.0
			}
			di++
		}
	}
	// find the optimal mapping between the two sets
	mplis := munkers_blackman(n1, n2, dists)
	if mplis == nil {
		return -1.0
	}
	// sum the set distance
	sum := float64(n2 - n1)
	for j := 0; j < n1; j++ {
		i := mplis[j]
		rs1 := rstrlis1[j]
		rs2 := rstrlis2[i]
		ln := len(rs1) + len(rs2)
		if ln > 0 {
			d := lev_edit_distance(rs1, rs2, 1)
			if d < 0 {
				return -1.0
			}
			sum += 2.0 * float64(d) / float64(ln)
		}
	}
	return sum
}

func setseq_common(strlis1, strlis2 []string, engine func(rstrlis1, rstrlis2 [][]rune) float64) float64 {
	n1 := len(strlis1)
	n2 := len(strlis2)
	if n1+n2 == 0 {
		return 1.0
	}
	if n1 == 0 {
		return float64(n2)
	}
	if n2 == 0 {
		return float64(n1)
	}
	rstrlis1 := make([][]rune, n1)
	rstrlis2 := make([][]rune, n2)
	for i, s := range strlis1 {
		rstrlis1[i] = []rune(s)
	}
	for i, s := range strlis2 {
		rstrlis2[i] = []rune(s)
	}
	r := engine(rstrlis1, rstrlis2)
	return r
}

// Compute similarity ratio of two sequences of strings.
//
// seqratio(string_sequence1, string_sequence2)
//
// This is like ratio(), but for string sequences.  A kind of ratio()
// is used to to measure the cost of item change operation for the
// strings.
//
// Examples:
//
// >>> seqratio(['newspaper', 'litter bin', 'tinny', 'antelope'],
// ...          ['caribou', 'sausage', 'gorn', 'woody'])
// 0.21517857142857144
//
func Seqratio(strlis1, strlis2 []string) float64 {
	lensum := len(strlis1) + len(strlis2)
	if lensum == 0 {
		return 1.0
	}
	r := setseq_common(strlis1, strlis2, lev_edit_seq_distance)
	if r < 0 {
		return -1.0
	}
	return float64((float64(lensum) - r) / float64(lensum))
}

// Compute similarity ratio of two strings sets (passed as sequences).
//
// setratio(string_sequence1, string_sequence2)
//
// The best match between any strings in the first set and the second
// set (passed as sequences) is attempted.  I.e., the order doesn't
// matter here.
//
// Examples:
//
// >>> setratio(['newspaper', 'litter bin', 'tinny', 'antelope'],
// ...          ['caribou', 'sausage', 'gorn', 'woody'])
// 0.281845
//
// No, even reordering doesn't help the tinny words to match the woody ones.
//
func Setratio(strlis1, strlis2 []string) float64 {
	lensum := len(strlis1) + len(strlis2)
	if lensum == 0 {
		return 1.0
	}
	r := setseq_common(strlis1, strlis2, lev_set_distance)
	if r < 0 {
		return -1.0
	}
	return float64((float64(lensum) - r) / float64(lensum))
}
