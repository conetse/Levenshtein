// Editops and other difflib-like stuff.
//  - Subtract an edit subsequence from a sequence
//  - Find sequence of edit operations transforming one string to another.
//  - Invert the sense of an edit operation sequence.
//  - Apply a sequence of edit operations to a string.
//  - Find identical blocks in two strings.
//  - Subtract an edit subsequence from a sequence.

package levenshtein

import (
	"fmt"
	"strings"
)

var (
	OpTypeNames = map[int]string{
		LEV_EDIT_KEEP:    "equal",
		LEV_EDIT_REPLACE: "replace",
		LEV_EDIT_INSERT:  "insert",
		LEV_EDIT_DELETE:  "delete",
		LEV_EDIT_LAST:    "last",
	}
)

// Edit operation (atomic).
// This is the `native' atomic edit operation.  It differs from the difflib
// one's because it represents a change of one character, not a block.  And
// we usually don't care about LEV_EDIT_KEEP, though the functions can handle
// them.  The positions are interpreted as at the left edge of a character.
//
type LevEditOp struct {
	OpType int // editing operation type
	Spos   int // source block position
	Dpos   int // destination position
}

// Edit operation (difflib-compatible).
// This is not `native', but conversion functions exist.  These fields exactly
// correspond to the codeops() tuples fields (and this method is also the
// source of the silly OpCode name).  Sequences must span over complete
// strings, subsequences are simply edit sequences with more (or larger)
// LEV_EDIT_KEEP blocks.
//
type LevOpCode struct {
	OpType     int // editing operation type
	Sbeg, Send int // source block begin, end
	Dbeg, Dend int // destination block begin, end
}

// Matching block (difflib-compatible).
type LevMatchingBlock struct {
	Spos int
	Dpos int
	Len  int
}

func getOpTypeName(optype int) string {
	if name, ok := OpTypeNames[optype]; ok {
		return name
	}
	return "unknown"
}

func FormatLevEditOp(op *LevEditOp) string {
	if op == nil {
		return "nil"
	}
	return "LevEditOp" + fmt.Sprintf("LevEditOp(OpType: %s, Spos: %d, Dpos: %d)",
		getOpTypeName(op.OpType), op.Spos, op.Dpos)
}

func FormatLevEditOps(ops []*LevEditOp) string {
	strs := make([]string, len(ops))
	for i, op := range ops {
		if op == nil {
			strs[i] = "nil"
			continue
		}
		strs[i] = fmt.Sprintf("(%s, %d, %d)",
			getOpTypeName(op.OpType), op.Spos, op.Dpos)
	}
	return "LevEditOp[" + strings.Join(strs, ", ") + "]"
}

func FormatLevOpCode(bop *LevOpCode) string {
	if bop == nil {
		return "nil"
	}
	return fmt.Sprintf("LevOpCode(OpType: %s, Dpos: %d, Send: %d, Dbeg: %d, Dend: %d)",
		getOpTypeName(bop.OpType), bop.Sbeg, bop.Send, bop.Dbeg, bop.Dend)
}

func FormatLevOpCodes(ops []*LevOpCode) string {
	strs := make([]string, len(ops))
	for i, bop := range ops {
		if bop == nil {
			strs[i] = "nil"
			continue
		}
		strs[i] = fmt.Sprintf("(%s, %d, %d, %d, %d)",
			getOpTypeName(bop.OpType), bop.Sbeg, bop.Send, bop.Dbeg, bop.Dend)
	}
	return "LevOpCode[" + strings.Join(strs, ", ") + "]"
}

func FormatLevMatchingBlock(mb *LevMatchingBlock) string {
	if mb == nil {
		return "nil"
	}
	return fmt.Sprintf("LevMatchingBlock(Spos: %d, Dpos: %d, Len: %d)",
		mb.Spos, mb.Dpos, mb.Len)
}

func FormatLevMatchingBlocks(mbs []*LevMatchingBlock) string {
	strs := make([]string, len(mbs))
	for i, mb := range mbs {
		if mb == nil {
			strs[i] = "nil"
			continue
		}
		strs[i] = fmt.Sprintf("(%d, %d, %d)", mb.Spos, mb.Dpos, mb.Len)
	}
	return "LevMatchingBlock[" + strings.Join(strs, ", ") + "]"
}

// ////////////////////////////////////////////////////
//
// Editops and other difflib-like stuff.
//
// ////////////////////////////////////////////////////

// lev_editops_find:
//
// Find an optimal edit sequence from @string1 to @string2.
//
// When there's more than one optimal sequence, a one is arbitrarily (though
// deterministically) chosen.
//
// Returns: The optimal edit sequence. It is normalized, i.e., keep operations are not included.
//
func lev_editops_find(rs1, rs2 []rune) []*LevEditOp {
	len1 := len(rs1)
	len2 := len(rs2)
	// strip common prefix
	var bi, ei int
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
	// initalize first row and column
	nlen1 := len1 - ei - bi
	nlen2 := len2 - ei - bi
	rs1new := rs1[bi:]
	rs2new := rs2[bi:]
	nlen1++
	nlen2++
	matrix := make([]int, nlen1*nlen2)
	for i := 0; i < nlen2; i++ {
		matrix[i] = i
	}
	for i := 1; i < nlen1; i++ {
		matrix[nlen2*i] = i
	}

	// find the costs and fill the matrix
	for i := 1; i < nlen1; i++ {
		previ := (i - 1) * nlen2
		pi := i * nlen2
		end := pi + nlen2 - 1
		x := i
		rs2i := 0
		pi++
		for pi <= end {
			c3 := matrix[previ] + btoi(rs1new[i-1] != rs2new[rs2i])
			previ++
			rs2i++
			x++
			if x > c3 {
				x = c3
			}
			c3 = matrix[previ] + 1
			if x > c3 {
				x = c3
			}
			matrix[pi] = x
			pi++
		}
	}
	// find the way back
	return editops_from_cost_matrix(rs1new[0:nlen1], rs2new[0:nlen2], bi, bi, matrix)
}

// editops_from_cost_matrix:
// off1: The offset where the matrix starts from the start of string1.
// off2: The offset where the matrix starts from the start of string2.
// matrix: The cost matrix.
//
// Reconstructs the optimal edit sequence from the cost matrix matrix.
//
// The matrix is freed.
//
// Returns: The optimal edit sequence, as a newly allocated array of
//          elementary edit operations.
//
func editops_from_cost_matrix(rs1, rs2 []rune, off1, off2 int, matrix []int) []*LevEditOp {
	len1 := len(rs1)
	len2 := len(rs2)
	n := matrix[len1*len2-1]
	ops := make([]*LevEditOp, n)
	i := len1 - 1
	j := len2 - 1
	var dir int
	pi := len1*len2 - 1
	pos := n
	count := 0
	for i != 0 || j != 0 {
		count++
		if count > len1*len2 {
			// err
			break
		}
		// prefer contiuning in the same direction
		if dir < 0 && j != 0 && matrix[pi] == matrix[pi-1]+1 {
			pos--
			j--
			ops[pos] = &LevEditOp{
				OpType: LEV_EDIT_INSERT,
				Spos:   i + off1,
				Dpos:   j + off2,
			}
			pi--
			continue
		}
		if dir > 0 && i != 0 && matrix[pi] == matrix[pi-len2]+1 {
			pos--
			i--
			ops[pos] = &LevEditOp{
				OpType: LEV_EDIT_DELETE,
				Spos:   i + off1,
				Dpos:   j + off2,
			}
			pi -= len2
			continue
		}
		if i != 0 && j != 0 && matrix[pi] == matrix[pi-len2-1] && rs1[i-1] == rs2[j-1] {
			// don't be stupid like difflib, don't store LEV_EDIT_KEEP
			i--
			j--
			pi -= len2 + 1
			dir = 0
			continue
		}
		if i != 0 && j != 0 && matrix[pi] == matrix[pi-len2-1]+1 {
			pos--
			i--
			j--
			ops[pos] = &LevEditOp{
				OpType: LEV_EDIT_REPLACE,
				Spos:   i + off1,
				Dpos:   j + off2,
			}
			pi -= len2 + 1
			dir = 0
			continue
		}
		// we cant't turn directly from -1 to 1, in this case it would be better
		// to go diagonally, but check it (dir == 0)
		if dir == 0 && j != 0 && matrix[pi] == matrix[pi-1]+1 {
			pos--
			j--
			ops[pos] = &LevEditOp{
				OpType: LEV_EDIT_INSERT,
				Spos:   i + off1,
				Dpos:   j + off2,
			}
			pi--
			dir = -1
			continue
		}
		if dir == 0 && i != 0 && matrix[pi] == matrix[pi-len2]+1 {
			pos--
			i--
			ops[pos] = &LevEditOp{
				OpType: LEV_EDIT_DELETE,
				Spos:   i + off1,
				Dpos:   j + off2,
			}
			pi -= len2
			dir = 1
			continue
		}
	}
	return ops
}

// lev_editops_to_opcodes:
// ops: An array of elementary edit operations.
// len1: The length of the source string.
// len2: The length of the destination string.
//
// Converts elementary edit operations to difflib block operation codes.
//
// Note the string lengths are necessary since difflib doesn't allow omitting
// keep operations.
//
// Returns: The converted block operation codes.
//
func lev_editops_to_opcodes(ops []*LevEditOp, len1, len2 int) []*LevOpCode {
	// compute the number of blocks
	n := len(ops)
	i := n
	oi := 0
	var spos, dpos int
	bops := make([]*LevOpCode, 0)
	var bop *LevOpCode
	for i != 0 {
		// simply pretend there are no keep blocks
		for ops[oi].OpType == LEV_EDIT_KEEP {
			i--
			if i == 0 {
				break
			}
			oi++
		}
		if i == 0 {
			break
		}
		bop = &LevOpCode{
			Sbeg: spos,
			Dbeg: dpos,
		}
		if spos < ops[oi].Spos || dpos < ops[oi].Dpos {
			bop.OpType = LEV_EDIT_KEEP
			bop.Send, bop.Dend = ops[oi].Spos, ops[oi].Dpos
			bops = append(bops, bop)
			spos, dpos = ops[oi].Spos, ops[oi].Dpos
			bop = &LevOpCode{
				Sbeg: spos,
				Dbeg: dpos,
			}
		}
		tp := ops[oi].OpType
		switch tp {
		case LEV_EDIT_REPLACE:
			for {
				spos++
				dpos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		case LEV_EDIT_DELETE:
			for {
				spos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		case LEV_EDIT_INSERT:
			for {
				dpos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		default:
			//
		}
		bop.OpType = tp
		bop.Send, bop.Dend = spos, dpos
		bops = append(bops, bop)
	}
	return bops
}

// lev_editops_invert:
// ops: An array of elementary edit operations.
//
// Inverts the sense of ops. It is modified in place.
//
// In other words, ops becomes a valid partial edit for the original source
// and destination strings with their roles exchanged.
//
func lev_editops_invert(ops []*LevEditOp) []*LevEditOp {
	for i, _ := range ops {
		ops[i].Dpos, ops[i].Spos = ops[i].Spos, ops[i].Dpos
		if ops[i].OpType&2 != 0 {
			ops[i].OpType ^= 1
		}
	}
	return ops
}

// lev_opcodes_invert:
// bops: An array of difflib block edit operation codes.
//
// Inverts the sense of bops.  It is modified in place.
//
// In other words, bops becomes a partial edit for the original source
// and destination strings with their roles exchanged.
//
func lev_opcodes_invert(bops []*LevOpCode) []*LevOpCode {
	for i, _ := range bops {
		bops[i].Dbeg, bops[i].Sbeg = bops[i].Sbeg, bops[i].Dbeg
		bops[i].Dend, bops[i].Send = bops[i].Send, bops[i].Dend
		if bops[i].OpType&2 != 0 {
			bops[i].OpType ^= 1
		}
	}
	return bops
}

// lev_editops_apply:
// ops: An array of elementary edit operations.
//
// Applies a partial edit ops from string1 to string2.
//
// NB: ops is not checked for applicability.
//
// Returns: The result of the partial edit as a newly allocated string.
//
func lev_editops_apply(rs1, rs2 []rune, ops []*LevEditOp) string {
	len1 := len(rs1)
	n := len(ops)
	dst := make([]rune, n+len1)
	var dposi, sposi int
	for _, op := range ops {
		j := op.Spos - sposi + btoi(op.OpType == LEV_EDIT_KEEP)
		if j != 0 {
			copy(dst[dposi:], rs1[sposi:(sposi+j)])
			sposi += j
			dposi += j
		}
		switch op.OpType {
		case LEV_EDIT_DELETE:
			sposi++
		case LEV_EDIT_REPLACE:
			sposi++
			fallthrough
		case LEV_EDIT_INSERT:
			dst[dposi] = rs2[op.Dpos]
			dposi++
		default:
			//
		}
	}
	j := len1 - sposi
	if j != 0 {
		copy(dst[dposi:(dposi+j)], rs1[sposi:(sposi+j)])
		sposi += j
		dposi += j
	}
	return string(dst[0:dposi])
}

// lev_opcodes_apply:
// bops: An array of difflib block edit operation codes.
//
// Applies a sequence of difflib block operations to a string.
//
// NB: @bops is not checked for applicability.
//
// Returns: The result of the edit as a newly allocated string.
//
func lev_opcodes_apply(rs1, rs2 []rune, bops []*LevOpCode) string {
	len1 := len(rs1)
	len2 := len(rs2)
	dst := make([]rune, len1+len2)
	dposi := 0
	for _, bop := range bops {
		switch bop.OpType {
		case LEV_EDIT_INSERT, LEV_EDIT_REPLACE:
			ln := bop.Dend - bop.Dbeg
			copy(dst[dposi:(dposi+ln)], rs2[bop.Dbeg:bop.Dend])
		case LEV_EDIT_KEEP:
			ln := bop.Send - bop.Sbeg
			copy(dst[dposi:(dposi+ln)], rs1[bop.Sbeg:bop.Send])
		default:
			//
		}
		dposi += bop.Dend - bop.Dbeg
	}
	dst_ret := make([]rune, dposi)
	copy(dst_ret, dst[0:dposi])
	return string(dst_ret)
}

// lev_editops_matching_blocks:
// len1: The length of the source string.
// len2: The length of the destination string.
// ops: An array of elementary edit operations.
// nmblocks: Where the number of matching block should be stored.
//
// Computes the matching block corresponding to an optimal edit ops.
//
// Returns: The matching blocks as a newly allocated array.
//
func lev_editops_matching_blocks(ops []*LevEditOp, len1, len2 int) []*LevMatchingBlock {
	n := len(ops)
	mblocks := make([]*LevMatchingBlock, 0)
	if n == 0 {
		return mblocks
	}
	var spos, dpos, oi int
	for i := n; i != 0; {
		// simply pretend there are no keep blocks
		for ops[oi].OpType == LEV_EDIT_KEEP {
			i--
			if i == 0 {
				break
			}
			oi++
		}
		if i == 0 {
			break
		}
		if spos < ops[oi].Spos || dpos < ops[oi].Dpos {
			mb := &LevMatchingBlock{
				Spos: spos,
				Dpos: dpos,
				Len:  ops[oi].Spos - spos,
			}
			mblocks = append(mblocks, mb)
			spos, dpos = ops[oi].Spos, ops[oi].Dpos
		}
		tp := ops[oi].OpType
		switch tp {
		case LEV_EDIT_REPLACE:
			for {
				spos++
				dpos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		case LEV_EDIT_DELETE:
			for {
				spos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		case LEV_EDIT_INSERT:
			for {
				dpos++
				i--
				oi++
				if !(i != 0 && ops[oi].OpType == tp && spos == ops[oi].Spos && dpos == ops[oi].Dpos) {
					break
				}
			}
		default:
			//
		}
	}
	if spos < len1 || dpos < len2 {
		if len1-spos != len2-dpos {
			// err
		}
		mb := &LevMatchingBlock{
			Spos: spos,
			Dpos: dpos,
			Len:  len1 - spos,
		}
		mblocks = append(mblocks, mb)
	}
	return mblocks
}

// lev_opcodes_matching_blocks:
// len1: The length of the source string.
// len2: The length of the destination string.
// bops: An array of difflib block edit operation codes.
// nmblocks: Where the number of matching block should be stored.
//
// Computes the matching block corresponding to an optimal edit @bops.
//
// Returns: The matching blocks as a newly allocated array.
//
func lev_opcodes_matching_blocks(bops []*LevOpCode, len1, len2 int) []*LevMatchingBlock {
	nb := len(bops)
	mblocks := make([]*LevMatchingBlock, 0)
	if nb == 0 {
		return mblocks
	}
	bi := 0
	var mb *LevMatchingBlock
	for i := nb; i != 0; i-- {
		if bops[bi].OpType == LEV_EDIT_KEEP {
			mb = &LevMatchingBlock{
				Spos: bops[bi].Sbeg,
				Dpos: bops[bi].Dbeg,
			}
			for i != 0 && bops[bi].OpType == LEV_EDIT_KEEP {
				i--
				bi++
			}
			if i == 0 {
				mb.Len = len1 - mb.Spos
				mblocks = append(mblocks, mb)
				break
			}
			mb.Len = bops[bi].Sbeg - mb.Spos
			mblocks = append(mblocks, mb)
		}
		bi++
	}
	return mblocks
}

// lev_editops_subtract:
// ops: An array of elementary edit operations.
// sub: A subsequence (ordered subset) of @ops
//
// Subtracts a subsequence of elementary edit operations from a sequence.
//
// The remainder is a sequence that, applied to result of application of @sub,
// gives the same final result as application of ops to original string.
//
// Returns: A newly allocated array of normalized edit operations.
//          It is always normalized, i.e, without any keep operations.
//          On failure, nil is returned.
//
func lev_editops_subtract(ops, sub []*LevEditOp) []*LevEditOp {
	var nr, nn int
	for _, op := range ops {
		if op.OpType != LEV_EDIT_KEEP {
			nr++
		}
	}
	for _, op := range sub {
		if op.OpType != LEV_EDIT_KEEP {
			nn++
		}
	}
	if nn > nr {
		return nil
	}
	nr -= nn
	// subtract
	// we could simply return NULL when nr == 0, but then it would be possible
	// to subtract *any* sequence of the right length to get an empty sequence
	// -- clrealy incorrectly; so we have to scan the list to check
	if nr == 0 {
		return nil
	}
	rem := make([]*LevEditOp, nr)
	n := len(ops)
	nn = 0
	var j, shift int
	var shifts = []int{0, 0, 1, -1}
	for _, op := range sub {
		for j < n && (ops[j].Spos != op.Spos || ops[j].Dpos != op.Dpos || ops[j].OpType != op.OpType) {
			if ops[j].OpType != LEV_EDIT_KEEP {
				rem[nn] = &LevEditOp{
					OpType: ops[j].OpType,
					Spos:   ops[j].Spos + shift,
					Dpos:   ops[j].Dpos,
				}
				nn++
			}
			j++
		}
		if j == n {
			return nil
		}
		shift += shifts[int(op.OpType)]
		j++
	}
	for j < n {
		if ops[j].OpType != LEV_EDIT_KEEP {
			rem[nn] = &LevEditOp{
				OpType: ops[j].OpType,
				Spos:   ops[j].Spos + shift,
				Dpos:   ops[j].Dpos,
			}
			nn++
		}
		j++
	}
	if nn != nr {
		// err
	}
	return rem
}

// Find sequence of edit operations transforming one string to another.
// editops(source_string, destination_string)
// editops(edit_operations, source_length, destination_length)
//
// The result is a list of triples (operation, spos, dpos), where
// operation is one of 'equal', 'replace', 'insert', or 'delete';  spos
// and dpos are position of characters in the first (source) and the
// second (destination) strings.  These are operations on signle
// characters.  In fact the returned list doesn't contain the 'equal',
// but all the related functions accept both lists with and without 'equal's.
//
// Examples:
// >>> editops('spam', 'park')
// [('delete', 0, 0), ('insert', 3, 2), ('replace', 3, 3)]
// The alternate form editops(opcodes, source_string, destination_string)
// can be used for conversion from opcodes (5-tuples) to editops (you can
// pass strings or their lengths, it doesn't matter).
//
func Editops(s1, s2 string) []*LevEditOp {
	// convert: we were called (bops, s1, s2)
	// todo
	//
	// find editops: we were called (s1, s2)
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	ops := lev_editops_find(rs1, rs2)
	return ops
}

// Find sequence of edit operations transforming one string to another.
// opcodes(source_string, destination_string)
// opcodes(edit_operations, source_length, destination_length)
//
// The result is a list of 5-tuples with the same meaning as in
// SequenceMatcher's get_opcodes() output.  But since the algorithms
// differ, the actual sequences from Levenshtein and SequenceMatcher
// may differ too.
//
// Examples:
// >>> for x in opcodes('spam', 'park'):
// ...     print(x)
// ...
// ('delete', 0, 1, 0, 0)
// ('equal', 1, 3, 0, 2)
// ('insert', 3, 3, 2, 3)
// ('replace', 3, 4, 3, 4)
//
// The alternate form opcodes(editops, source_string, destination_string)
// can be used for conversion from editops (triples) to opcodes.
//
func Opcodes(s1, s2 string) []*LevOpCode {
	// convert: we were called (ops, s1, s2)
	// todo
	//
	// find opcodes: we were called (s1, s2)
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	ops := lev_editops_find(rs1, rs2)
	bops := lev_editops_to_opcodes(ops, len(rs1), len(rs2))
	return bops
}

// Invert the sense of an edit operation sequence.
//
// inverse(edit_operations)
//
// In other words, it returns a list of edit operations transforming the
// second (destination) string to the first (source).  It can be used
// with both editops and opcodes.
//
// Examples:
//
// >>> inverse(editops('spam', 'park'))
// [('insert', 0, 0), ('delete', 2, 3), ('replace', 3, 3)]
// >>> editops('park', 'spam')
// [('insert', 0, 0), ('delete', 2, 3), ('replace', 3, 3)]
//
func Inverse(ops []*LevEditOp) []*LevEditOp {
	return Inverse_editops(ops)
}

func Inverse_editops(ops []*LevEditOp) []*LevEditOp {
	return lev_editops_invert(ops)
}

func Inverse_opcodes(bops []*LevOpCode) []*LevOpCode {
	return lev_opcodes_invert(bops)
}

// Apply a sequence of edit operations to a string.
//
// apply_edit(edit_operations, source_string, destination_string)
//
// In the case of editops, the sequence can be arbitrary ordered subset
// of the edit sequence transforming source_string to destination_string.
//
// Examples:
//
// >>> e = editops('man', 'scotsman')
// >>> apply_edit(e, 'man', 'scotsman')
// 'scotsman'
// >>> apply_edit(e[:3], 'man', 'scotsman')
// 'scoman'
//
// The other form of edit operations, opcodes, is not very suitable for
// such a tricks, because it has to always span over complete strings,
// subsets can be created by carefully replacing blocks with 'equal'
// blocks, or by enlarging 'equal' block at the expense of other blocks
// and adjusting the other blocks accordingly.
//
// Examples:
// >>> a, b = 'spam and eggs', 'foo and bar'
// >>> e = opcodes(a, b)
// >>> apply_edit(inverse(e), b, a)
// 'spam and eggs'
// >>> e[4] = ('equal', 10, 13, 8, 11)
// >>> a, b, e
// >>> apply_edit(e, a, b)
// 'foo and ggs'
//
func Apply_edit(ops []*LevEditOp, s1, s2 string) string {
	return Apply_edit_editops(ops, s1, s2)
}

func Apply_edit_editops(ops []*LevEditOp, s1, s2 string) string {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_editops_apply(rs1, rs2, ops)
}

func Apply_edit_opcodes(bops []*LevOpCode, s1, s2 string) string {
	rs1 := []rune(s1)
	rs2 := []rune(s2)
	return lev_opcodes_apply(rs1, rs2, bops)
}

// Find identical blocks in two strings.
//
// matching_blocks(edit_operations, source_length, destination_length)
//
// The result is a list of triples with the same meaning as in
// SequenceMatcher's get_matching_blocks() output.  It can be used with
// both editops and opcodes.  The second and third arguments don't
// have to be actually strings, their lengths are enough.
//
// Examples:
//
// >>> a, b = 'spam', 'park'
// >>> matching_blocks(editops(a, b), a, b)
// [(1, 0, 2), (4, 4, 0)]
// >>> matching_blocks(editops(a, b), len(a), len(b))
// [(1, 0, 2), (4, 4, 0)]
//
// The last zero-length block is not an error, but it's there for
// compatibility with difflib which always emits it.
//
// One can join the matching blocks to get two identical strings:
// >>> a, b = 'dog kennels', 'mattresses'
// >>> mb = matching_blocks(editops(a,b), a, b)
// >>> ''.join([a[x[0]:x[0]+x[2]] for x in mb])
// 'ees'
// >>> ''.join([b[x[1]:x[1]+x[2]] for x in mb])
// 'ees'
//
func Matching_blocks(ops []*LevEditOp, len1, len2 int) []*LevMatchingBlock {
	if len(ops) == 0 || len1 < 0 || len2 < 0 {
		return nil
	}
	return Matching_blocks_editops(ops, len1, len2)
}

func Matching_blocks_editops(ops []*LevEditOp, len1, len2 int) []*LevMatchingBlock {
	if len(ops) == 0 || len1 < 0 || len2 < 0 {
		return nil
	}
	mbs := lev_editops_matching_blocks(ops, len1, len2)
	mb := &LevMatchingBlock{
		Spos: len1,
		Dpos: len2,
		Len:  0,
	}
	// The last zero-length block is not an error, but it's there for
	// compatibility with difflib which always emits it
	mbs = append(mbs, mb)
	return mbs
}

func Matching_blocks_opcodes(bops []*LevOpCode, len1, len2 int) []*LevMatchingBlock {
	if len(bops) == 0 || len1 < 0 || len2 < 0 {
		return nil
	}
	mbs := lev_opcodes_matching_blocks(bops, len1, len2)
	mb := &LevMatchingBlock{
		Spos: len1,
		Dpos: len2,
		Len:  0,
	}
	mbs = append(mbs, mb)
	return mbs
}

// Subtract an edit subsequence from a sequence.
//
// subtract_edit(edit_operations, subsequence)
//
// The result is equivalent to
// editops(apply_edit(subsequence, s1, s2), s2), except that is
// constructed directly from the edit operations.  That is, if you apply
// it to the result of subsequence application, you get the same final
// string as from application complete edit_operations.  It may be not
// identical, though (in amibuous cases, like insertion of a character
// next to the same character).
//
// The subtracted subsequence must be an ordered subset of edit_operations.
//
// Note this function does not accept difflib-style opcodes as no one in
// his right mind wants to create subsequences from them.
//
// Examples:
//
// >>> e = editops('man', 'scotsman')
// >>> e1 = e[:3]
// >>> bastard = apply_edit(e1, 'man', 'scotsman')
// >>> bastard
// 'scoman'
// >>> apply_edit(subtract_edit(e, e1), bastard, 'scotsman')
// 'scotsman'
//
func Subtract_edit(ops, osub []*LevEditOp) []*LevEditOp {
	return lev_editops_subtract(ops, osub)
}
