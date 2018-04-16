// StringMatcher.py is an example SequenceMatcher-like class built
// on the top of Levenshtein. It misses some SequenceMatcher’s
// functionality, and has some extra OTOH.

package levenshtein

import (
	"fmt"
)

type StringMatcher interface {
	GetSeqs() (string, string)
	SetSeqs(string, string)
	SetSeq1(string)
	SetSeq2(s2 string)
	Distance() int
	Ratio() float64
	QuickRatio() float64
	RealQuickRatio() float64
	GetEditops() []*LevEditOp
	GetOpcodes() []*LevOpCode
	GetMatchingblocks() []*LevMatchingBlock
}

// A SequenceMatcher-like class built on the top of Levenshtein
type stringMatcher struct {
	s1, s2              string
	len1, len2          int
	_distance           int
	_distance_ok        bool
	_ratio              float64
	_ratio_ok           bool
	_editops            []*LevEditOp
	_editops_ok         bool
	_opcodes            []*LevOpCode
	_opcodes_ok         bool
	_matching_blocks    []*LevMatchingBlock
	_matching_blocks_ok bool
}

func NewStringMatcher(s1, s2 string) StringMatcher {
	m := &stringMatcher{
		s1:   s1,
		s2:   s2,
		len1: len([]rune(s1)),
		len2: len([]rune(s2)),
	}
	m.reset_cache()
	return m
}

func (m *stringMatcher) reset_cache() {
	m._distance_ok = false
	m._ratio_ok = false
	m._editops_ok = false
	m._editops = nil
	m._opcodes_ok = false
	m._opcodes = nil
	m._matching_blocks_ok = false
	m._matching_blocks = nil
}

func (m *stringMatcher) GetSeqs() (string, string) {
	return m.s1, m.s2
}

func (m *stringMatcher) SetSeqs(s1, s2 string) {
	m.s1, m.s2 = s1, s2
	m.len1, m.len2 = len([]rune(s1)), len([]rune(s2))
	m.reset_cache()
	return
}

func (m *stringMatcher) SetSeq1(s1 string) {
	m.s1 = s1
	m.len1 = len([]rune(s1))
	m.reset_cache()
	return
}

func (m *stringMatcher) SetSeq2(s2 string) {
	m.s2 = s2
	m.len2 = len([]rune(s2))
	m.reset_cache()
	return
}

func (m *stringMatcher) Distance() int {
	if !m._distance_ok {
		m._distance = Distance(m.s1, m.s2)
		m._distance_ok = true
	}
	return m._distance
}

func (m *stringMatcher) Ratio() float64 {
	if !m._ratio_ok {
		m._ratio = Ratio(m.s1, m.s2)
		m._ratio_ok = true
	}
	return m._ratio
}

func (m *stringMatcher) QuickRatio() float64 {
	return m.Ratio()
}

func (m *stringMatcher) RealQuickRatio() float64 {
	if m.len1+m.len2 == 0 {
		return 1.0
	}
	return 2.0 * float64(min2int(m.len1, m.len2)) / float64(m.len1+m.len2)
}

func (m *stringMatcher) GetEditops() []*LevEditOp {
	if !m._editops_ok {
		m._editops = Editops(m.s1, m.s2)
		m._editops_ok = true
	}
	fmt.Println(len(m._editops))
	return m._editops
}

func (m *stringMatcher) GetOpcodes() []*LevOpCode {
	if !m._opcodes_ok {
		m._opcodes = Opcodes(m.s1, m.s2)
		m._opcodes_ok = true
	}
	return m._opcodes
}

func (m *stringMatcher) GetMatchingblocks() []*LevMatchingBlock {
	if !m._matching_blocks_ok {
		m._matching_blocks = Matching_blocks_opcodes(m.GetOpcodes(), m.len1, m.len2)
		m._matching_blocks_ok = true
	}
	return m._matching_blocks
}
