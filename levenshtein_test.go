package levenshtein

import (
	"fmt"
	"testing"
)

func TestLevenshtein(t *testing.T) {
	var s1, s2 string
	var hdist int
	var ratio, dist float64
	s1 = "Levenshtein"
	s2 = "Lenvinsten"
	//s1 = "你好一二三四五"
	//s2 = "世界abcde"
	hdist = Distance(s1, s2)
	// 4
	fmt.Println(hdist)
	s1 = "Hello world!"
	s2 = "Holly grail!"
	ratio = Ratio(s1, s2)
	// 0.5833333333333334
	fmt.Println(ratio)
	//
	hdist = Hamming(s1, s2)
	// 7
	fmt.Println(hdist)
	s1 = "Thorkel"
	s2 = "Thorgier"
	dist = Jaro(s1, s2)
	// 0.7797619047619048
	fmt.Println(dist)
	dist = Jaro_winkler(s1, s2)
	// 0.8678571428571429
	fmt.Println(dist)

	var medstr string
	var strlis, strlis1, strlis2, fixme []string
	strlis = []string{"SpSm", "mpamm", "Spam", "Spa", "Sua", "hSam"}
	fixme = []string{"Levnhtein", "Leveshein", "Leenshten", "Leveshtei", "Lenshtein", "Lvenstein", "Levenhtin", "evenshtei"}
	medstr = Median(strlis, nil)
	// Spam
	fmt.Println(medstr)
	medstr = Median(fixme, nil)
	// Levenshtein
	fmt.Println(medstr)
	s1 = "spam"
	medstr = Median_improve(s1, fixme, nil)
	// enhtein
	fmt.Println(medstr)
	medstr = Quickmedian(fixme, nil)
	// Levnshein
	fmt.Println(medstr)
	strlis = []string{"ehee", "cceaes", "chees", "chreesc", "chees", "cheesee", "cseese", "chetese"}
	medstr = Setmedian(strlis, nil)
	// chees
	fmt.Println(medstr)
	strlis1 = []string{"newspaper", "litter bin", "tinny", "antelope"}
	strlis2 = []string{"caribou", "sausage", "gorn", "woody"}
	ratio = Seqratio(strlis1, strlis2)
	// 0.21517857142857144
	fmt.Println(ratio)
	ratio = Setratio(strlis1, strlis2)
	// 0.281845
	fmt.Println(ratio)
}

func TestDifflib(t *testing.T) {
	var a, b, c, bastard string
	var ops []*LevEditOp
	a = "spam"
	b = "park"
	ops = Editops(a, b)
	// [('delete', 0, 0), ('insert', 3, 2), ('replace', 3, 3)]
	fmt.Println(FormatLevEditOps(ops))
	//
	a = "Levenshtein"
	b = "Lenvinsten"
	ops2 := Editops(a, b)
	// [('insert', 2, 2), ('replace', 3, 4), ('delete', 6, 7), ('delete', 9, 9)]
	fmt.Println(FormatLevEditOps(ops2))
	//
	a = "spam"
	b = "park"
	bops := Opcodes(a, b)
	//[('delete', 0, 1, 0, 0),
	// ('equal', 1, 3, 0, 2),
	// ('insert', 3, 3, 2, 3),
	// ('replace', 3, 4, 3, 4)]
	fmt.Println(FormatLevOpCodes(bops))
	ops1 := Editops(b, a)
	Inverse(ops1)
	// [('insert', 0, 0), ('delete', 2, 3), ('replace', 3, 3)]
	fmt.Println(FormatLevEditOps(ops))
	bops2 := Opcodes(b, a)
	Inverse_opcodes(bops2)
	//
	fmt.Println(FormatLevOpCodes(bops2))
	a = "spam and eggs"
	b = "foo and bar"
	ed := Editops(a, b)
	c = Apply_edit(Inverse_editops(ed), b, a)
	// spam and eggs
	fmt.Println(c)
	e := Opcodes(a, b)
	c = Apply_edit_opcodes(Inverse_opcodes(e), b, a)
	// spam and eggs
	fmt.Println(c)
	a = "spam"
	b = "park"
	ms := Matching_blocks(Editops(a, b), len([]rune(a)), len([]rune(b)))
	// [(1, 0, 2), (4, 4, 0)]
	fmt.Println(FormatLevMatchingBlocks(ms))
	ms = Matching_blocks_opcodes(Opcodes(a, b), len([]rune(a)), len([]rune(b)))
	fmt.Println(FormatLevMatchingBlocks(ms))
	a = "man"
	b = "scotsman"
	ee := Editops(a, b)
	fmt.Println(FormatLevEditOps(ee))
	e1 := ee[:3]
	fmt.Println(FormatLevEditOps(e1))
	bastard = Apply_edit(e1, a, b)
	// scoman
	fmt.Println(bastard)
	e2 := Subtract_edit(ee, e1)
	fmt.Println(FormatLevEditOps(e2))
	c = Apply_edit(Subtract_edit(ee, e1), bastard, b)
	// scotsman
	fmt.Println(c)
}

func TestStringMatcher(t *testing.T) {
	var a, b string
	var r float64
	a = "Levenshtein"
	b = "Lenvinsten"
	m := NewStringMatcher(a, b)
	d := m.Distance()
	// 4
	fmt.Println(d)
	m.SetSeqs("Hello world!", "Holly grail!")
	r = m.Ratio()
	// 0.5833333333333334
	fmt.Println(r)
	r = m.RealQuickRatio()
	fmt.Println(r)
	a = "spam"
	b = "park"
	m.SetSeqs(a, b)
	e := m.GetEditops()
	fmt.Println(FormatLevEditOps(e))
	ep := m.GetOpcodes()
	fmt.Println(FormatLevOpCodes(ep))
	ms := m.GetMatchingblocks()
	fmt.Println(FormatLevMatchingBlocks(ms))
}
