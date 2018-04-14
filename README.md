## levenshtein

It's a Levenshtein Distance Golang package with the similar and complete API with Python-Levensthtein.

## Features

A Golang package for fast computation of:
 - Levenshtein (edit) distance and edit sequence manipulation
 - string similarity
 - approximate median strings, and generally string averaging
 - string sequence and set similarity

## Installation and Testing

### Install the Go package:
```
go get -u -v github.com/conetse/levenshtein
```
or, you can install it as follow:
```
cd $GOPATH/src/github.com/conetse/
git clone https://github.com/conetse/levenshtein.git
```

### Run the Tests
```
cd $GOPATH/src/github.com/conetse/levenshtein
go test
```

## Example and Usage

```
package main

import (
    "fmt"
    lev "github.com/conetse/levenshtein"
)

func test_levenshtein() {
    var s1, s2 string
    var hdist int
    var ratio, dist float64
    s1 = "Levenshtein"
    s2 = "Lenvinsten"
    //s1 = "你好一二三四五"
    //s2 = "世界abcde"
    hdist = lev.Distance(s1, s2)
    // 4
    fmt.Println(hdist)
    s1 = "Hello world!"
    s2 = "Holly grail!"
    ratio = lev.Ratio(s1, s2)
    // 0.5833333333333334
    fmt.Println(ratio)
    //
    hdist = lev.Hamming(s1, s2)
    // 7
    fmt.Println(hdist)
    s1 = "Thorkel"
    s2 = "Thorgier"
    dist = lev.Jaro(s1, s2)
    // 0.7797619047619048
    fmt.Println(dist)
    dist = lev.Jaro_winkler(s1, s2)
    // 0.8678571428571429
    fmt.Println(dist)
}

func main() {
    test_levenshtein()
}
```


## Links

The similar API Python module is here [Python-Levensthtein](https://pypi.python.org/pypi/python-Levenshtein).

## License
This project is under the MIT License.

