#!/bin/sh
nevt=1000
pt=-100.0

t0=000
./EXKalTest $nevt $t0 $pt
mv h.root p${pt}t${t0}.root

t0=007
./EXKalTest $nevt $t0 $pt
mv h.root p${pt}t${t0}.root

t0=014
./EXKalTest $nevt $t0 $pt
mv h.root p${pt}t${t0}.root

t0=021
./EXKalTest $nevt $t0 $pt
mv h.root p${pt}t${t0}.root
