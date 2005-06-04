#!/bin/sh
nevt=10000
t0=000

#pt=-1.0
#./EXKalTest $nevt $t0 $pt
#mv h.root p${pt}t${t0}.10k.root

pt=-100.0
./EXKalTest $nevt $t0 $pt
mv h.root p${pt}t${t0}.10k.root

