# Matlab solver of the Hasegawa-Wakatani non linear turbulent system.

![](density.gif)

(Simulation Demo, N=256, L=40, a=0.1, mu=1e-3)

How to use it : 

easiest) go to wk/ and launch run.m with matlab

customable) go to wk/ and check setup.m for the parameters, run it, then run run.m, then analysis.m

more manually) open wk/fort.90, setup your parameters, make from ../, run ../bin/helaz fort.90, then wk/analysis.m

Ref: go check Ammar Hakim's journal http://ammar-hakim.org/sj/je/je17/je17-hasegawa-wakatani.html
