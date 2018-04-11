function T = getBRTrafo(sys,B)
sys.B = B;
[sysr,g,T,Ti] = balreal(sys);







