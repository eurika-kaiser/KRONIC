function sysr = getBRSys(sys,B)
sys.B = B;
[sysr,g,T,Ti] = balreal(sys);







