import tellurium as te
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

mehr_93_simple = """
    # model with a single seeding wave of BM cells
    
    # Species
    var X #proliferating cells
    var Y #resting cells
    var N := X + Y
    
    # Reactions
    -> X; lambda*X*(1-(X + Y)/K);
    X -> ; mu*X;
    X -> Y; sigma*X;
    Y -> ; mu*Y;
    
    #X' = lambda*X*(1-(X + Y)/K) - sigma*X - mu*X
    #Y' = sigma*X - mu*Y
    
    # Species initializations:
    X = 600;
    Y = 0;
    
    # Variable initializations:
    K = 1E5;
    lambda = 1/40;
    sigma = 0.01;   
    mu = 0.001;
    
"""

mehr_95 = """
    # Species
    var N;
    var P;
    var M4;
    var M8;
    var Z := N + P + M4 + M8;
    
    # Reactions
    -> N; (1-N/Kn)*(b + rn*N)
    N ->; dn*N;
    N -> P; sn*N;
    -> P; (1-Z/K)*rp*P;
    P ->; dp*P;
    P -> M4; s4*P;
    P -> M8; s8*P;
    -> M4; (1-Z/K)*r4*M4;
    M4 ->; (d4 + s04)*M4;
    -> M8; (1-Z/K)*r8*M8; 
    M8 ->; (d8 + s08)*M8
    
    
    # Variable initialization
    K = 1E8;
    Kn = 2.E7;
    b = 100;
    rn = 1.1;
    dn = 0.0;
    sn = 1.0;
    rp = 1.0;
    dp = 0.7;
    s4 = 0.04;
    s8 = 0.02;
    r4 = 0.02;
    d4 = 0.0;
    s04 = 0.4;
    r8 = 0.02;
    d8 = 0.0;
    s08 = 0.4;   

"""

mehr_96 = """
    # Species
    var N;
    var P;
    var Ps;
    var M4;
    var M8;
    var S4;
    
    # Formula
    T := N + P + Ps + M4 + M8 + S4;
    K4 := 0.05*T;
    Fcd4 := (S4 + M4)/(S4 + M4 +  K4);
    Frn := 1 + Hrn*Fcd4;
    Fdn := 1 + Hdn*Fcd4;
    Fsn := 1 + Hsn*Fcd4;
    Frp := 1 + Hrp*Fcd4;
    Fdp := 1 + Hdp*Fcd4;
    Fsp := 1 + Hsp*Fcd4;
    Frps := 1 + Hrps*Fcd4;
    Fdps := 1 + Hdps*Fcd4;
    Fs4 := 1 + Hs4*Fcd4;
    Fs8 := 1 + Hs8*Fcd4;
    Fr4 := 1 + Hr4*Fcd4;
    Fd4 := 1 + Hd4*Fcd4;
    Fr8 := 1 + Hr8*Fcd4;
    Fd8 := 1 + Hd8*Fcd4;
    
    
    # Reactions
    -> N; (1-N/Kn)*Frn*rn*N;
    N ->; Fdn*dn*N;
    N -> P; Fsn*sn*N;
    -> P; (1-T/K)*Frp*rp*P;
    P ->; Fdp*dp*P;
    P -> Ps; Fsp*sp*P;
    -> Ps; (1-T/K)*Frps*rps*Ps;
    Ps ->; Fdps*dps*Ps;
    Ps -> M4; Fs4*s4*Ps;
    Ps -> M8; Fs8*s8*Ps;
    -> M4; (1-T/K)*Fr4*r4*M4;
    M4 ->; Fd4*d4*M4;
    -> M8; (1-T/K)*Fr8*r8*M8;
    M8 ->; Fd8*d8*M8;
    S4 ->; ds*S4; 
    
    # Parameter initialization
    K = 2.E5; Kn = 4.E4;
    rn = 1.5; dn = 0.0; sn = 1.0;
    rp = 1.25; dp = 0.6; sp = 0.5;
    rps = 0.0; dps = 0.8;
    r4 = 0.0; d4 = 0.6; s4 = 0.3;
    r8 = 0.0; d8 = 0.6; s8 = 1.0;
    ds = 0.12;
    Hrn = 0.0;
    Hdn = 0.0;
    Hsn = 0.0;
    Hrp = -0.8;
    Hdp = 0.0;
    Hsp = 0.5;
    Hrps = 0.0;
    Hdps = 0.0;
    Hs4 = 1.;
    Hs8 = 0.0;
    Hr4 = 0.0;
    Hd4 = -0.6;
    Hr8 = 0.0;
    Hd8 = 0.0;
    
    # Species initialization
    N = 0.8E5;
    P = 0.0;
    Ps = 0.0;
    M4 = 0.0;
    M8 = 0.0;
    S4 = 0.0;
        
"""

r = te.loada(mehr_96)
res = r.simulate(0, 15, 100)
r.plot()
