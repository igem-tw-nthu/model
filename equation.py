def lotVolRevised(
    # initial condition
    y,
    # time interval array
    t,
    # parameters 
    ## Vibrio and AHL 
    rV, Vm, a, kA, lmA, thes,
    ## Ecoli and Nisin 
    rE, Em, b, kN, lmN,
    ## mechincal 
    D ):
 
    V , A, E, N = y

    if V<thes: kA=0

    dV_dt = (rV-D)*V - (rV/Vm)*pow(V,2) - a*V*N
    dA_dt = kA*V - lmA*A
    dE_dt = (rE-D)*V - (rE/Em)*pow(E,2) - b*E*A
    dN_dt = kN*E - lmN*N

    dydt = [dV_dt, dA_dt, dE_dt, dN_dt]

    return dydt


def lotVol(
    # initial condition
    y,
    # time interval array
    t,
    # parameters 
    a,b,c,d,lm1,lm2,k1,k2,mt):

    m1,ahl,m2,med = y

    if m1<mt: k1=0
    if ahl<0: ahl=0 
    if med<0: med=0 
    if m2<0: m2=0

    dm1_dt = a*m1 - b*m1*med
    dahl_dt = k1*dm1_dt - lm1*ahl
    dm2_dt = c*m2*ahl - d*m2
    dmed_dt = k2*dm2_dt - lm2*med

    dydt = [ dm1_dt, dahl_dt, dm2_dt, dmed_dt ]

    return dydt