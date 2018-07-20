def lotVolRevised(
    t,
    # initial condition
    y,
    # required keyword parameters 
    ## Vibrio and AHL
    *, rV, Vm, a, kA, lmA, 
    ## Ecoli and Nisin and Suicide(GP2)
    b, c, lmE, kN, lmN, kS, lmS, thesV,
    ## mechincal 
    D ):

    """
    Solve the Vibio - Ecoli cometition model.

    This equation this mainly used for scipy.integrate.solve_ivp
    which is a new version of scipy.integrate.odeint 
    """


    V , A, E, N, S = y

    if V>=thesV: kS=0
    else: kN=0

    if V<0: V=0
    if A<0: A=0
    if E<0: E=0
    if N<0: N=0
    if S<0: S=0
    
    dV_dt = (rV-D)*V - (rV/Vm)*pow(V,2) - a*V*N
    dA_dt = kA*V - lmA*A
    dE_dt = b*E*A - c*E*S - lmS*E - D*E
    dN_dt = kN*E - lmN*N
    dS_dt = kS*E - lmS*S


    dydt = [dV_dt, dA_dt, dE_dt, dN_dt, dS_dt]
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