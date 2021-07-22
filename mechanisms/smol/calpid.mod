TITLE CaL pid
COMMENT
    CaL channel with builtin PID controller (does that work)
ENDCOMMENT
NEURON {
    SUFFIX calpid
    USEION ca WRITE ica, cao VALENCE 2 ? Assuming valence = 2 (Ca ion); TODO check this!!
    RANGE gion
    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!
    RANGE conductance                       : parameter
    RANGE g                                 : exposure
    RANGE fopen                             : exposure
    RANGE k_instances                       : parameter
    RANGE k_tau                             : exposure
    RANGE k_inf                             : exposure
    RANGE k_rateScale                       : exposure
    RANGE k_fcond                           : exposure
    RANGE k_steadyState_rate                : parameter
    RANGE k_steadyState_midpoint            : parameter
    RANGE k_steadyState_scale               : parameter
    RANGE k_steadyState_x                   : exposure
    RANGE k_timeCourse_tau                  : parameter
    RANGE k_timeCourse_t                    : exposure
    RANGE l_instances                       : parameter
    RANGE l_tau                             : exposure
    RANGE l_inf                             : exposure
    RANGE l_rateScale                       : exposure
    RANGE l_fcond                           : exposure
    RANGE l_steadyState_rate                : parameter
    RANGE l_steadyState_midpoint            : parameter
    RANGE l_steadyState_scale               : parameter
    RANGE l_steadyState_x                   : exposure
    RANGE l_timeCourse_TIME_SCALE           : parameter
    RANGE l_timeCourse_VOLT_SCALE           : parameter
    RANGE l_timeCourse_t                    : exposure
    RANGE k_tauUnscaled                     : derived variable
    RANGE l_timeCourse_V                    : derived variable
    RANGE l_tauUnscaled                     : derived variable
    RANGE conductanceScale                  : derived variable
    RANGE fopen0                            : derived variable
}
UNITS {
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
}
PARAMETER {
    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!
    conductance = 1.0E-5 (uS)
    k_instances = 3 
    k_steadyState_rate = 1 
    k_steadyState_midpoint = -61 (mV)
    k_steadyState_scale = 4.2 (mV)
    k_timeCourse_tau = 1 (ms)
    l_instances = 1 
    l_steadyState_rate = 1 
    l_steadyState_midpoint = -85.5 (mV)
    l_steadyState_scale = -8.5 (mV)
    l_timeCourse_TIME_SCALE = 1 (ms)
    l_timeCourse_VOLT_SCALE = 1 (mV)
}
ASSIGNED {
    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel
    v (mV)
    celsius (degC)
    temperature (K)
    eca (mV)
    k_steadyState_x                        : derived variable
    k_timeCourse_t (ms)                    : derived variable
    k_rateScale                            : derived variable
    k_fcond                                : derived variable
    k_inf                                  : derived variable
    k_tauUnscaled (ms)                     : derived variable
    k_tau (ms)                             : derived variable
    l_steadyState_x                        : derived variable
    l_timeCourse_V                         : derived variable
    l_timeCourse_t (ms)                    : derived variable
    l_rateScale                            : derived variable
    l_fcond                                : derived variable
    l_inf                                  : derived variable
    l_tauUnscaled (ms)                     : derived variable
    l_tau (ms)                             : derived variable
    conductanceScale                       : derived variable
    fopen0                                 : derived variable
    fopen                                  : derived variable
    g (uS)                                 : derived variable
    rate_k_q (/ms)
    rate_l_q (/ms)
    rate_vstd (/ms)
    rate_gmax (/ms)
}
STATE {
    k_q
    l_q
    gmax_actual
    vmean
    vstd
    timer
}
INITIAL {
    eca = 120.0
    temperature = celsius + 273.15
    rates(v, )
    k_q = k_inf
    l_q = l_inf
    gmax_actual = gmax
    vmean = v
    vstd = 0.4
    timer = 3
}
BREAKPOINT {
    SOLVE states METHOD cnexp
    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=cal type=ionChannelHH), from conductanceScaling; null
    ? Path not present in component, using factor: 1
    conductanceScale = 1 
    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=cal type=ionChannelHH), from gates; Component(id=k type=gateHHtauInf)
    ? multiply applied to all instances of fcond in: <gates> ([Component(id=k type=gateHHtauInf), Component(id=l type=gateHHtauInf)]))
    fopen0 = k_fcond * l_fcond ? path based, prefix = 
    fopen = conductanceScale  *  fopen0 ? evaluable
    g = conductance  *  fopen ? evaluable
    gion = gmax_actual * fopen 
    ica = gion * (v - eca)
    if (vstd>0.6) {
        vstd = 0.6
    }
    if (vstd<0.0) {
        vstd = 0
    }
    if (timer < 0) {
        timer = 1
        if (vstd < 0.4) {
            gmax_actual = gmax_actual - 0.003
        }
    }
}
DERIVATIVE states {
    rates(v, )
    k_q' = rate_k_q
    l_q' = rate_l_q
    vmean' = 0.005*(v - vmean)
    vstd' = 0.005*rate_vstd
    ?gmax_actual' = rate_gmax
    timer' = -0.003
}
PROCEDURE rates(v, ) {
    k_steadyState_x = k_steadyState_rate  / (1 + exp(0 - (v -  k_steadyState_midpoint )/ k_steadyState_scale )) ? evaluable
    k_timeCourse_t = k_timeCourse_tau ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=k type=gateHHtauInf), from q10Settings; null
    ? Path not present in component, using factor: 1
    k_rateScale = 1 
    k_fcond = k_q ^ k_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=k type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)
    k_inf = k_steadyState_x ? path based, prefix = k_
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=k type=gateHHtauInf), from timeCourse; Component(id=null type=fixedTimeCourse)
    k_tauUnscaled = k_timeCourse_t ? path based, prefix = k_
    k_tau = k_tauUnscaled  /  k_rateScale ? evaluable
    l_steadyState_x = l_steadyState_rate  / (1 + exp(0 - (v -  l_steadyState_midpoint )/ l_steadyState_scale )) ? evaluable
    l_timeCourse_V = v /  l_timeCourse_VOLT_SCALE ? evaluable
    l_timeCourse_t = l_timeCourse_TIME_SCALE *((20 * exp(( l_timeCourse_V  + 160) / 30) / (1 + exp(( l_timeCourse_V  + 84) / 7.3))) +35) ? evaluable
    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=l type=gateHHtauInf), from q10Settings; null
    ? Path not present in component, using factor: 1
    l_rateScale = 1 
    l_fcond = l_q ^ l_instances ? evaluable
    ? DerivedVariable is based on path: steadyState/x, on: Component(id=l type=gateHHtauInf), from steadyState; Component(id=null type=HHSigmoidVariable)
    l_inf = l_steadyState_x ? path based, prefix = l_
    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=l type=gateHHtauInf), from timeCourse; Component(id=null type=CaL_tau)
    l_tauUnscaled = l_timeCourse_t ? path based, prefix = l_
    l_tau = l_tauUnscaled  /  l_rateScale ? evaluable
    rate_k_q = ( k_inf  -  k_q ) /  k_tau ? Note units of all quantities used here need to be consistent!
    rate_l_q = ( l_inf  -  l_q ) /  l_tau ? Note units of all quantities used here need to be consistent!
    rate_vstd = ((vmean - v)*(vmean-v))^0.5 - vstd
    ?rate_gmax = -0.00001*(vstd-0.4) ? + (gmax-gmax_actual)? take into acount gate state
}
