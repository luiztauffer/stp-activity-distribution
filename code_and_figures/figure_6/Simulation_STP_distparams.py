'''
This script simulates the effects on LIF neurons caused by input spike trains
structured or scattered under noisy conditions (synaptic bombardment)

- The target neuron can be set to spike or have a free membrane Vm (large Vth)
- The firing rate and membrane voltage of the neurons are compared
- Presynaptic fibers release are modulated by STP rules
- Local Inhibitory circuit is driven by depressive synapses from the main Poisson
population and have its firing rate linearly modulated by it Irate = a*PRRs2
- Inhibitory spikes to target cells are static

sim_type = 0   -   sparse | dynamic s1 | dynamic s2
sim_type = 1   -   sparse | dynamic s1 | static s2
sim_type = 2   -   sparse | static s1  | dynamic s2
sim_type = 3   -   dense  | dynamic s1 | dynamic s2
'''

from brian2 import *
import numpy as np
import scipy.io as sio


def runsim(filename, sim_type=0):

    nBlocks = 3
    tBlocks = 40
    nTrials = 50  # 3010
    tTime = nTrials * nBlocks * tBlocks   # ms
    simdt = 0.1
    defaultclock.dt = simdt * ms

    # --------------------------------------------------------------------------
    # Biophysical model parameters
    # --------------------------------------------------------------------------

    ref = 2. * ms		      # refractory period
    Cm = 250. * pF        # membrane capacitance

    Vth = -50. * mV             # spike threshold
    Vreset = -55. * mV         # reset value
    Vl =-70. * mV            # resting potential
    Va = 0. * mV          # AMPA inverse potential
    Vg = -75. * mV        # GABA inverse potential
    uVm = -53. * mV       # mean Vm during pure noise

    taua = .5 * ms             # tau AMPA
    taug = 2. * ms             # tau GABA

    gl = 0. * nS        # leak conductance
    Be = 25. * nS          # excitatory max conductance
    Ben = 25.          # excitatory max conductance (for gaussian sampling)
    Bi = 2. * nS          # inhibitory max conductance

    # v_clamp = -55*mV
    # i_clamp = gl*(Vrest-v_clamp)

    # --------------------------------------------------------------------------
    # Equations
    # --------------------------------------------------------------------------

    # Target neuron equations
    eqs_neuron = '''
    dv/dt = (ge*(Va-v) + gi*(Vg-v) + gl*(Vl-v) )/Cm  : volt (unless refractory)

    dge/dt = (ye-ge)/taua : siemens
    dye/dt = -ye/taua  : siemens
    dgi/dt = (yi-gi)/taug : siemens
    dyi/dt = -yi/taug  : siemens
    '''

    # Inhibitory neurons equations
    eqs_inhibitory = '''
        dv/dt = (ye-v)/dt : 1
        dye/dt = -ye/dt : 1
    '''

    # Synaptic models
    model_syn1 = '''
        u : 1
        x : 1
        dpost: siemens
        w : siemens
        tauf: second
        taud: second
        U: 1
        lastupdate: second
    '''

    eqs_syn1 = '''
        u = u*exp(-(t-lastupdate)/tauf)
        u += U*(1-u)
        x = 1 - (1-x)*exp(-(t-lastupdate)/taud)
        dpost = w*u*x
        ye_post += dpost
        x -= u*x
        lastupdate = t
    '''

    model_syn2 = '''
        u : 1
        x : 1
        dpost: 1
        w : 1
        tauf: second
        taud: second
        U: 1
        lastupdate: second
    '''

    eqs_syn2 = '''
        u = u*exp(-(t-lastupdate)/tauf)
        u += U*(1-u)
        x = 1 - (1-x)*exp(-(t-lastupdate)/taud)
        dpost = w*u*x
        ye_post += dpost
        x -= u*x
        lastupdate = t
    '''

    # --------------------------------------------------------------------------
    # Network activity
    # --------------------------------------------------------------------------

    # Mean (u) ans Std (s) of STP parameters
    tauf1u = 0.200 * second
    taud1u = 0.050 * second
    U1u = .1
    tauf1s = tauf1u / 5.
    taud1s = taud1u / 5.
    U1s = U1u / 5.

    tauf2u = .05 * second
    taud2u = 0.2 * second
    U2u = .7
    tauf2s = tauf2u / 5.
    taud2s = taud2u / 5.
    U2s = U2u / 5.

    # Read PRR file and assign values network values ---------------------------
    fio = sio.loadmat(filename)
    Ne = fio['N'][0][0].astype(int)
    Rne = fio['Rn'][0][0]
    # Rne = .5 * Hz
    ExtraRate = fio['Re'][0][0]
    # ExtraRate = 8000. * Hz
    Nopt = fio['Nopt'][0][0]

    # Background activity to keep mean Vm  -------------------------------------
    uinf1 = U1u * (1 + Rne * tauf1u * Hz) / (1 + U1u * Rne * tauf1u * Hz)
    xinf1 = 1 / (1 + uinf1 * Rne * taud1u * Hz)
    AvgRelease1 = uinf1 * xinf1

    De = (Va - uVm) * Be * taua
    Di = (Vg - uVm) * Bi * taug
    DL = (Vl - uVm) * gl

    # Background Inhibition coefficient, Inh_rate = a*PRRs2
    uinf2 = U2u * (1 + Rne * tauf2u * Hz) / (1 + U2u * Rne * tauf2u * Hz)
    xinf2 = 1 / (1 + uinf2 * Rne * taud2u * Hz)
    AvgRelease2 = uinf2 * xinf2

    a = -(De * AvgRelease1) / (Di * AvgRelease2)

    # --------------------------------------------------------------------------
    # Create neurons
    # --------------------------------------------------------------------------

    # Target neuron - spiking
    Target = NeuronGroup(
        1,
        model=eqs_neuron,
        threshold='v > Vth',
        reset='v=Vreset',
        refractory=ref,
        method='euler'
    )
    Target.v = uVm

    # Target neuron - free membrane
    Target_fm = NeuronGroup(
        1,
        model=eqs_neuron,
        threshold='v > -1000*Vth',
        reset='v=Vreset',
        refractory=ref,
        method='euler'
    )
    Target_fm.v = uVm

    # Inhibitory pool
    Ni = 400
    Inhpool = NeuronGroup(
        Ni,
        model=eqs_inhibitory,
        threshold='v>rand()',
        method='euler'
    )

    # --------------------------------------------------------------------------
    # Background activity - connects to both target neurons and local inhibitory pool
    # --------------------------------------------------------------------------

    # Weight dispersion (factor for standard deviation)
    w_disp = 10

    if Ne > 16000:
        Nee = 16000
        Background = PoissonGroup((Ne - Nee), rates=Rne * Hz)

        # Synapse type 1 - facilitatory to target spiking neuron
        SB1 = Synapses(Background, Target, model=model_syn1, on_pre=eqs_syn1)
        SB1.connect()
        SB1.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SB1.w)) / exp(-1))
        SB1.u[:] = uinf1
        SB1.x[:] = xinf1
        SB1.tauf[:] = tauf1u
        SB1.taud[:] = taud1u
        SB1.U[:] = U1u

        # Synapse type 1 - facilitatory to target free membrane neuron
        SB2 = Synapses(Background, Target_fm, model=model_syn1, on_pre=eqs_syn1)
        SB2.connect()
        SB2.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SB2.w)) / exp(-1))
        SB2.u[:] = uinf1
        SB2.x[:] = xinf1
        SB2.tauf[:] = tauf1u
        SB2.taud[:] = taud1u
        SB2.U[:] = U1u

        # Synapse type 2 - depressing to local inhibitory neurons
        SB3 = Synapses(Background, Inhpool, model=model_syn2, on_pre=eqs_syn2)
        spn = int((Ne - Nee) / Ni)   # Synapses per neuron
        for ki in range(0, Ni):
            find = ki*spn
            SB3.connect(i=range(find, find + spn), j=ki)
        SB3.w[:] = np.abs(np.random.normal(a, a/w_disp, len(SB3.w)))
        SB3.u[:] = uinf2
        SB3.x[:] = xinf2
        SB3.tauf[:] = tauf2u
        SB3.taud[:] = taud2u
        SB3.U[:] = U2u
    else:
        Nee = np.copy(Ne)

    # --------------------------------------------------------------------------
    # Extra activity - connects to both target neurons and local inhibitory pool
    # --------------------------------------------------------------------------

    # Sparse presynaptic activity
    if sim_type != 3:
        rate_matrix = np.ones([int(tTime / tBlocks), Nee]) * Rne
        burst_blocks = np.arange(nBlocks / 2, tTime / tBlocks, nBlocks).astype(int)
        for j in burst_blocks:
            # selected bursty units
            chosen_units = np.random.choice(Nee, size=Nopt, replace=False)
            rate_matrix[j, chosen_units] += ExtraRate / Nopt
        StimE = TimedArray(rate_matrix, dt=tBlocks * ms)
        PE = PoissonGroup(Nee, rates='StimE(t,i)*Hz')
    # Dense presynaptic activity
    else:
        rate_matrix = np.ones(int(tTime / tBlocks)) * Rne
        burst_blocks = np.arange(nBlocks / 2, tTime / tBlocks, nBlocks).astype(int)
        for j in burst_blocks:
            # selected bursty units
            chosen_units = np.random.choice(Nee, size=Nopt, replace=False)
            rate_matrix[j] += float(ExtraRate) / Nee
        StimE = TimedArray(rate_matrix, dt=tBlocks * ms)
        PE = PoissonGroup(Nee, rates='StimE(t)*Hz')

    # Synapse type 1 - facilitatory to target spiking neuron
    tauf1array = np.abs(np.random.normal(tauf1u, tauf1s, Nee)) * second
    taud1array = np.abs(np.random.normal(taud1u, taud1s, Nee)) * second
    U1array = np.abs(np.random.normal(U1u, U1s, Nee))
    ui1array = U1array * (1 + Rne * tauf1array * Hz) / (1 + U1array * Rne * tauf1array * Hz)
    xi1array = 1 / (1 + ui1array * Rne * taud1array * Hz)

    # Synapse type 2 - depressing to to local inhibitory neurons
    tauf2array = np.abs(np.random.normal(tauf2u, tauf2s, Nee)) * second
    taud2array = np.abs(np.random.normal(taud2u, taud2s, Nee)) * second
    U2array = np.abs(np.random.normal(U2u, U2s, Nee))
    ui2array = U2array * (1 + Rne * tauf2array * Hz) / (1 + U2array * Rne * tauf2array * Hz)
    xi2array = 1 / (1 + ui2array * Rne * taud2array * Hz)

    # sparse | dynamic s1 | dynamic s2
    if sim_type == 0:
        # S1 - signal population to target spiking neuron
        SE = Synapses(PE, Target, model=model_syn1, on_pre=eqs_syn1)
        SE.connect()
        SE.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SE.w)) / exp(-1))
        SE.u[:] = ui1array
        SE.x[:] = xi1array
        SE.tauf[:] = tauf1array
        SE.taud[:] = taud1array
        SE.U[:] = U1array

        # S1 - signal population to target free membrane neuron
        SEfm = Synapses(PE, Target_fm, model=model_syn1, on_pre=eqs_syn1)
        SEfm.connect()
        SEfm.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SEfm.w)) / exp(-1))
        SEfm.u[:] = ui1array
        SEfm.x[:] = xi1array
        SEfm.tauf[:] = tauf1array
        SEfm.taud[:] = taud1array
        SEfm.U[:] = U1array

        # S2 - signal population to local inhibitory neurons
        SEI = Synapses(PE, Inhpool, model=model_syn2, on_pre=eqs_syn2)
        spn = int(Nee / Ni)   # Synapses per neuron
        for ki in range(0, Ni):
            find = ki * spn
            SEI.connect(i=range(find, find + spn), j=ki)
        SEI.w[:] = np.abs(np.random.normal(a, a/w_disp, len(SEI.w)))
        SEI.u[:] = ui2array
        SEI.x[:] = xi2array
        SEI.tauf[:] = tauf2array
        SEI.taud[:] = taud2array
        SEI.U[:] = U2array

    # sparse | dynamic s1 | static s2
    elif sim_type == 1:
        # S1 - signal population to target spiking neuron
        SE = Synapses(PE, Target, model=model_syn1, on_pre=eqs_syn1)
        SE.connect()
        SE.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SE.w)) / exp(-1))
        SE.u[:] = ui1array
        SE.x[:] = xi1array
        SE.tauf[:] = tauf1array
        SE.taud[:] = taud1array
        SE.U[:] = U1array

        # S1 - signal population to target free membrane neuron
        SEfm = Synapses(PE, Target_fm, model=model_syn1, on_pre=eqs_syn1)
        SEfm.connect()
        SEfm.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SEfm.w)) / exp(-1))
        SEfm.u[:] = ui1array
        SEfm.x[:] = xi1array
        SEfm.tauf[:] = tauf1array
        SEfm.taud[:] = taud1array
        SEfm.U[:] = U1array

        # S2 - signal population to local inhibitory neurons
        SEI = Synapses(PE, Inhpool, model='avgrel : 1', on_pre='ye_post+=a*avgrel')
        spn = int(Nee / Ni)   # Synapses per neuron
        for ki in range(0, Ni):
            find = ki * spn
            SEI.connect(i=range(find, find + spn), j=ki)
        SEI.avgrel[:] = AvgRelease2

    # sparse | static s1 | dynamic s2
    elif sim_type == 2:
        # S1 - signal population to target spiking neuron
        SE = Synapses(PE, Target, model='avgrel : 1', on_pre='ye_post+=(Be*avgrel)/exp(-1)')
        SE.connect()
        SE.avgrel[:] = AvgRelease1

        # S1 - signal population to target free membrane neuron
        SEfm = Synapses(PE, Target_fm, model='avgrel : 1', on_pre='ye_post+=(Be*avgrel)/exp(-1)')
        SEfm.connect()
        SEfm.avgrel[:] = AvgRelease1

        # S2 - signal population to local inhibitory neurons
        SEI = Synapses(PE, Inhpool, model=model_syn2, on_pre=eqs_syn2)
        spn = int(Nee / Ni)   # Synapses per neuron
        for ki in range(0, Ni):
            find = ki * spn
            SEI.connect(i=range(find, find + spn), j=ki)
        SEI.w[:] = np.abs(np.random.normal(a, a/w_disp, len(SEI.w)))
        SEI.u[:] = ui2array
        SEI.x[:] = xi2array
        SEI.tauf[:] = tauf2array
        SEI.taud[:] = taud2array
        SEI.U[:] = U2array

    # dense | dynamic s1 | dynamic s2
    elif sim_type == 3:
        # S1 - signal population to target spiking neuron
        SE = Synapses(PE, Target, model=model_syn1, on_pre=eqs_syn1)
        SE.connect()
        SE.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SE.w)) / exp(-1))
        SE.u[:] = ui1array
        SE.x[:] = xi1array
        SE.tauf[:] = tauf1array
        SE.taud[:] = taud1array
        SE.U[:] = U1array

        # S1 - signal population to target free membrane neuron
        SEfm = Synapses(PE, Target_fm, model=model_syn1, on_pre=eqs_syn1)
        SEfm.connect()
        SEfm.w[:] = nS * np.abs(np.random.normal(Ben, Ben/w_disp, len(SEfm.w)) / exp(-1))
        SEfm.u[:] = ui1array
        SEfm.x[:] = xi1array
        SEfm.tauf[:] = tauf1array
        SEfm.taud[:] = taud1array
        SEfm.U[:] = U1array

        # S2 - signal population to local inhibitory neurons
        SEI = Synapses(PE, Inhpool, model=model_syn2, on_pre=eqs_syn2)
        spn = int(Nee / Ni)   # Synapses per neuron
        for ki in range(0, Ni):
            find = ki * spn
            SEI.connect(i=range(find, find + spn), j=ki)
        SEI.w[:] = np.abs(np.random.normal(a, a/w_disp, len(SEI.w)))
        SEI.u[:] = ui2array
        SEI.x[:] = xi2array
        SEI.tauf[:] = tauf2array
        SEI.taud[:] = taud2array
        SEI.U[:] = U2array

    del rate_matrix

    # --------------------------------------------------------------------------
    # Inhibitory activity - connects to both target neurons
    # --------------------------------------------------------------------------

    SI = Synapses(Inhpool, Target, on_pre='yi_post+=Bi/exp(-1)')
    SI.connect()

    SIfm = Synapses(Inhpool, Target_fm, on_pre='yi_post+=Bi/exp(-1)')
    SIfm.connect()

    # --------------------------------------------------------------------------
    # Set up monitors and run simulation ---------------------------------------
    # --------------------------------------------------------------------------

    # Create monitors
    stateTarget = StateMonitor(Target_fm, 'v', record=True)
    GiTarget = StateMonitor(Target_fm, 'gi', record=True)
    GeTarget = StateMonitor(Target_fm, 'ge', record=True)
    spikeTarget = SpikeMonitor(Target)
    # spikeInhpool = SpikeMonitor(Inhpool)
    # stateSE = StateMonitor(SE, ['dpost'], record=True)
    # stateSI = StateMonitor(SI, ['dpost'], record=range(100))

    # Run simulation
    run(tTime * ms)

    # Output variables
    spkt = spikeTarget.t[np.where(spikeTarget.i == 0)[0]] / ms
    # spkt_Inh = spikeInhpool.t[np.where(spikeInhpool.i==0)[0]]/ms
    VmP = stateTarget.v[0]
    GI = GiTarget.gi[0]
    GE = GeTarget.ge[0]

    return spkt, VmP, GI, GE
