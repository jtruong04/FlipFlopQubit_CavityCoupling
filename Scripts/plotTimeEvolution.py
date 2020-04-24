import matplotlib.ticker as ticker
import datetime
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import rc

UBColors = ("#002f56", "#990000", "#6da04b", "#ad841f", "#005bbb", "#e56a54",
            "#ffc72c", "#ebec00", "#00a69c", "#006570", "#2f9fd0", "#002f56")

def plotCavityPopulation(t,rho_cav, Nm = 1, name = ''):
    plt.rc("font", size=24)
    matplotlib.rcParams['font.serif'] = "Times New Roman"
    matplotlib.rcParams['font.family'] = "serif"
    rc('text', usetex=True)

    fig, ax = plt.subplots(figsize=(12, 8))
    for occ in np.arange(rho_cav.shape[1]):
        if(Nm > 1):
            label = rf'$M = {occ/Nm}, N = {occ%Nm}$'
        else:
            label = rf'$N = {occ}$'
        ax.plot(t/1000000, rho_cav[:, occ, occ],
                label=label, color=UBColors[occ%12])
    ax.set_xlabel(r'Time (ms)')
    ax.set_ylabel('')
    ax.set_title('Cavity Occupation')
    plt.grid(alpha=0.25)
    ax.set_ylim([-0.01, 1.01])
    plt.legend(loc=1)
    if (len(name) == 0):
        timeString = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        fileName = 'CavityPopulation_'+timeString+'.png'
    else:
        fileName = name
    plt.savefig('_figures/'+fileName, dpi=600)
    plt.show()
    return fileName


def base10toN(num, base):
    """Change ``num'' to given base
    Upto base 36 is supported."""

    converted_string, modstring = "", ""
    currentnum = num
    if not 1 < base < 37:
        raise ValueError("base must be between 2 and 36")
    if not num:
        return '0'
    while currentnum:
        mod = currentnum % base
        currentnum = currentnum // base
        converted_string = chr(48 + mod + 7*(mod > 10)) + converted_string
    return converted_string

def plotFlipFlopPopulation(t, rho_ff, numDonors, states=[], name = ''):
    plt.rc("font", size=24)
    matplotlib.rcParams['font.serif'] = "Times New Roman"
    matplotlib.rcParams['font.family'] = "serif"
    rc('text', usetex=True)

    fig, ax = plt.subplots(figsize=(12, 8))
    if len(states) == 0:
        for state in np.arange(rho_ff.shape[1]):
            statelabel = base10toN(state, 4).rjust(numDonors,'0')
            ax.plot(t/1000000, rho_ff[:, state, state], label=rf'$P({statelabel})$', color=UBColors[state % 12])
    else:
        for i, state in enumerate(states):
            statelabel = base10toN(state, 4).rjust(numDonors, '0')
            ax.plot(t/1000000, rho_ff[:, state, state],
                    label=rf'$P({statelabel})$', color=UBColors[i%12])
    ax.set_xlabel(r'Time (ms)')
    ax.set_ylabel('')
    ax.set_title('Flip Flop Population')
    plt.grid(alpha=0.25)
    ax.set_ylim([-0.01, 1.01])
    plt.legend(loc=1)
    if(len(name) == 0) :
        timeString = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        fileName = 'FlipFlopPopulation_'+timeString+'.png'
    else:
        fileName = name
    plt.savefig('_figures/'+fileName, dpi=600)
    plt.show()
    return fileName

def saveParameters(filename, parameters_qubits, parameters_cavity,parameters_noise,detuning,approxInteraction,approxNoise):
    eps = []
    wB = []
    Vt = []
    wc = []
    gc = []
    wn = []
    delta = []
    Nd = len(parameters_qubits)
    Np = len(parameters_cavity)
    Nn = len(parameters_noise)
    # print(parameters_cavity)
    # print(parameters_qubits)
    for donor in range(Nd):
        # print(parameters_qubits[donor])
        eps.append(parameters_qubits[donor]['eps'])
        wB.append(parameters_qubits[donor]['wB'])
        Vt.append(parameters_qubits[donor]['Vt'])
    for mode in range(Np):
        wc.append(parameters_cavity[mode]['wc'])
        gc.append(parameters_cavity[mode]['gc'])
        delta.append(detuning)
    for noise in range(Nn):
        wn.append(parameters_noise[noise]['wn'])
    with open('_figures/Exp_Parameters.csv', 'a+') as f:
        f.write(
            f'{filename},\t"{eps}",\t"{wB}",\t"{Vt}",\t"{wc}",\t"{delta}",\t"{gc}",\t"{wn}",\t"{approxInteraction}",\t"{approxNoise}"\n'
        )
