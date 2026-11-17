"""Independent event-driven evidence, counterexamples, and publication figures."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import json
import hashlib
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from purposeful.assembly import Assembly, success_interval


def simulate(model, alpha, z, samples=30000, seed=1731):
    """Independent six-channel Gillespie implementation with a finite inventory."""
    rng = np.random.default_rng(seed)
    a = float(model.assembly_rate)*float(alpha)
    r = float(model.repair_rate)*(1-float(alpha))*float(z)
    hrate = float(model.damage_rate)
    cache = {}
    def outgoing(state):
        if state in cache:
            return cache[state]
        h, d, b = state
        free, waste = model.pairs-h-d, model.fuel-b
        changes = [(1,0,-1), (-1,0,1), (-1,1,0), (1,-1,0), (1,-1,-1), (-1,1,1)]
        rates = [a*free**2*b, a*float(z)/model.fuel_factor*h*waste,
                 hrate*h, hrate/model.damage_factor*d, r*h*d*b,
                 r*model.damage_factor/model.fuel_factor*h*(h-1)*waste]
        allowed = [(j, tuple(s+c for s,c in zip(state,change)), rate)
                   for j,(change,rate) in enumerate(zip(changes,rates)) if rate > 0]
        cache[state] = allowed
        return allowed
    deadline, lifetime = 1., 4.
    formed, successful = 0, 0
    assembly_times, failures, successful_repairs = [], [], 0
    example = None
    max_energy_residual = 0.
    for _ in range(samples):
        state, now, first = (0,0,model.fuel), 0., None
        net_assembly = net_repair = repair_events = 0
        heat = 0.
        trace = [[0., *state]]
        while True:
            channels = outgoing(state)
            total = sum(edge[2] for edge in channels)
            if not total:
                break
            delay = rng.exponential(1/total)
            horizon = deadline if first is None else first+lifetime
            if now+delay > horizon:
                if first is not None:
                    successful += 1
                    successful_repairs += repair_events
                    if repair_events and example is None:
                        example = trace + [[horizon-first, *state]]
                break
            now += delay
            pick = rng.random()*total
            cumulative = 0.
            for channel, target, rate in channels:
                cumulative += rate
                if pick <= cumulative:
                    break
            old = state
            state = target
            if channel in [0, 1]:
                net_assembly += 1 if channel == 0 else -1
            if channel in [4, 5]:
                net_repair += 1 if channel == 4 else -1
                repair_events += channel == 4
            h,d,b = state
            assert net_assembly == h+d
            assert net_assembly+net_repair == model.fuel-b
            def energy(s):
                h,d,b = s
                return np.log(float(z))*(h+d)-np.log(model.damage_factor)*d+np.log(model.fuel_factor)*b
            heat += energy(old)-energy(state)
            max_energy_residual = max(max_energy_residual, abs(heat-(energy((0,0,model.fuel))-energy(state))))
            if first is None:
                first = now
                formed += 1
                assembly_times.append(first)
                trace = [[0., *state]]
            else:
                trace.append([now-first, *state])
            if h == 0:
                failures.append(now-first)
                break
    probability = successful/samples
    zz = 1.959963984540054
    center = (probability+zz**2/(2*samples))/(1+zz**2/samples)
    radius = zz*np.sqrt(probability*(1-probability)/samples+zz**2/(4*samples**2))/(1+zz**2/samples)
    return {'seed': seed, 'samples': samples, 'formed': formed, 'successful': successful,
            'estimate': probability, 'wilson_95': [center-radius,center+radius],
            'forward_repairs_on_successful_paths': successful_repairs,
            'max_path_energy_residual': max_energy_residual,
            'assembly_times': assembly_times, 'failure_times': failures, 'example': example}


def main():
    model, alpha, z = Assembly(fuel=4), F(33,64), F(1,4)
    simulation = simulate(model,alpha,z)
    exact = model.formation_and_maintenance(1,4,alpha,z)
    se = np.sqrt(exact*(1-exact)/simulation['samples'])
    assert abs(simulation['estimate']-exact) < 5*se
    assert simulation['max_path_energy_residual'] < 1e-12
    assert simulation['forward_repairs_on_successful_paths'] > 0
    bg, ink, muted, teal, amber = '#0d141d','#e6edf7','#adc0d6','#72dbc1','#edbd75'
    plt.rcParams.update({'figure.facecolor':bg,'axes.facecolor':bg,'text.color':ink,
        'axes.labelcolor':ink,'axes.edgecolor':muted,'xtick.color':muted,'ytick.color':muted,
        'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    directory=ROOT/'results/figures'
    def save(fig,name):
        fig.savefig(directory/(name+'.png'),dpi=180)
        fig.savefig(directory/(name+'.pdf'))
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    axes[0].axis('off')
    labels=[('P + Q + fuel','H + waste',.8),('Active complex H','Damaged complex D',.5),('H + D + fuel','2 H + waste',.2)]
    for left,right,y in labels:
        axes[0].text(.02,y,left,va='center',color=teal)
        axes[0].text(.98,y,right,va='center',ha='right',color=amber)
        axes[0].annotate('',xy=(.59,y),xytext=(.41,y),arrowprops={'arrowstyle':'<->','color':ink})
    axes[0].text(.5,.69,'Assembly consumes fuel; disassembly returns it',ha='center',fontsize=9)
    axes[0].text(.5,.39,'Damage and thermal restoration exchange heat',ha='center',fontsize=9)
    axes[0].set_title('Assembly creates the catalyst needed for repair')
    axes[1].axis('off')
    axes[1].text(.05,.82,'Finite inventory',fontsize=18,color=teal)
    axes[1].text(.05,.60,'P + H + D = Q + H + D = 2\nFuel + waste = 4',linespacing=1.8)
    axes[1].text(.05,.32,r'$E/(k_BT)=-\theta(h+d)-\delta d+\epsilon b$',fontsize=12)
    axes[1].text(.05,.10,'Every enabled channel retains its physical reverse.\nThe experiment stops at the first loss of all active complexes.',fontsize=10)
    save(fig,'assembly-network')
    grid=np.linspace(0,1,81)
    surface=np.array([[model.formation_and_maintenance(1,4,a,zz) for a in grid] for zz in grid])
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    image=axes[0].imshow(surface,origin='lower',extent=(0,1,0,1),aspect='auto',cmap='viridis')
    axes[0].scatter([float(alpha)],[float(z)],color=amber,marker='x',s=90)
    axes[0].set(xlabel='Assembly capacity fraction alpha',ylabel='Inverse binding z = exp(-theta)',title='Assemble by time 1, then operate for 4')
    fig.colorbar(image,ax=axes[0],label='Success probability')
    for a in [.25,.5,.75,1]:
        axes[1].plot(grid,[model.formation_and_maintenance(1,4,a,zz) for zz in grid],label=f'alpha = {a:g}')
    axes[1].set(xlabel='Inverse binding z',ylabel='Success probability',title='Intermediate binding helps in this kinetic family')
    axes[1].legend(facecolor=bg,labelcolor=ink)
    save(fig,'assembly-frontier')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    trace=np.array(simulation['example'])
    for column,label,color in [(1,'Active complexes',teal),(2,'Damaged complexes',amber),(3,'Fuel',muted)]:
        axes[0].step(trace[:,0],trace[:,column],where='post',label=label,color=color)
    axes[0].set(xlabel='Time since first assembly',ylabel='Molecule count',title='A seeded successful path with paid repair')
    axes[0].legend(facecolor=bg,labelcolor=ink)
    times=np.linspace(0,4,81)
    empirical=np.array([sum(t>u for t in simulation['failure_times'])+simulation['successful'] for u in times])/simulation['formed']
    # At u=4, failure paths have already failed and successful paths remain.
    exact_survival=[model.formation_and_maintenance(1,t,alpha,z)/(-np.expm1(-float(alpha)*16)) for t in times]
    axes[1].plot(times,exact_survival,color=teal,label='Killed generator')
    axes[1].plot(times,empirical,'--',color=amber,label='Independent Gillespie paths')
    axes[1].set(xlabel='Time since first assembly',ylabel='Conditional persistence',title='30,000 seeded assembly attempts')
    axes[1].legend(facecolor=bg,labelcolor=ink)
    save(fig,'assembly-trajectories')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    theta=np.linspace(0,6,101)
    for pairs in [1,2]:
        current=Assembly(pairs=pairs,fuel=4)
        axes[0].plot(theta,[current.formation_and_maintenance(1,4,alpha,np.exp(-t)) for t in theta],label=f'{pairs} monomer pair(s)')
    axes[0].set(xlabel='Binding strength theta / (kBT)',ylabel='Success probability',title='One complex has no cooperative repair')
    axes[0].legend(facecolor=bg,labelcolor=ink)
    fuels=list(range(9))
    exact_fuel=[Assembly(fuel=n).formation_and_maintenance(1,4,alpha,z) for n in fuels]
    axes[1].plot(fuels,exact_fuel,'o-',color=teal)
    axes[1].errorbar([4],[simulation['estimate']],yerr=[[simulation['estimate']-simulation['wilson_95'][0]],[simulation['wilson_95'][1]-simulation['estimate']]],fmt='s',color=amber,capsize=4,label='Monte Carlo 95% interval')
    axes[1].set(xlabel='Initial fuel molecules',ylabel='Success probability',title='Fixed kinetic design across finite reserves')
    axes[1].legend(facecolor=bg,labelcolor=ink)
    save(fig,'assembly-resources')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    formation_times=np.sort(simulation['assembly_times'])
    axes[0].plot(formation_times,np.arange(1,len(formation_times)+1)/simulation['samples'],color=amber,label='Empirical formation distribution')
    ft=np.linspace(0,1,201)
    axes[0].plot(ft,-np.expm1(-float(alpha)*16*ft),'--',color=teal,label='Exponential law')
    axes[0].set(xlabel='Time from unassembled launch',ylabel='Probability formed',title='The first active complex has an exact arrival law')
    axes[0].legend(facecolor=bg,labelcolor=ink,fontsize=9)
    af=np.linspace(.005,1,121)
    axes[1].plot(af,[model.joint(4,a,z) for a in af],color=amber,label='Assembled and alive at absolute time 4')
    axes[1].plot(af,[model.formation_and_maintenance(1,4,a,z) for a in af],color=teal,label='Form by 1 and persist for 4')
    axes[1].set(xlabel='Assembly capacity fraction alpha',ylabel='Event probability',title='A common deadline can reward late formation')
    axes[1].legend(facecolor=bg,labelcolor=ink,fontsize=9)
    save(fig,'assembly-formation')
    simulation.pop('assembly_times')
    simulation.pop('failure_times')
    evidence={'generated_utc':datetime.now(timezone.utc).isoformat(),'design':{'alpha':str(alpha),'z':str(z)},
              'exact_success':exact,'rational_enclosure':list(map(str,success_interval(model,1,4,alpha,z))),
              'simulation':simulation,'fuel_sweep':{'fuel':fuels,'success':exact_fuel},
              'surface':{'alpha':grid.tolist(),'z':grid.tolist(),'success':surface.tolist()},
              'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                               for name in ['src/purposeful/assembly.py','scripts/assembly_evidence.py']}}
    (ROOT/'results/assembly-evidence.json').write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8')
    print('Independent simulation and five figure sets passed.',flush=True)
    print({key:value for key,value in simulation.items() if key!='example'},flush=True)


if __name__=='__main__':
    main()
