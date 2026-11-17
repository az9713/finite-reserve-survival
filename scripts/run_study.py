"""Generate numerical evidence, standalone figures and offline atlas data."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import sys
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from purposeful.model import (Sensor, sensor_for_accuracy, frontier_power,
    hazard_survival, hazard_decay, hazard_generator, damage_survival,
    damage_generator, simulate_survival, four_state_thermodynamics)


def main(samples):
    out = ROOT/'results'; out.mkdir(exist_ok=True)
    figures = out/'figures'; figures.mkdir(exist_ok=True)
    start = time.monotonic()
    sensors = {'slow': sensor_for_accuracy(.1, .4, .2),
               'fast': sensor_for_accuracy(.1, 5.8, .2)}
    result = {'model': 'binary thermal sensor plus explicitly effective failure layer',
              'units': 'k_B T = 1; arbitrary but shared time unit',
              'horizon':20., 'hazard':.5, 'damage_stages':4, 'damage_duration':2.,
              'calibration':'both laws have mean lethal uninterrupted mismatch duration 2',
              'initialization':'unfailed stationary sensor; damage clock starts fresh',
              'novelty':'candidate synthesis; priority unestablished', 'sensors':{}}
    for idx, (name,s) in enumerate(sensors.items()):
        initial = [1-s.q,s.q]
        damage_initial = [1-s.q,s.q,0,0,0]
        result['sensors'][name] = {
            'r':s.r,'k':s.k,'affinity':s.affinity,'a':s.a,'b':s.b,
            'q':s.q,'d':s.d,'power':s.power,
            'mutual_information_nats':s.mutual_information,
            'learning_rate':s.learning_rate,'efficiency':s.learning_rate/s.power,
            'hazard_survival':hazard_survival(s,20,h=.5),
            'hazard_decay':hazard_decay(s,h=.5),
            'damage_survival':damage_survival(s,20),
            'hazard_mc':simulate_survival(hazard_generator(s,h=.5),initial,20,samples,801+idx),
            'damage_mc':simulate_survival(damage_generator(s),damage_initial,20,samples,901+idx),
        }
    # Independent full-state energy balance across a reproducible random sample.
    rng = np.random.default_rng(20261116)
    checks=[]
    for _ in range(200):
        s=Sensor(10**rng.uniform(-2,1),10**rng.uniform(-2,2),rng.uniform(0,9))
        flux=four_state_thermodynamics(s)
        checks.append(max(abs(flux[k]-s.power) for k in ['heat','work','entropy_production']))
    result['max_flux_residual_200_cases']=max(checks)
    # Robustness: report reversals as counts in this explicit grid, not probabilities.
    transfer=[]
    for q in [.1,.2,.3,.4]:
        minimum_k=.1*(1-2*q)/q
        slow=sensor_for_accuracy(.1,minimum_k*1.4,q)
        fast=sensor_for_accuracy(.1,minimum_k*20,q)
        for horizon in [5.,20.,50.]:
            for stages in [1,2,4,8]:
                for initial_q in [q,0.,1.]:
                    hs=hazard_survival(slow,horizon,h=.5,initial_q=initial_q)
                    hf=hazard_survival(fast,horizon,h=.5,initial_q=initial_q)
                    ds=damage_survival(slow,horizon,stages=stages,initial_q=initial_q)
                    df=damage_survival(fast,horizon,stages=stages,initial_q=initial_q)
                    transfer.append(dict(q=q,horizon=horizon,stages=stages,initial_q=initial_q,
                        slow_k=slow.k,fast_k=fast.k,hazard_slow=hs,hazard_fast=hf,
                        damage_slow=ds,damage_fast=df,
                        reversal=bool(hs>hf+1e-9 and df>ds+1e-9)))
    with (out/'transfer.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(transfer[0]));writer.writeheader();writer.writerows(transfer)
    result['transfer']={'total':len(transfer),'reversals':sum(x['reversal'] for x in transfer),
        'interpretation':'counts over the stated grid; not a probability or universal claim',
        'fresh_damage_clock':True}
    # Atlas: 625 designs at each external switching rate, with several horizons.
    atlas={'r_values':[.03,.1,.3],'horizons':[5.,20.,50.],
           'affinities':np.linspace(.1,6,25).tolist(),
           'capacities':np.geomspace(.05,10,25).tolist(), 'slices':[]}
    for r in atlas['r_values']:
        points=[]
        for k in atlas['capacities']:
            for affinity in atlas['affinities']:
                s=Sensor(r,k,affinity)
                points.append({'k':k,'affinity':affinity,'q':s.q,'power':s.power,
                    'information':s.mutual_information,
                    'hazard':[hazard_survival(s,t,h=.5) for t in atlas['horizons']],
                    'damage':[damage_survival(s,t) for t in atlas['horizons']]})
        atlas['slices'].append({'r':r,'points':points})
    (out/'atlas-data.json').write_text(json.dumps(atlas,separators=(',',':')))
    # Figure 1: same accuracy/information, opposite survival rankings.
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':150,'savefig.bbox':'tight'})
    colors={'slow':'#9258ba','fast':'#007f82'}
    ts=np.linspace(0,50,250)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    for name,s in sensors.items():
        axes[0].plot(ts,[hazard_survival(s,t,h=.5) for t in ts],label=name,color=colors[name])
        axes[1].plot(ts,[damage_survival(s,t) for t in ts],label=name,color=colors[name])
    axes[0].plot(ts,np.exp(-.1*ts),'--',color='#888',label='mean-hazard shortcut')
    for ax,title in zip(axes[:2],['Instantaneous mismatch hazard','Four-stage damage with reset']):
        ax.set(title=title,xlabel='Time',ylabel='Survival probability',ylim=(0,1));ax.legend()
    names=list(sensors)
    axes[2].bar(names,[sensors[n].power for n in names],color=[colors[n] for n in names])
    axes[2].set(title='Stationary sensor work rate',ylabel='Power / (k_B T per time)')
    fig.suptitle('Equal 80% accuracy; survival depends on the damage mechanism',y=1.04)
    for ext in ['png','svg','pdf']:fig.savefig(figures/f'ranking-reversal.{ext}')
    plt.close(fig)
    # Figure 2: bound and kinetic ceiling.
    fig,ax=plt.subplots(figsize=(6.5,4.2));qs=np.linspace(.04,.49,300)
    for cap in [.4,1,5.8,100]:
        x=[];y=[]
        for q in qs:
            try:p=frontier_power(.1,cap,q)
            except ValueError:continue
            x.append(1-q);y.append(p)
        ax.plot(x,y,label=f'K ceiling = {cap:g}')
    ax.set(xlabel='Required stationary accuracy',ylabel='Minimum sensor power / (k_B T per time)',
           title='Exact frontier within the binary sensor family',ylim=(0,1.2))
    ax.legend();fig.savefig(figures/'power-frontier.png');fig.savefig(figures/'power-frontier.svg');plt.close(fig)
    result['elapsed_seconds']=time.monotonic()-start
    (out/'summary.json').write_text(json.dumps(result,indent=2))
    tracked=[p for d in ['src','scripts','tests','docs','research'] for p in (ROOT/d).rglob('*')
             if p.is_file() and p.suffix not in ['.pyc','.pdf'] and 'sources' not in p.parts]
    manifest={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,
              'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
              'command':f'python scripts/run_study.py --samples {samples}',
              'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(tracked)}}
    (out/'run-manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--samples',type=int,default=30000)
    main(parser.parse_args().samples)
