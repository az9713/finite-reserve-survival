"""Run the preregistered finite-fuel organization comparison and draw its figures."""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from purposeful.organization import Organization, build_system


def draw(data):
    bg,ink,muted='#0d141d','#e6edf7','#adc0d6'
    colors={'spectator':'#d899e8','fixed':'#72dbc1','passive':'#edbd75'}
    labels={'spectator':'Fuel-spectator organization','fixed':'Best sampled fixed intensity','passive':'Passive organization only'}
    plt.rcParams.update({'figure.facecolor':bg,'axes.facecolor':bg,'text.color':ink,
        'axes.labelcolor':ink,'axes.edgecolor':muted,'xtick.color':muted,'ytick.color':muted,
        'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    times=np.array(data['times']); fine=data['refinements'][-1]
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for name,curve in fine['survival'].items():
        axes[0].plot(times,curve,label=labels[name],color=colors[name],linestyle='--' if name=='fixed' else '-')
    axes[0].set(xlabel='Time (model units)',ylabel='Probability of no first collapse',title='Identical dispersed start and initial fuel')
    axes[0].legend(facecolor=bg,labelcolor=ink,fontsize=8)
    for row in data['refinements']:
        delta=np.array(row['survival']['spectator'])-row['survival']['fixed']
        axes[1].plot(times,100*delta,label=f'n = {row["n"]}')
    axes[1].axhline(0,color=muted,linewidth=.8)
    axes[1].set(xlabel='Time (model units)',ylabel='Spectator minus fixed (percentage points)',title='Small policy gap across spatial refinements')
    axes[1].legend(facecolor=bg,labelcolor=ink,fontsize=8)
    save(fig,'organization-survival')
    obs=data['observables']; t=np.array(data['observation_times'])
    fig,axes=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
    for name,values in obs.items():
        mass=np.asarray(values['mass'])
        axes[0].plot(t,np.asarray(values['tether'])/mass,label=labels[name],color=colors[name])
        axes[1].plot(t,np.asarray(values['gate'])/mass,color=colors[name])
        axes[2].plot(t,values['stopped_fuel'],color=colors[name])
    axes[0].set(xlabel='Time',ylabel='Tether closed, conditional on survival',title='Organization emerges from q = 0')
    axes[1].set(xlabel='Time',ylabel='Mean contact activity, conditional on survival',title='More contact need not mean more survival')
    axes[2].set(xlabel='Time',ylabel='Expected remaining fuel at first collapse or t',title='Repair and organization share one inventory')
    axes[0].legend(facecolor=bg,labelcolor=ink,fontsize=7)
    save(fig,'organization-resources')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    curve=data['fixed_sweep']
    axes[0].plot(curve['alpha'],100*np.array(curve['survival']),color=colors['fixed'],marker='.',label='Fixed family at T = 2')
    axes[0].axhline(100*curve['spectator'],color=colors['spectator'],linestyle='--',label='Spectator candidate, same grid')
    axes[0].set(xlabel='Constant paid-channel prefactor alpha',ylabel='Survival (%)',title=f'Full fixed-intensity interval (plot: n={curve["n"]})')
    axes[0].legend(facecolor=bg,labelcolor=ink,fontsize=8)
    certificates=data['certificate_summary']
    ns=[c['n'] for c in certificates]
    low=np.array([c['gap'][0] for c in certificates])*100
    high=np.array([c['gap'][1] for c in certificates])*100
    axes[1].errorbar(ns,(low+high)/2,yerr=(high-low)/2,fmt='o',capsize=4,color=colors['spectator'])
    axes[1].axhline(0,color=muted,linewidth=.8)
    axes[1].set(xlabel='Full lattice points per axis n',ylabel='Certified candidate minus best fixed (percentage points)',title='Finite-grid certificate; no continuum enclosure')
    save(fig,'organization-fixed-certificate')


def save(fig,name):
    for extension in ['png','pdf']:
        fig.savefig(ROOT/f'results/figures/{name}.{extension}',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grids',nargs='+',type=int,default=[16,24,32,48,64])
    parser.add_argument('--observables-grid',type=int,default=32)
    args=parser.parse_args()
    verification_path=ROOT/'results/organization-verification.json'
    verification_bytes=verification_path.read_bytes()
    verification=json.loads(verification_bytes)
    model=Organization(); times=np.linspace(0,2,201)
    fixed_alpha=verification['certificates'][-1]['fixed']['attaining_sample_alpha']
    policies={'spectator':'spectator','fixed':fixed_alpha,'passive':0.}
    refinements=[]
    for n in args.grids:
        started=time.perf_counter();system=build_system(model,n)
        curves={};diagnostics={}
        for name,policy in policies.items():
            record={};values=system.survival(policy,times,diagnostics=record)
            assert values.min()>=-1e-10 and values.max()<=1+1e-10
            assert np.diff(values).max()<1e-10
            curves[name]=values.tolist();diagnostics[name]=record
        refinements.append({'n':n,'states':len(system.initial()),'survival':curves,
            'diagnostics':diagnostics,'seconds':time.perf_counter()-started})
        print(json.dumps({'n':n,'endpoints':{k:v[-1] for k,v in curves.items()},
                          'seconds':refinements[-1]['seconds']}),flush=True)
    obs_times=np.linspace(0,2,21);system=build_system(model,args.observables_grid)
    observables={};ledgers={}
    for name,policy in policies.items():
        values=system.observables(policy,obs_times)
        fuel=values['stopped_fuel']+values['paid_forward']-values['paid_reverse']+values['repair_forward']-values['repair_reverse']-model.fuel
        tether=values['stopped_tether']-values['paid_forward']+values['paid_reverse']-values['passive_forward']+values['passive_reverse']
        ledgers[name]={'maximum_stopped_fuel_error':float(np.abs(fuel).max()),
                       'maximum_stopped_tether_error':float(np.abs(tether).max())}
        assert max(ledgers[name].values())<1e-9
        observables[name]={k:v.tolist() for k,v in values.items()}
    alphas=np.linspace(0,1,21)
    sweep=[float(system.survival(float(a),[2])[0]) for a in alphas]
    certificate_summary=[{'n':c['n'],'fixed_lower':c['fixed']['lower'],'fixed_upper':c['fixed']['upper'],
        'attaining_sample_alpha':c['fixed']['attaining_sample_alpha'],
        'spectator_interval':c['spectator_interval'],'gap':c['spectator_minus_best_fixed_interval']}
        for c in verification['certificates']]
    inputs=['src/purposeful/organization.py','src/purposeful/spatial_repair.py','scripts/run_organization.py',
            'research/organization-preregistration.md','research/organization-originality-audit.md',
            'research/organization-search-receipt.json','research/organization-model-design.md']
    record={'generated_utc':datetime.now(timezone.utc).isoformat(),'model':asdict(model),
        'horizon':2,'initial':{'h':1,'fuel':model.fuel,'q':0,'position':'uniform periodic relative position'},
        'fixed_scope':'All constant paid-channel prefactors alpha in [0,1], common other parameters; not all fixed physical architectures.',
        'policy_used_for_fixed_curves':fixed_alpha,'times':times.tolist(),'refinements':refinements,
        'observation_times':obs_times.tolist(),'observables_grid':args.observables_grid,'observables':observables,
        'stopped_ledger_checks':ledgers,'fixed_sweep':{'n':args.observables_grid,'alpha':alphas.tolist(),
        'survival':sweep,'spectator':float(system.survival('spectator',[2])[0])},
        'certificate_summary':certificate_summary,'verification_sha256':hashlib.sha256(verification_bytes).hexdigest(),
        'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'continuum_error_enclosed':False}
    (ROOT/'results/organization.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    draw(record)


if __name__=='__main__':main()
