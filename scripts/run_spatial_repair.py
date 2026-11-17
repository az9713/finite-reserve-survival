"""Reproduce the two-complex spatial maintenance experiment and its figures.

The lattice equations are deterministic approximations to a continuous Brownian
model. Independent continuous-space path checks live in verify_spatial_repair.py.
"""
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import asdict
import argparse
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
from purposeful.spatial_repair import (SpatialRepair, build_grid, calibrate,
    count_survival, spatial_survivals)


def hashes():
    paths = [ROOT/'src/purposeful/spatial_repair.py',
             ROOT/'scripts/run_spatial_repair.py',
             ROOT/'research/spatial-repair-originality-audit.md',
             ROOT/'research/spatial-repair-search-receipt.json']
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths}


def figures(data):
    bg, ink, muted = '#0d141d', '#e6edf7', '#adc0d6'
    colors = {'uniform':'#72dbc1', 'localized':'#d899e8', 'reaction':'#edbd75'}
    plt.rcParams.update({'figure.facecolor':bg, 'axes.facecolor':bg,
        'text.color':ink, 'axes.labelcolor':ink, 'axes.edgecolor':muted,
        'xtick.color':muted, 'ytick.color':muted, 'font.size':10,
        'axes.spines.top':False, 'axes.spines.right':False})
    times = np.array(data['times'])
    fine = data['refinements'][-1]
    fig, axes = plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax, initial in zip(axes, ['HD','HH']):
        case=fine['initial_conditions'][initial]
        for name,values in case['spatial'].items():
            label={'reaction':'Gate-weighted position law'}.get(name, name.capitalize()+' start')
            ax.plot(times,values,color=colors[name],label=label)
        ax.plot(times,case['count_edgewise'],'--',color='white',label='Count: each forward mean matched')
        ax.plot(times,case['count_reference'],':',color='#91b5ff',label='Count: reference mean matched')
        ax.set(title=f'{initial} initial chemistry; {data["parameters"]["fuel"]} initial fuel',
               xlabel='Time (model units)',ylabel='Probability of no first collapse',ylim=(0,1))
        ax.legend(facecolor=bg,labelcolor=ink,fontsize=8)
    for extension in ['png','pdf']:
        fig.savefig(ROOT/f'results/figures/spatial-repair-survival.{extension}',dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for refinement in data['refinements']:
        case=refinement['initial_conditions']['HD']
        axes[0].plot(times,np.array(case['spatial']['localized'])-case['count_edgewise'],
                     label=f'n = {refinement["n"]}')
    axes[0].axhline(0,color=muted,linewidth=.7)
    axes[0].set(xlabel='Time (model units)',ylabel='Spatial minus calibrated count survival',
                title='Localized-start effect across spatial refinements')
    axes[0].legend(facecolor=bg,labelcolor=ink)
    rates=np.array(fine['calibration_rates'])
    b=np.arange(1,len(rates)+1)
    axes[1].plot(b,rates,'o-',color=colors['uniform'],label='Forward-mean matched coefficient')
    axes[1].axhline(fine['intrinsic_average_rate'],color=colors['reaction'],linestyle='--',
                    label='Fast-mixing spatial-average coefficient')
    axes[1].set(xlabel='Fuel before forward repair b',ylabel='Coefficient per fuel molecule',
                title='Calibration retains equilibrium ratios')
    axes[1].legend(facecolor=bg,labelcolor=ink,fontsize=8)
    for extension in ['png','pdf']:
        fig.savefig(ROOT/f'results/figures/spatial-repair-refinement.{extension}',dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--grids',type=int,nargs='+',default=[32,48,64,96])
    parser.add_argument('--steps',type=int,default=1600)
    parser.add_argument('--output',default='results/spatial-repair.json')
    args=parser.parse_args()
    if len(args.grids)<3 or args.grids != sorted(set(args.grids)):
        parser.error('Supply at least three distinct increasing even grids.')
    model=SpatialRepair()
    times=np.linspace(0,4,args.steps+1)
    results=[]
    for n in args.grids:
        start=time.perf_counter()
        grid=build_grid(model,n)
        calibration_diagnostics={}
        rates=calibrate(model,grid,diagnostics=calibration_diagnostics)
        reference=np.full(model.fuel+1,rates[model.fuel])
        diagnostics={}
        all_curves=spatial_survivals(model,grid,times,[(1,model.fuel),(2,model.fuel)],diagnostics=diagnostics)
        cases={}
        for name,h in [('HD',1),('HH',2)]:
            initial=(h,model.fuel)
            curves=all_curves[initial]
            edge=count_survival(model,times,rates,initial=initial)
            ref=count_survival(model,times,reference,initial=initial)
            summaries={}
            for law,values in curves.items():
                gap=np.abs(np.array(values)-edge)
                # Each survival derivative belongs to [-eta,0]. Consequently
                # their difference is eta-Lipschitz, not 2*eta-Lipschitz.
                mesh=float(gap.max())
                summaries[law]={'at_time_four':float(values[-1]),
                    'sampled_max_abs_gap':mesh,
                    'signed_gap_at_sampled_max':float((np.asarray(values)-edge)[gap.argmax()]),
                    'time_of_sampled_max':float(times[gap.argmax()]),
                    'time_continuum_gap_upper_at_this_spatial_grid':
                        min(1.,mesh+model.eta*float(np.diff(times).max())/2)}
            cases[name]={'spatial':{k:np.asarray(v).tolist() for k,v in curves.items()},
                'count_edgewise':edge.tolist(),'count_reference':ref.tolist(),
                'comparisons_to_edgewise':summaries}
        record={'n':n,'spatial_orbits':len(grid.weight),
            'mean_gate':float(grid.weight@grid.gate),
            'intrinsic_average_rate':float(model.activity*(grid.weight@grid.gate)),
            'exponential_solver_diagnostics':diagnostics,
            'calibration_diagnostics':calibration_diagnostics,
            'calibration_rates':rates[1:].tolist(),
            'initial_conditions':cases,'seconds':time.perf_counter()-start}
        results.append(record)
        print(json.dumps({'completed_grid':n,'seconds':record['seconds'],
                          'HD':cases['HD']['comparisons_to_edgewise']}),flush=True)
    refinement=[]
    for coarse,fine in zip(results[:-1],results[1:]):
        differences={name:{law:float(np.max(np.abs(np.array(
            coarse['initial_conditions'][name]['spatial'][law])-
            fine['initial_conditions'][name]['spatial'][law])))
            for law in ['uniform','localized','reaction']} for name in ['HD','HH']}
        refinement.append({'grids':[coarse['n'],fine['n']],
                           'max_curve_differences':differences})
    record={'generated_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Two preassembled complexes; continuous Brownian maintenance model, approximated by symmetry-lumped lattice equations. No self-assembly or new theorem claim.',
        'parameters':asdict(model),'units':'Dimensionless model units; no biological calibration.',
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'calibration':'Each reversible fuel edge uses the inverse uniform-start forward reaction mean divided by b. Reverse partner uses the same coefficient; no survival fitting.',
        'times':times.tolist(),'refinements':results,'refinement_differences':refinement,
        'source_sha256':hashes()}
    (ROOT/args.output).write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    figures(record)


if __name__=='__main__':main()
