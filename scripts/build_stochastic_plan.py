"""Calculate the diffusion feasibility screen and render the implementation plan."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
from math import pi
import os
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from purposeful.encounter import (shell_radius, mean_capture_time, uniform_capture_probability,
                                 capture_before_closure, partial_reactivity_mean)
from build_research_docs import page


def main():
    diffusion, radius = 7.7, .005
    volumes = [.001,.01,1.,1000.]
    rows=[]
    for volume in volumes:
        mean=mean_capture_time(volume,diffusion,radius)
        smol=volume/(4*pi*diffusion*radius)
        exact=uniform_capture_probability(volume,diffusion,radius,1.)
        rows.append({'volume_um3':volume,'outer_radius_um':shell_radius(volume,radius),
                     'mean_capture_seconds':mean,'smoluchowski_seconds':smol,
                     'crossing_scale_seconds':volume**(2/3)/(6*diffusion),
                     'probability_before_1_per_second_closure':exact,
                     'matched_mean_exponential_probability':1/(1+mean)})
    gamma=np.logspace(-3,3,121)
    exact=np.array([uniform_capture_probability(1,diffusion,radius,g) for g in gamma])
    matched=1/(1+gamma*mean_capture_time(1,diffusion,radius))
    near=np.array([capture_before_closure(radius+.001,1,diffusion,radius,g) for g in gamma])
    bg,ink,muted,teal,amber='#0d141d','#e6edf7','#adc0d6','#72dbc1','#edbd75'
    plt.rcParams.update({'figure.facecolor':bg,'axes.facecolor':bg,'text.color':ink,
        'axes.labelcolor':ink,'axes.edgecolor':muted,'xtick.color':muted,'ytick.color':muted,
        'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    axes[0].loglog(gamma,exact,color=teal,label='Exact uniform-start search')
    axes[0].loglog(gamma,matched,'--',color=amber,label='Matched-mean exponential')
    axes[0].loglog(gamma,near,color='#d899e8',label='Start 1 nm outside capture surface')
    axes[0].set(xlabel='Independent window closure rate (1/s)',ylabel='Probability of capture before closure',title='Spatial initialization can matter more than mean search time')
    axes[0].legend(facecolor=bg,labelcolor=ink,fontsize=9)
    axes[1].semilogx(gamma,abs(exact-matched),color=teal)
    axes[1].axhline(.01,color=amber,linestyle=':',label='Proposed 1 percentage point decision threshold')
    axes[1].set(xlabel='Independent window closure rate (1/s)',ylabel='Absolute probability difference',title='Uniform-start comparison after matching the mean')
    axes[1].legend(facecolor=bg,labelcolor=ink,fontsize=8)
    fig.savefig(ROOT/'results/figures/encounter-feasibility.png',dpi=180)
    fig.savefig(ROOT/'results/figures/encounter-feasibility.pdf')
    plt.close(fig)
    ds=np.logspace(-1,2,81)
    vs=np.logspace(-3,3,81)
    means=np.array([[mean_capture_time(v,d,radius) for d in ds] for v in vs])
    fig,ax=plt.subplots(figsize=(8,5.5),layout='constrained')
    image=ax.pcolormesh(ds,vs,np.log10(means),shading='auto',cmap='viridis')
    contour=ax.contour(ds,vs,means,levels=[.01,1,100],colors=amber,linewidths=1)
    ax.clabel(contour,fmt='%g s',fontsize=10)
    ax.scatter([7.7],[1],marker='x',s=90,color='white')
    ax.set(xscale='log',yscale='log',xlabel='Diffusion coefficient (um squared/s)',ylabel='Accessible compartment volume (um cubed)',title='One searching particle and a fixed 5 nm target')
    fig.colorbar(image,ax=ax,label='log10 exact mean capture time (seconds)')
    fig.savefig(ROOT/'results/figures/encounter-timescales.png',dpi=180)
    fig.savefig(ROOT/'results/figures/encounter-timescales.pdf')
    plt.close(fig)
    inputs=['src/purposeful/encounter.py','tests/test_encounter.py','scripts/build_stochastic_plan.py',
            'docs/stochastic-calculus-plan.md','research/stochastic-calculus-prior-art.md']
    result={'generated_utc':datetime.now(timezone.utc).isoformat(),
            'status':'Classical benchmark and implementation plan, with no new novelty claim',
            'assumptions':{'diffusion_um2_per_s':diffusion,'capture_radius_um':radius,
                           'geometry':'fixed central target and reflecting outer sphere',
                           'initial_distribution':'uniform over accessible spherical shell',
                           'closure_rate_per_s':1.,
                           'parameter_provenance':'D is a GFP tracer measurement anchor. Capture radius, volumes and closure rates are assumed sensitivity inputs.'},
            'volume_cases':rows,
            'parameter_sensitivity':{
                'capture_radius_cases':[{'a_um':a,'mean_seconds':mean_capture_time(1,diffusion,a),
                                         'probability_at_gamma_one':uniform_capture_probability(1,diffusion,a,1)}
                                        for a in [.002,.005,.01]],
                'diffusion_anchor_cases':[{'D_um2_per_s':d,'mean_seconds':mean_capture_time(1,d,radius),
                                          'probability_at_gamma_one':uniform_capture_probability(1,d,radius,1)}
                                         for d in [5.2,7.7,10.2]]},
            'baseline':{'diffusion_limited_rate_um3_per_s':4*pi*diffusion*radius,
                        'rate_M_inverse_s_inverse':4*pi*diffusion*radius*6.02214076e8,
                        'one_particle_concentration_nM':1e24/6.02214076e23,
                        'perfect_capture_probability':uniform_capture_probability(1,diffusion,radius,1),
                        'near_start_probability_at_gamma_one':capture_before_closure(.006,1,diffusion,radius,1),
                        'mean_with_kappa_1000_um_per_s':partial_reactivity_mean(1,diffusion,radius,1000),
                        'mean_with_kappa_10_um_per_s':partial_reactivity_mean(1,diffusion,radius,10),
                        'target_scale_brownian_step_seconds':radius**2/(6*diffusion)},
            'probability_sweep':{'gamma_per_s':gamma.tolist(),'exact_uniform':exact.tolist(),
                                 'matched_exponential':matched.tolist(),'localized_start':near.tolist()},
            'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs}}
    (ROOT/'results/encounter-feasibility.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    run=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s','tests','-v'],cwd=ROOT,
                       env={**os.environ,'PYTHONPATH':str(ROOT/'src')},capture_output=True,text=True)
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'return_code':run.returncode,
             'command':'python -B -m unittest discover -s tests -v','test_output':run.stdout+run.stderr,
             'new_test_scope':'MFPT quadrature, independent boundary-value solver, scaling, small-rate limit and invalid-input checks',
             'data_sha256':hashlib.sha256((ROOT/'results/encounter-feasibility.json').read_bytes()).hexdigest()}
    (ROOT/'results/encounter-verification.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    if run.returncode:raise RuntimeError(run.stdout+run.stderr)
    table='\n\n## Calculated feasibility cases\n\n| Accessible volume (um cubed) | Exact mean search (s) | Smoluchowski estimate (s) | Exact capture before closure | Matched-mean exponential |\n|---|---:|---:|---:|---:|\n'
    for row in rows:
        table+=f"| {row['volume_um3']:g} | {row['mean_capture_seconds']:.6g} | {row['smoluchowski_seconds']:.6g} | {row['probability_before_1_per_second_closure']:.8g} | {row['matched_mean_exponential_probability']:.8g} |\n"
    content=(ROOT/'docs/stochastic-calculus-plan.md').read_text(encoding='utf-8')+table
    content+='\n\n[Machine-readable calculations](../results/encounter-feasibility.json) · [Full-suite verification receipt](../results/encounter-verification.json)\n'
    page('STOCHASTIC-CALCULUS-PLAN.html',content,'Implementation plan — diffusion, repair, and stochastic calculus')
    print(json.dumps({'volume_cases':rows,'baseline':result['baseline'],'tests':run.stderr[-120:]},indent=2))


if __name__=='__main__':main()
