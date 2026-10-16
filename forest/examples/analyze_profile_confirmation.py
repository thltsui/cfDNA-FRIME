"""Audit fresh confirmations and assess conditional sampling sensitivity.

Read bootstrap and pseudo-data recovery use a fixed, complete Beta ABC bank.
Neither is a test of a universal exponent or a calibrated biological interval.
"""
from pathlib import Path
import argparse,csv,hashlib,json,os,sys
from datetime import datetime,timezone
from time import perf_counter
import numpy as np

EPS=(.05,.10,.15)

def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvout(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def quantiles(x):
    x=np.asarray(x);x=x[np.isfinite(x)]
    return np.quantile(x,[.025,.5,.975]).tolist() if len(x) else None
def cdfs(reads):
    reads=np.asarray(reads);depths=reads.sum(axis=1);good=depths>0
    out=np.zeros(reads.shape,dtype=float);out[good]=reads[good].cumsum(axis=1)/depths[good,None]
    return out,good


def refit(bank_reads, targets, alphas):
    """Exact fixed-bin rejection for each target; unavailable fits remain NaN."""
    bank,bgood=cdfs(bank_reads);target,tgood=cdfs(targets)
    shape=(len(target),len(EPS))
    result={k:np.full(shape,np.nan) for k in ('negative_fraction','alpha025','median','alpha975')}
    result['accepted']=np.zeros(shape,dtype=np.int64)
    for first in range(0,len(target),16):
        d=np.max(np.abs(target[first:first+16,None,:]-bank[None,:,:]),axis=2)
        d[:,~bgood]=np.inf;d[~tgood[first:first+16]]=np.inf
        for j,row in enumerate(d):
            for e,epsilon in enumerate(EPS):
                a=alphas[row<=epsilon];i=first+j
                result['accepted'][i,e]=len(a)
                if len(a):
                    result['negative_fraction'][i,e]=np.mean(a<0)
                    q=np.quantile(a,[.025,.5,.975])
                    for k,v in zip(('alpha025','median','alpha975'),q):result[k][i,e]=v
    return result


def variance_components(histograms,reads,depth):
    population,good=cdfs(histograms);read,readgood=cdfs(reads)
    assert np.array_equal(good,readgood)
    population=population[good];read=read[good]
    if len(population)<2:return None
    latent=population.var(axis=0,ddof=1)
    sampling=np.mean(population*(1-population)/depth,axis=0)
    return dict(latent=latent,sampling=sampling,total=latent+sampling,empirical=read.var(axis=0,ddof=1))


def setup(root):
    plan=load(root/'plan.json');ready=load(root/'prepare-complete.json')
    assert sha(root/'plan.json')==ready['plan_sha256']
    for name,h in plan['input_sha256'].items():assert sha(Path(name))==h,name
    screen=Path(plan['screen_directory']);spec=load(screen/'plan.json')
    bank=Path(plan['bank_directory']);rows=load(bank/'plan.json')['rows']
    for name,h in spec['bank_sha256'].items():assert sha(bank/name)==h,name
    ix=np.array([r['case'] in ('positive_beta','negative_beta') for r in rows])
    assert int(ix.sum())==2000
    assert sum(r['case']=='positive_beta' for r in rows)==1000
    alpha=np.array([r['theta']['alpha'] for r in rows])[ix]
    with np.load(screen/'screen.npz',allow_pickle=False) as z:
        reads=z['reads'][:,ix];names=z['profiles'].tolist();obs=z['observed'];bins=z['bins_bp']
    return plan,reads,names,obs,alpha,bins


def bootstrap(root):
    started=datetime.now(timezone.utc).isoformat();start=perf_counter()
    plan,bank,names,observed,alpha,bins=setup(root)
    forest=root.parents[1];sys.path.insert(0,str(forest/'examples'))
    from frime_forest.trees import rng_for
    out=root/'read-bootstrap';out.mkdir(exist_ok=False)
    nboot=plan['analysis']['read_bootstrap_repeats'];all_results=[];all_reads=[];summaries=[]
    for profile in plan['profiles']:
        j=names.index(profile);obs=observed[j];n=int(obs.sum())
        reads=rng_for(plan['seed'],j,902).multinomial(n,obs/n,size=nboot)
        result=refit(bank[j],reads,alpha);all_results.append(result);all_reads.append(reads)
        for e,epsilon in enumerate(EPS):
            median=result['median'][:,e];valid=np.isfinite(median)
            summaries.append(dict(profile=profile,depth=n,epsilon=epsilon,repeats=nboot,
                available=int(valid.sum()),no_acceptance=int((~valid).sum()),sparse=int(np.sum((result['accepted'][:,e]>0)&(result['accepted'][:,e]<20))),
                median_alpha025_50_975=quantiles(median),negative_median_count=int(np.sum(median<0)),
                positive_median_count=int(np.sum(median>0)),negative_fraction025_50_975=quantiles(result['negative_fraction'][:,e]),
                accepted_count025_50_975=quantiles(result['accepted'][:,e])))
    arrays={k:np.array([r[k] for r in all_results]) for k in all_results[0]}
    group=[]
    for e,epsilon in enumerate(EPS):
        med=arrays['median'][:,:,e];valid=np.all(np.isfinite(med),axis=0)
        gap=med[:2].mean(axis=0)-med[2:].mean(axis=0)
        group.append(dict(epsilon=epsilon,repeats=nboot,available=int(valid.sum()),unavailable=int((~valid).sum()),
            positive_gap_count=int(np.sum(gap[valid]>0)),gap025_50_975=quantiles(gap),
            interpretation='Conditional read-bootstrap distribution, not a biological confidence interval or hypothesis test.'))
    np.savez_compressed(out/'bootstrap.npz',**arrays,reads=np.array(all_reads),profiles=np.array(plan['profiles']))
    dump(out/'summary.json',summaries);dump(out/'group-gap.json',group)
    dump(out/'manifest.json',dict(script_sha256=sha(Path(__file__)),plan_sha256=sha(root/'plan.json'),
        started_utc=started,wall_seconds=perf_counter()-start,numpy=np.__version__,
        seed=plan['seed'],component=902,repeats_per_profile=nboot,
        interpretation='Resampling observed reads only, conditional on empirical histograms and fixed simulation bank. Not new FRIME populations.',
        artifacts_sha256={name:sha(out/name) for name in ('bootstrap.npz','summary.json','group-gap.json')}))
    print(json.dumps(dict(stage='read-bootstrap',groups=group)),flush=True)


def audit(root,plan,reads,hist,indices,distances):
    done=load(root/'complete.json')
    assert not (root/'FAILED.json').exists()
    for name,h in done['artifacts_sha256'].items():assert sha(root/name)==h,name
    expected=[(c['cell'],j) for c in plan['cells'] for j in range(plan['repeats'])]
    assert indices.tolist()==[list(x) for x in expected]
    seeds=set();streams=set();families=nodes=0;bound=0.;count=0
    with (root/'runs.jsonl').open() as f:
        for i,line in enumerate(f):
            r=json.loads(line);count+=1;cell=plan['cells'][r['cell']]
            assert (r['cell'],r['repeat'])==expected[i]
            assert r['profile']==cell['profile'] and r['model']==cell['model'] and r['proposal']==cell['id']
            assert r['certificate']==plan['certificates'][r['cell']]['certificate']
            assert np.array_equal(r['reads'],reads[i]) and np.array_equal(r['histogram'],hist[i])
            assert hist[i].sum()==r['observable_count']
            obs=np.array(plan['observed'][r['profile']]);depth=int(reads[i].sum())
            assert depth in (0,int(obs.sum()))
            d=float(np.max(np.abs(reads[i].cumsum()/depth-obs.cumsum()/obs.sum()))) if depth else np.inf
            assert d==distances[i] and d==(r['distance'] if r['distance'] is not None else np.inf)
            seeds.add(r['seed']);streams.add(tuple(r['observation_stream']))
            families+=r['immigrants'];nodes+=r['generated_nodes'];bound+=r['certificate']['lambda_bound']
    assert count==len(expected)==done['populations']==len(seeds)==len(streams)
    assert families==done['immigrants'] and nodes==done['generated_nodes'] and bound==done['history_union_bound']
    old=Path(plan['bank_directory'])
    oldseeds=set();oldstreams=set()
    for name in ('screen-runs.jsonl','confirm-runs.jsonl'):
        with (old/name).open() as f:
            for line in f:
                r=json.loads(line);oldseeds.add(r['seed']);oldstreams.add(tuple(r['observation_stream']))
    assert not seeds & oldseeds and not streams & oldstreams
    return dict(populations=count,unique_new_population_streams=len(seeds),unique_new_read_streams=len(streams),
        disjoint_from_old_screen_and_confirmation=True,all_histograms_reads_distances_models_certificates_checked=True,
        immigrants=families,generated_nodes=nodes,history_union_bound=bound)


def confirm(root):
    started=datetime.now(timezone.utc).isoformat();start=perf_counter()
    plan,bank,names,observed,alpha,bins=setup(root)
    with np.load(root/'confirm.npz',allow_pickle=False) as z:
        hist=z['histograms'];reads=z['reads'];indices=z['indices'];distances=z['distances']
    checked=audit(root,plan,reads,hist,indices,distances)
    bmanifest=load(root/'read-bootstrap/manifest.json')
    assert bmanifest['plan_sha256']==sha(root/'plan.json')
    for name,h in bmanifest['artifacts_sha256'].items():assert sha(root/'read-bootstrap'/name)==h
    out=root/'analysis';out.mkdir(exist_ok=False)
    summaries=load(root/'summary.json');candidate_rows=[];recovery_rows=[];variance_rows=[];recoveries=[]
    variances={}
    for cell in plan['cells']:
        ix=indices[:,0]==cell['cell'];j=names.index(cell['profile']);n=int(observed[j].sum())
        s=next(s for s in summaries if s['cell']==cell['cell']);d=distances[ix]
        f,good=cdfs(reads[ix]);gap=np.abs(f[good].mean(axis=0)-observed[j].cumsum()/n) if np.any(good) else None
        if gap is not None:assert np.isclose(gap.max(),s['mean_cdf_max_gap'],rtol=0,atol=1e-14)
        row=dict(profile=cell['profile'],case=cell['case'],rank=cell['rank'],proposal=cell['id'],alpha=cell['theta']['alpha'],
            mean_cdf_gap=s.get('mean_cdf_max_gap'),median_distance=float(np.median(d[np.isfinite(d)])) if np.any(np.isfinite(d)) else None,
            no_reads=int(np.sum(~good)))
        for e in EPS:
            item=s['successes'][str(e)];assert item['count']==int(np.sum(d<=e))
            key=f'{e:.2f}'
            row.update({f'count_{key}':item['count'],f'rate_{key}':item['proportion'],f'low_{key}':item['wilson95'][0],f'high_{key}':item['wilson95'][1]})
        candidate_rows.append(row)
        vc=variance_components(hist[ix],reads[ix],n);variances[cell['cell']]=vc
        if vc is not None:
            worst=int(np.argmax(gap));v=vc['total'][worst]
            variance_rows.append(dict(profile=cell['profile'],case=cell['case'],rank=cell['rank'],proposal=cell['id'],
                evaluated_at_bp=float(bins[1:][worst]-.5),latent_variance=float(vc['latent'][worst]),
                read_variance=float(vc['sampling'][worst]),predicted_variance=float(v),
                empirical_readout_variance=float(vc['empirical'][worst]),
                read_fraction=float(vc['sampling'][worst]/v) if v>0 else None))
        recovery=refit(bank[j],reads[ix],alpha);recoveries.append(recovery)
        truth=cell['theta']['alpha']
        for e,epsilon in enumerate(EPS):
            med=recovery['median'][:,e];valid=np.isfinite(med)
            correct=(med<0) if truth<0 else (med>0)
            covers=(recovery['alpha025'][:,e]<=truth)&(recovery['alpha975'][:,e]>=truth)
            recovery_rows.append(dict(profile=cell['profile'],case=cell['case'],rank=cell['rank'],proposal=cell['id'],true_alpha=truth,
                epsilon=epsilon,repeats=len(med),available=int(valid.sum()),no_acceptance=int((~valid).sum()),
                sparse=int(np.sum((recovery['accepted'][:,e]>0)&(recovery['accepted'][:,e]<20))),
                correct_median_sign_count=int(correct.sum()),opposite_median_sign_count=int(np.sum(valid&~correct)),
                median_alpha025=quantiles(med)[0] if valid.any() else None,
                median_alpha_median=quantiles(med)[1] if valid.any() else None,
                median_alpha975=quantiles(med)[2] if valid.any() else None,
                true_alpha_in_sample_quantile_range_count=int(covers.sum())))
    csvout(out/'candidates.csv',candidate_rows);csvout(out/'recovery.csv',recovery_rows);csvout(out/'variance.csv',variance_rows)
    np.savez_compressed(out/'recovery.npz',**{k:np.array([r[k] for r in recoveries]) for k in recoveries[0]})
    primary=[r for r in candidate_rows if r['rank']==1]
    dump(out/'primary.json',primary);dump(out/'audit.json',dict(**checked,script_sha256=sha(Path(__file__)),started_utc=started,wall_seconds=perf_counter()-start))
    make_plots(root,out,plan,summaries,primary,bins)
    print(json.dumps(dict(stage='analysis',audit=checked,primary=primary)),flush=True)


def make_plots(root,out,plan,summaries,primary,bins):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/frime-mpl-config')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    colors={'positive_beta':'#B65D1B','negative_beta':'#246A9F'}
    fig,ax=plt.subplots(figsize=(10,5.5))
    for j,profile in enumerate(plan['profiles']):
        for case,offset in [('positive_beta',-.12),('negative_beta',.12)]:
            r=next(r for r in primary if r['profile']==profile and r['case']==case)
            q=r['rate_0.10'];low=r['low_0.10'];high=r['high_0.10']
            ax.errorbar(q,j+offset,xerr=[[max(0.,q-low)],[max(0.,high-q)]],fmt='o',color=colors[case],capsize=3)
            ax.text(1.03,j+offset,f"{r['count_0.10']}/100; α={r['alpha']:.2f}",va='center',fontsize=9)
    ax.set_yticks(range(6),plan['profiles']);ax.invert_yaxis();ax.set_xlim(-.03,1.02)
    ax.set_xlabel('Probability of a fresh readout having CDF distance D ≤ 0.10')
    ax.set_title('Fresh populations: candidates fixed by screening rank',pad=48)
    ax.legend([Line2D([0],[0],color=colors[c],marker='o') for c in colors],['Positive α candidate','Negative α candidate'],loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2,frameon=False)
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
    fig.text(.5,.015,'100 independent populations per candidate. Marginal Wilson 95% intervals; conditional fit probabilities, not biological confidence intervals.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,.82,1));fig.savefig(out/'fresh-fit-probabilities.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True,sharey=True)
    x=bins[1:]-.5
    for ax,profile in zip(axes.flat,plan['profiles']):
        obs=np.array(plan['observed'][profile]);ax.step(x,obs.cumsum()/obs.sum(),where='post',color='black',lw=1.8,label='Observed')
        for case in colors:
            s=next(r for r in summaries if r['profile']==profile and r['case']==case and r['rank']==1)
            if 'mean_cdf' in s:
                ax.plot(x,s['mean_cdf'],color=colors[case],label=f"α={s['theta']['alpha']:.2f}")
                ax.fill_between(x,s['pointwise025'],s['pointwise975'],color=colors[case],alpha=.15)
        ax.set_xscale('log');ax.set_xlim(40,16569);ax.set_ylim(0,1)
        ax.set_title(f'{profile} (n={obs.sum():,})');ax.grid(alpha=.12);ax.legend(fontsize=8)
    fig.supxlabel('Fragment length (bp)',y=.05);fig.supylabel('CDF')
    fig.suptitle('Fresh predictive profiles: original screening rank 1 in each sign',fontsize=15)
    fig.text(.5,.01,'Shading: pointwise 95% predictive intervals, not simultaneous bands or calibrated goodness-of-fit tests.',ha='center',fontsize=10)
    fig.tight_layout(rect=(.02,.08,1,.95));fig.savefig(out/'fresh-predictive-cdfs.png',dpi=180,bbox_inches='tight');plt.close(fig)
    bs=load(root/'read-bootstrap/summary.json')
    fig,ax=plt.subplots(figsize=(10,5.5))
    for j,profile in enumerate(plan['profiles']):
        r=next(r for r in bs if r['profile']==profile and r['epsilon']==.10)
        q=r['median_alpha025_50_975']
        if q:
            ax.plot([q[0],q[2]],[j,j],color='#6B4093',lw=2);ax.scatter(q[1],j,color='#6B4093')
        ax.text(1.02,j,f"{r['available']}/1000 fits",transform=ax.get_yaxis_transform(),va='center',fontsize=9)
    ax.set_yticks(range(6),plan['profiles']);ax.invert_yaxis();ax.axvline(0,color='black',ls=':',lw=1)
    ax.set_xlim(-1.5,1.5);ax.set_xlabel('ABC median α after resampling observed reads (D ≤ 0.10)')
    ax.set_title('Sensitivity to observed read sampling, conditional on the fixed ABC bank')
    ax.grid(axis='x',alpha=.15);ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.015,'Points/lines: median and 2.5–97.5% resampling quantiles. Not confidence intervals for biological α; excludes process and bank uncertainty.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,.84,1));fig.savefig(out/'read-resampling-alpha.png',dpi=180,bbox_inches='tight');plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--stage',choices=('bootstrap','confirm'),required=True)
    a=p.parse_args();(bootstrap if a.stage=='bootstrap' else confirm)(a.directory.resolve())

if __name__=='__main__':main()
