"""Export every new rational f in one index, in the article's conventions."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
CASES=[351,375,405,441,476,540,550,575,703,1344]
def load(path):
    spec=importlib.util.spec_from_file_location('audit_export',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def graph_bits(n,mask):
    return ''.join(str((mask>>(j*(j-1)//2+i))&1) for i in range(n) for j in range(i+1,n))
graphs=[]
for line in (ROOT/'article_replay/pure405/atlas/graphs.txt').read_text().splitlines():
    i,n,mask=map(int,line.split());assert i==len(graphs);graphs.append((n,mask))
index=[];all_f=[]
for v in CASES:
    folder=ROOT/'fresh'/str(v);c=json.loads((folder/'certificate.json').read_text());audit=json.loads((folder/'verification.json').read_text())
    assert audit.get('verified') is True or audit.get('status')=='PASS'
    structural=v in [375,441,575];p=7 if v in [375,405,441,575] else 6 if v in [351,476] else 5
    vectors=[];aux=[]
    if p<7:
        for t in c['terms']:
            f=t['f'];vectors.append(dict(H_bits=t['H'],anchor_size=(1+int((1+8*len(t['H']))**.5))//2,degree=t['ell'],f=f,eta=t['eta']))
        aux=[dict(kind='nonnegative_labelled_count',graph_index=i,weight=t) for i,t in enumerate(c['z']) if F(t)]
        delta=F(1);margin=None;count_sources=len(aux);moment_sources=0
    else:
        delta=F(c['exact_margin']);margin=float(delta);v0,k,lam,mu=c['parameters'];rho=F(k+3,v);d={375:34,405:30,441:44,575:46}[v]
        if structural:
            mod=load(ROOT/'article_replay/legacy_ten/verify_psd_completion.py');a=mod.Audit(v)
            core=[0]+json.loads((ROOT/f'article_replay/legacy_ten/rebuild_{v}/exact_elimination_{v}.json').read_text())['free_core_ids']
        else:core=json.loads((ROOT/'article_replay/pure405/data/degree2/core_recipes_405.json').read_text())['cores']
        count_sources=moment_sources=0
        for t in c['terms']:
            s=t['source'];weight=F(t['weight'])/F(t['scale'])/delta
            if s['kind']=='PSD':
                n,mask=graphs[s['anchor']];ell=s['degree'];D=d-n
                if structural:profiles=a.geometry(s['anchor'])[1]
                else:profiles=list(product((0,1),repeat=n))
                f=[0]*(2**n)
                for pattern,value in zip(profiles,s['direction']):f[sum(int(bit)*2**(n-1-i) for i,bit in enumerate(pattern))]=value
                chi=[F(1),3-rho,(3-rho)**2*F(D-1,D)][ell] if structural else F(1)
                eta=weight*F(2**(n*(n+1)//2),prod(range(v-n,v+1)))*chi
                assert eta>0
                vectors.append(dict(H_bits=graph_bits(n,mask),anchor_size=n,atlas_id=s['anchor'],atlas_mask=mask,degree=ell,f=list(map(str,f)),eta=str(eta)))
            else:
                if s['kind']=='nonnegative_count':count_sources+=1
                else:moment_sources+=1
                aux.append(dict(kind=s['kind'],source=s,normalized_weight=str(weight),nonnegative_form='minus full_row dot y',full_row=t['full_row']))
        # Convert the exact bounded-residual comparison to the normalized -1 identity.
        combined=list(map(F,c['combined_coefficients']));check=[-q/delta for q in combined]
        for j,q in enumerate(combined[1:],1):
            if not q:continue
            if q>0:
                aux.append(dict(kind='residual_nonnegative_scaled_count',atlas_id=core[j],weight=str(q/delta)))
                check[j]+=q/delta
            else:
                U=F(c['bounds'][j]);aux.append(dict(kind='residual_upper_bound_slack',atlas_id=core[j],U=str(U),weight=str(-q/delta),nonnegative_form='U - y_F'))
                check[0]+=-q*U/delta;check[j]+=q/delta
        assert check==[F(-1)]+[F(0)]*(len(check)-1)
    form=dict(parameters=c['parameters'],count_order=p,structural_prerequisites=structural,
              profile_order='binary lexicographic, 0..2^m-1; omitted profiles extended by zero',
              H_bits_order='01,02,...,0(m-1),12,...',normalization='article S; labelled N_F; sum nonnegative sources = -1',
              scaled_count_definition='y_F = 2^(h*(h-1)/2) N_F / ((v)_h), h=|F|; (v)_h=v(v-1)...(v-h+1)',
              kernel_terms=vectors,auxiliary_nonnegative_terms=aux,original_exact_margin=str(delta),
              boundary='For structural cases, the article root/complete-shell theorems justify source equations; not a formal proof-assistant encoding.')
    (folder/'article_form.json').write_text(json.dumps(form,indent=2)+'\n');all_f.append(form)
    index.append(dict(parameters=c['parameters'],status='EXACT_VERIFIED',fresh_cone_sdp=True,count_order=p,
                      proof_type='structural' if structural else 'core',kernel_terms=len(vectors),count_sources=count_sources,
                      integer_moment_sources=moment_sources,positive_margin_decimal=margin,
                      certificate=f'fresh/{v}/certificate.json',article_form=f'fresh/{v}/article_form.json',verification=f'fresh/{v}/verification.json'))
(ROOT/'all_weight_functions.json').write_text(json.dumps(all_f,indent=2)+'\n')
(ROOT/'certificate_index.json').write_text(json.dumps(index,indent=2)+'\n')
print(json.dumps(index,indent=2))
