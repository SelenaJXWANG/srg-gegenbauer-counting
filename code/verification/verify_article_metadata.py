#!/usr/bin/env python3
"""Cross-check the article's parameter and structural summary data.

Standard library only; no inputs or reports are changed. This is an alignment
check, supplementary to the independent exact source-row reconstructions.
It reads the actual TeX parameter/inventory/quartic tables and strict margins,
then compares them with the declared cases and retained exact certificate data.
It does not certify arbitrary mathematical prose or the structural theorem.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import sys

# Exact certificate integers exceed Python 3.11's default decimal limit.
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = [
    (351,140,73,44), (375,102,45,21), (405,132,63,33),
    (441,88,35,13), (476,133,60,28), (540,147,66,30),
    (550,162,75,36), (575,112,45,16), (703,182,81,35),
    (1344,221,88,26),
]
STRUCTURAL = {375,441,575}
PREVIOUS = {351,550,703,1344}


def require(value, message):
    if not value:
        raise ValueError(message)


def unique_block(source, label):
    blocks = re.findall(r'\\begin\{table\}.*?\\end\{table\}', source, re.S)
    blocks = [block for block in blocks if r'\label{'+label+'}' in block]
    require(len(blocks)==1, 'Missing or duplicated table '+label)
    return blocks[0]


def cells(line):
    return [cell.strip().strip('$') for cell in
            re.sub(r'\\\\\s*$', '', line.strip()).split('&')]


def parameter_table(source):
    rows=[]
    for line in unique_block(source,'tab:parameters').splitlines():
        if not re.match(r'\s*\$?\d+',line):
            continue
        row=cells(line)
        require(len(row)==8,'Wrong parameter-table column count')
        match=re.fullmatch(r'(\d+)(\^\\dagger)?',row[0])
        require(match is not None,'Invalid parameter-table order')
        v=int(match[1]); par=tuple([v]+list(map(int,row[1:4])))
        r,d,p=map(int,row[4:7]); _,k,lam,mu=par
        require((v-k-1)*mu==k*(k-lam-1),'Invalid SRG parameter relation')
        require(k==3*r+mu and r==lam-mu+3,'Wrong restricted eigenvalue')
        require(Q(3*(v-1)-k,r+3)==d,'Wrong eigenspace dimension')
        require(p==(7 if v in STRUCTURAL|{405} else 6 if v in {351,476} else 5),
                'Wrong displayed certificate order')
        require(bool(match[2])==(v in PREVIOUS),'Changed Liao–Yu marker')
        expected=(r'Used (Section~\ref{sec:structural})' if v in STRUCTURAL else 'Not needed')
        require(row[7]==expected,'Wrong additional-structure classification')
        rows.append(par)
    require(rows==PARAMETERS,'Parameter table is missing, duplicated, or reordered')
    # Cross-check available explicit parameter tuples in the exact inputs.
    for par in PARAMETERS:
        v=par[0]
        if v in STRUCTURAL or v==1344:
            continue
        if v==405:
            path=ROOT/'data/pure405/data/degree2/pure_exact_certificate.json'
        else:
            path=ROOT/'data/legacy_ten/six_exact_data'/('_'.join(map(str,par))+'.raw.json')
        require(tuple(json.loads(path.read_text())['parameters'])==par,
                'Parameter-table tuple differs from certificate data')
    return rows


def structural_data(source):
    folder=ROOT/'data/legacy_ten'
    graphs={i:(n,mask) for i,n,mask in
            (map(int,line.split()) for line in (folder/'atlas_graphs.txt').read_text().splitlines())}
    records={v:json.loads((folder/f'rebuild_{v}/recovered_exact_certificate.json').read_text())
             for v in sorted(STRUCTURAL)}
    printed=[]
    for line in unique_block(source,'tab:struct-certificates').splitlines():
        if re.match(r'\s*\d+\s*&',line):
            row=list(map(int,cells(line)))
            require(len(row)==5,'Wrong structural inventory column count')
            printed.append(row)
    actual=[];quartics=[];details={}
    for v,record in records.items():
        require(record['v']==v,'Wrong structural certificate order')
        inventory=Counter(row['source']['kind'] for row in record['rows'])
        require(set(inventory)=={'nonnegative_count','PSD','integer_moment'},
                'Unexpected structural source kind')
        require(all(Q(row['weight'])>0 and Q(row['normalization'])>0
                    for row in record['rows']), 'Nonpositive source multiplier or row scale')
        actual.append([v,inventory['nonnegative_count'],inventory['PSD'],
                       inventory['integer_moment'],len(record['rows'])])
        moments=json.loads((folder/f'rebuild_{v}/integer_moments.json').read_text())
        for row in record['rows']:
            if row['source']['kind']!='integer_moment':
                continue
            index=row['source']['row_index']
            matches=[w for w in moments['witnesses'] if w['kind']=='ineq' and w['row_index']==index]
            require(bool(matches),'Quartic has no integer-moment witness')
            w=matches[0]; a,aa,b,bb=w['nodes']
            require(aa==a+1 and bb==b+1 and w['sign']==1,'Wrong quartic root pattern')
            n,mask=graphs[w['anchor']]
            require(n==3,'Quartic anchor does not have three vertices')
            edges=bin(mask).count('1')
            require(edges in (0,2,3),'Unexpected quartic anchor type')
            name={0:r'\overline{K_3}',2:r'P_3',3:r'K_3'}[edges]
            quartics.append((v,name,a,b))
        coefficients=list(map(Q,record['combined_coefficients']))
        upper=record['variable_upper_bound']
        require(upper==(1024 if v==375 else 256),'Wrong count upper bound')
        require(len(coefficients)==record['dimension']==111,'Wrong structural free dimension')
        margin=coefficients[0]+upper*sum((c for c in coefficients[1:] if c<0),Q(0))
        require(margin==Q(record['exact_margin']),'Saved structural margin mismatch')
        pattern=(r'\\Delta_\{'+str(v)+r'\}\s*>\s*\\frac'
                 r'(?:\{(\d+)\}|(\d+))\{10\^\{?(\d+)\}?\}')
        bounds=re.findall(pattern,source)
        require(len(bounds)==1,'Missing or duplicated structural margin for '+str(v))
        numerator,short,power=bounds[0];bound=Q(int(numerator or short),10**int(power))
        require(0<bound<margin,'Printed strict structural margin is false')
        details[str(v)]={'inventory':dict(inventory),'sources':len(record['rows']),
                         'printed_lower_bound':str(bound),'exact_margin':str(margin)}
    require(printed==actual,'Printed structural source inventory differs from data')
    # Identify the quartic table by its rows, independently of surrounding prose.
    candidate=[]
    for block in re.findall(r'\\begin\{tabular\}.*?\\end\{tabular\}',source,re.S):
        rows=[]
        for line in block.splitlines():
            row=cells(line)
            if len(row)==4 and row[0].isdigit() and int(row[0]) in STRUCTURAL:
                require(row[2].isdigit() and row[3].isdigit(),'Invalid printed quartic nodes')
                rows.append((int(row[0]),row[1],int(row[2]),int(row[3])))
        if rows:
            candidate.append(rows)
    require(len(candidate)==1 and candidate[0]==quartics,
            'Missing, duplicated, or changed structural quartic table')
    return details,quartics


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tex',type=Path,default=ROOT/'paper/srg_counting_gegenbauer_article_1014.tex')
    args=parser.parse_args()
    source=args.tex.read_text(encoding='utf-8')
    parameters=parameter_table(source)
    structural,quartics=structural_data(source)
    print(json.dumps({'verified':True,'input':args.tex.name,
          'scope':'Ten parameter-table rows, spectra, orders and structural classification; three structural inventories, six quartics and strict rational margins.',
          'parameters':parameters,'structural':structural,'quartics':quartics,
          'arithmetic':'fractions.Fraction; standard library only'},indent=2))


if __name__=='__main__':
    main()
