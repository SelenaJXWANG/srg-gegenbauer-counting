#!/usr/bin/env python3
"""Verify certificate data ACTUALLY PRINTED in the manuscript and 405 supplement.

Standard library only. The default manuscript is in the paper/ directory; --tex selects another copy. The required 405 supplement
defaults to srg_counting_gegenbauer_article_1014_405_supplement.tex beside the
chosen manuscript; --supplement selects another copy. Checks all six-point
integer weights and directions in the manuscript, all 405 source weights and
directions in the supplement, an optional main-text worked example, and both
documents' printed strict bounds using exact reconstructed coefficients.
Duplicated 405 tables or source indices are rejected. Both actual input paths
are reported. This supplements, rather than replaces, the independent
atlas/kernel replay. No input or existing output file is modified.
"""
import argparse
from fractions import Fraction as Q
import importlib.util
import json
from math import gcd, lcm
from pathlib import Path
import re
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent
BASE = 10**40
BIG_PATTERN = r'\\srgBigInt\{(\d+)\}\{(\d+)\}\{(\d+)\}'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def table(source, label):
    blocks = re.findall(r'\\begin\{longtable\}.*?\\end\{longtable\}', source, re.S)
    hits = [block for block in blocks if r'\label{' + label + '}' in block]
    require(len(hits) == 1, 'Missing or duplicated table: ' + label)
    return hits[0]


def integer_stack(value):
    match = re.fullmatch(r'\s*\$?\s*' + BIG_PATTERN + r'\s*\$?\s*', value)
    require(match is not None, 'Malformed three-part integer: ' + value[:80])
    a, b, c = map(int, match.groups())
    require(all(0 <= n < BASE for n in (a,b,c)), 'Integer block exceeds 10^40')
    return a*BASE**2 + b*BASE + c


def fraction(value):
    match = re.fullmatch(r'\s*\$?\s*\\(?:frac|dfrac|tfrac)\{(-?\d+)\}\{(\d+)\}\s*\$?\s*', value)
    require(match is not None, 'Malformed fraction: ' + value[:80])
    require(int(match[2]) > 0, 'Fraction denominator not positive')
    return Q(int(match[1]), int(match[2]))


def without_row_end(line):
    return re.sub(r'\\\\(?:\[[^\]]*\])?\s*$', '', line.strip())


def array_values(fragment):
    marker = r'\begin{array}'
    require(fragment.count(marker) == 1, 'Expected exactly one direction array')
    start = fragment.index(marker) + len(marker)
    while fragment[start].isspace():
        start += 1
    require(fragment[start] == '{', 'Array has no column specification')
    depth, end = 1, start+1
    while depth:
        require(end < len(fragment), 'Unclosed array column specification')
        depth += (fragment[end] == '{') - (fragment[end] == '}')
        end += 1
    finish = fragment.index(r'\end{array}', end)
    content = fragment[end:finish]
    require(re.fullmatch(r'[\s\d&\\+-]*', content) is not None,
            'Unsupported token in printed direction array')
    return list(map(int, re.findall(r'-?\d+', content)))


def sixpoint(source):
    module_path = HERE.parent/'data/legacy_ten/six_exact_data/verify_raw_certificate.py'
    spec = importlib.util.spec_from_file_location('expanded_raw_verifier', module_path)
    raw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(raw)
    result = {}
    for v, den, stem in ((351,30,'351_140_73_44'), (476,10,'476_133_60_28')):
        data = json.loads((module_path.parent/(stem+'.raw.json')).read_text())
        raw_report = raw.verify(data)
        require(raw_report['verified'], 'Raw identity did not verify')
        terms = data['nonnegative_terms']
        kernels = [term for term in terms if term['kind'] == 'sos']
        counts = [term for term in terms if term['kind'] == 'count']
        require(len(kernels)+len(counts) == 24, 'Six-point source inventory changed')

        constants = re.findall(r'L_\{'+str(v)+r'\}\s*=\s*('+BIG_PATTERN+r')', source)
        require(len(constants) == 1, 'Missing or duplicated L_'+str(v))
        L = integer_stack(constants[0][0])
        require(L > 0, 'Nonpositive common scale')
        printed_kernel, printed_count = {}, {}
        block = table(source, 'tab:six'+str(v)+'-weights')
        for line in block.splitlines():
            if not re.match(r'^\s*\d+\s*&', line):
                continue
            fields = without_row_end(line).split('&')
            require(len(fields) == 5, 'Wrong number of six-point weight columns')
            i = int(fields[0]); require(i not in printed_kernel, 'Duplicate integer kernel multiplier')
            printed_kernel[i] = integer_stack(fields[1])
            if fields[2].strip():
                j = int(fields[2]); require(j not in printed_count, 'Duplicate integer count multiplier')
                bits = re.fullmatch(r'\s*\\texttt\{([01]{15})\}\s*', fields[3])
                require(bits is not None, 'Invalid printed six-vertex graph')
                printed_count[j] = (bits[1], integer_stack(fields[4]))
            else:
                require(not fields[3].strip() and not fields[4].strip(), 'Unindexed count data')
        require(set(printed_kernel) == set(range(1,len(kernels)+1)), 'Missing kernel integer weights')
        require(set(printed_count) == set(range(1,len(counts)+1)), 'Missing count integer weights')
        ratios = []
        for i,term in enumerate(kernels,1):
            aut = raw.aut(term['anchor_size'],raw.read_graph(term['anchor_size'],term['anchor_bits']))
            target = Q(term['weight'])/aut
            require(printed_kernel[i] > 0 and Q(printed_kernel[i],L) == target,
                    'Printed integer kernel multiplier fails for '+str(v)+', row '+str(i))
            ratios.append(target)
        for j,term in enumerate(counts,1):
            aut = raw.aut(term['size'],raw.read_graph(term['size'],term['bits']))
            bits,weight = printed_count[j]
            target = Q(term['weight'])/aut
            require(bits == term['bits'] and weight > 0 and Q(weight,L) == target,
                    'Printed integer count multiplier fails for '+str(v)+', row '+str(j))
            ratios.append(target)
        # Also confirm the actual printed L is the specified primitive common scale.
        expected_L = lcm(*(q.denominator for q in ratios))
        div = gcd(expected_L, *(int(q*expected_L) for q in ratios))
        require(L == expected_L//div, 'Printed common scale L differs from exact source')

        directions = {}
        block = table(source,'tab:direct-directions-'+str(v))
        require(re.search(r'divide each displayed vector by \$'+str(den)+r'\$',block),
                'Wrong denominator in Appendix B caption')
        for line in block.splitlines():
            if not re.match(r'^\s*\d+\s*&',line):
                continue
            fields = without_row_end(line).split('&')
            require(len(fields) == 4, 'Wrong direction columns')
            i = int(fields[0]); require(i not in directions,'Duplicate six-point direction')
            vector = re.search(r'\$\(([-\d,\s]+)\)\$', fields[3])
            require(vector is not None, 'Malformed printed six-point vector')
            directions[i] = (fields[1].strip(),int(fields[2]),
                              [Q(int(x),den) for x in vector[1].split(',')])
        require(set(directions) == set(range(1,len(kernels)+1)), 'Missing six-point directions')
        for i,term in enumerate(kernels,1):
            require(directions[i] == (term['anchor_bits'],term['degree'],list(map(Q,term['vector']))),
                    'Printed six-point direction fails at '+str(v)+', row '+str(i))
        result[str(v)] = {'integer_multipliers':24,'directions':len(kernels),
                          'common_scale_checked':True,'raw_identity_verified':True}
    return result


def pure405(main_source, supplement_source):
    source = main_source + '\n' + supplement_source
    # The complete data belong in the supplement, exactly once across both inputs.
    for label in ('tab:405-countdata', 'tab:405-kerneldata'):
        table(supplement_source, label)
        table(source, label)
    root = HERE.parent/'data/pure405'
    record = json.loads((root/'data/degree2/pure_exact_certificate.json').read_text())
    graphs = [tuple(map(int,line.split())) for line in (root/'atlas/graphs.txt').read_text().splitlines()]
    require([g[0] for g in graphs] == list(range(1252)), 'Atlas indices changed')
    terms = record['terms']
    require(len(terms) == 174, 'Unexpected 405 source count')
    printed = {}
    for line in table(source,'tab:405-countdata').splitlines():
        if not re.match(r'^\s*\d+\s*&',line):
            continue
        fields = without_row_end(line).split('&')
        require(len(fields) == 6, 'Wrong number of 405 count-table columns')
        for fields in (fields[:3],fields[3:]):
            if not fields[0].strip():
                require(not fields[1].strip() and not fields[2].strip(), 'Unindexed 405 count source')
                continue
            i = int(fields[0]); require(i not in printed,'Duplicate 405 source index')
            graph = re.fullmatch(r'\s*\$\((\d+),(\d+)\)\$\s*',fields[1])
            require(graph is not None,'Invalid 405 graph pair')
            printed[i] = {'kind':'nonnegative_count','graph':tuple(map(int,graph.groups())),
                          'weight':fraction(fields[2])}
    require(set(printed) == set(range(1,68)), 'Missing printed 405 count source')
    for line in table(source,'tab:405-kerneldata').splitlines():
        if not re.match(r'^\s*\d+\s*&',line):
            continue
        fields = without_row_end(line).split('&',2)
        require(len(fields) == 3, 'Wrong 405 kernel-table columns')
        i = int(fields[0]); require(i not in printed,'Duplicate 405 source index')
        graph = re.fullmatch(r'\s*\$\((\d+),(\d+);(\d+)\)\$\s*',fields[1])
        require(graph is not None,'Invalid printed 405 anchor/degree')
        m,b,ell = map(int,graph.groups())
        vector = array_values(fields[2])
        rationals = re.findall(r'\\(?:frac|dfrac|tfrac)\{(-?\d+)\}\{(\d+)\}',fields[2])
        require(len(rationals)==1,'Expected one exact rational weight in 405 source row')
        weight = Q(int(rationals[0][0]),int(rationals[0][1]))
        require(len(vector) == 2**m,'Incomplete printed 405 profile vector')
        printed[i] = {'kind':'PSD','graph':(m,b),'degree':ell,
                      'weight':weight,'direction':vector}
    require(set(printed) == set(range(1,175)), 'Missing printed 405 sources')
    for i,term in enumerate(terms,1):
        source_data = term['source']; row = printed[i]
        require(row['kind'] == source_data['kind'], '405 source-kind mismatch '+str(i))
        gid = source_data['graph_id'] if row['kind']=='nonnegative_count' else source_data['anchor']
        require(row['graph'] == graphs[gid][1:], 'Printed 405 graph/anchor mismatch '+str(i))
        require(row['weight'] == Q(term['weight'])/Q(term['scale']) and row['weight']>0,
                'Printed 405 rational coefficient mismatch '+str(i))
        if row['kind']=='PSD':
            require(row['degree'] == source_data['degree'] and row['direction'] == source_data['direction'],
                    'Printed 405 direction/degree mismatch '+str(i))

    # Check the illustrative degree-two 76th source in the main proof as well.
    marker = 'For example, source $i=76$'
    require(source.count(marker) <= 1,'Duplicated main 405 worked example')
    if marker in source:
        require(marker in main_source, 'Worked 405 example must be in the main manuscript')
        example = main_source.split(marker,1)[1].split('This one source',1)[0]
        require(array_values(example) == printed[76]['direction'],'Main-proof f_76 differs from appendix')
        example_weight = re.search(r'w_\{76\}\s*=\s*(\\frac\{\d+\}\{\d+\})',example)
        require(example_weight is not None and fraction(example_weight[1]) == printed[76]['weight'],
                'Main-proof w_76 differs from appendix')
    require(printed[76]['graph'] == (5,222) and printed[76]['degree']==2,'K23 example mismatch')

    # Recombine the exact independently reconstructed raw rows with PRINTED weights.
    raw_recipes = json.loads((root/'data/degree2/core_recipes_405.json').read_text())
    ids = raw_recipes['cores']
    recipes = {int(i):{int(j):Q(q) for j,q in row.items()}
               for i,row in raw_recipes['expressions'].items()}
    require(len(ids)==174 and all(recipes[gid]=={j:Q(1)} for j,gid in enumerate(ids)),
            '405 free rows must be an identity')
    total = [Q(0)]*174
    for i,term in enumerate(terms,1):
        reduced = [Q(0)]*174
        for gid,q in term['full_row'].items():
            for j,a in recipes[int(gid)].items():
                reduced[j] += Q(q)*a
        for j,q in enumerate(reduced):
            total[j] += printed[i]['weight']*q
    require(total == list(map(Q,record['combined_core_coefficients'])),
            'Printed 405 combination does not reconstruct saved exact coefficients')
    bounds = [2**(graphs[gid][1]*(graphs[gid][1]-1)//2) for gid in ids]
    require(bounds == record['core_density_upper_bounds'],'405 count upper bounds differ')
    c0 = total[0]
    error = -sum((q*b for q,b in zip(total[1:],bounds[1:]) if q<0),Q(0))
    margin = c0-error
    replay = json.loads((root/'data/degree2/pure_exact_replay.json').read_text())
    require(replay['status']=='PASS' and margin==Q(replay['exact_margin'])==Q(record['exact_margin']),
            '405 margin differs from independent exact replay')
    require(c0>Q(570,10**10) and 0<=error<Q(11,10**10) and margin>Q(559,10**10),
            'A printed strict rational bound is false')

    for name, document in (('manuscript', main_source), ('405 supplement', supplement_source)):
        c0_bounds = re.findall(r'c_0\s*>\s*\\frac\{(\d+)\}\{10\^\{(\d+)\}\}', document)
        e_bounds = re.findall(r'0\s*\\le\s*E(?:_\{405\})?\s*<\s*\\frac\{(\d+)\}\{10\^\{(\d+)\}\}', document)
        require(len(c0_bounds)==1 and Q(int(c0_bounds[0][0]),10**int(c0_bounds[0][1]))==Q(570,10**10),
                'Missing, duplicated, or changed printed c0 lower bound in ' + name)
        require(len(e_bounds)==1 and Q(int(e_bounds[0][0]),10**int(e_bounds[0][1]))==Q(11,10**10),
                'Missing, duplicated, or changed printed E upper bound in ' + name)
        require(len(re.findall(r'>\s*\\frac\{559\}\{10\^\{10\}\}', document))==1,
                'Expected exactly one checked lower bound on the margin in ' + name)
    return {'source_coefficients':174,'count_sources':67,'full_directions':107,
            'source_76_checked':True,'main_example_76_present':marker in source,'exact_coefficients_recombined':True,
            'exact_replay_margin_matched':True,'c0_gt':'570/10^10','E_lt':'11/10^10',
            'margin_gt':'559/10^10'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tex',type=Path,default=HERE.parent/'paper/srg_counting_gegenbauer_article_1014.tex')
    parser.add_argument('--supplement',type=Path,
                        help='405 supplement TeX; defaults to the named sibling of --tex')
    args=parser.parse_args()
    manuscript = args.tex.resolve()
    supplement = (args.supplement if args.supplement is not None else
                  manuscript.with_name('srg_counting_gegenbauer_article_1014_405_supplement.tex')).resolve()
    require(manuscript.is_file(), 'Required manuscript TeX does not exist: ' + str(manuscript))
    require(supplement.is_file(), 'Required 405 supplement TeX does not exist: ' + str(supplement))
    require(manuscript != supplement, 'Manuscript and supplement must be different files')
    source=manuscript.read_text(encoding='utf-8')
    supplement_source=supplement.read_text(encoding='utf-8')
    print('Manuscript TeX: ' + str(manuscript), flush=True)
    print('405 supplement TeX: ' + str(supplement), flush=True)
    six=sixpoint(source)
    print('PASS actual TeX: 48 integer multipliers, 2 common scales, and 33 six-point directions.',flush=True)
    pure=pure405(source,supplement_source)
    print('PASS actual TeX: 174 rational source weights, all 107 full 405 directions, and source 76.',flush=True)
    print('PASS exact replay: c0 > 570/10^10; 0 <= E < 11/10^10; margin > 559/10^10.',flush=True)
    print(json.dumps({'verified':True,'manuscript':str(manuscript),'supplement':str(supplement),
                      'six_point':six,'405':pure,
                      'arithmetic':'fractions.Fraction; standard library only'},indent=2))

if __name__=='__main__':
    main()
